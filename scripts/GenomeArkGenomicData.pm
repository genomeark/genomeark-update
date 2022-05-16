package GenomeArkGenomicData;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(loadSummaries writeSummaries estimateRawDataScaling);

use strict;
use warnings;

use GenomeArkUtility;

#use Time::Local;
#use List::Util;
use File::Basename;

my $aws          = "aws --no-sign-request";
my $seqrequester = "./seqrequester/build/bin/seqrequester";
my $bamsummarize = "./seqrequester/build/bin/bamsummarize";

my $downloadRAW = 1;

my %dataBytes;   #  Number of bytes in a genomic_data sequence file.
my %dataQuant;   #  Amount of raw file summarized (in gigabytes, or "all").
my %dataBases;   #  Number of bases in piece summarized.
my %dataReads;   #  Number of reads in piece summarized.

#
#  Load a single data summary, which MUST exist.  If the summary is for a
#  partial file, scale the summarized values to approximate the whole file.
#
sub loadSummaryFromFile ($$$$) {
    my $size    = shift @_;
    my $file    = shift @_;
    my $amount  = shift @_;
    my $summary = shift @_;

    $dataBytes{$file} = $size;
    $dataQuant{$file} = $amount;
    $dataBases{$file} = 0;
    $dataReads{$file} = 0;

    print STDERR "  Load summary from '$summary'\n";

    open(ST, "< $summary") or die;
    while (<ST>) {
        $dataBases{$file} = $1   if (m/^G=(\d+)\s/);
        $dataReads{$file} = $1   if (m/^001.000x\s+(\d+)\s/)
    }
    close(ST);

    if ($amount =~ m/(\d+)GB/) {
        my $scale = $size / ($1 * 1024.0 * 1024.0 * 1024.0);

        $dataBases{$file} = int($dataBases{$file} * $scale);
        $dataReads{$file} = int($dataReads{$file} * $scale);
    }
    if ($amount =~ m/(\d+)MB/) {
        my $scale = $size / ($1 * 1024.0 * 1024.0);

        $dataBases{$file} = int($dataBases{$file} * $scale);
        $dataReads{$file} = int($dataReads{$file} * $scale);
    }

    if ($dataBases{$file} == 0) {    #  If an invalid estimate, delete the
        unlink $summary;             #  summary file and clear our data.

        delete $dataBytes{$file};
        delete $dataQuant{$file};
        delete $dataBases{$file};
        delete $dataReads{$file};
    }
}

#
#  Load a summary, if it exists.
#
sub loadSummary (@) {
    my $size = shift @_;
    my $file = shift @_;

    #  If the summary exists, and it is for all the data, do nothing.
    if (exists($dataQuant{$file}) && $dataQuant{$file} eq "all") {
    }

    #  If the 'all' summary exists, load it.
    elsif (-e "$file.summary") {
        loadSummaryFromFile($size, $file, "all", "$file.summary");
    }

    #  If any subset summaries exist, load them.
    elsif (-e "$file.1GB.summary") {
        loadSummaryFromFile($size, $file, "1GB", "$file.1GB.summary");
    }
    elsif (-e "$file.9MB.summary") {
        loadSummaryFromFile($size, $file, "9MB", "$file.9MB.summary");
    }
}


#
#  Load existing data summaries from 'species/$name/genomic_data.summary',
#  then scan the list of seqFiles to find any new summaries that aren't in
#  the global file.
#
sub loadSummaries ($$) {
    my $name     = shift @_;
    my $seqFiles = shift @_;

    undef %dataBytes;
    undef %dataQuant;
    undef %dataBases;
    undef %dataReads;

    if (-e "species/$name/genomic_data.summary") {
        print STDERR "  Load summaries from species/$name/genomic_data.summary\n";

        open(ST, "< species/$name/genomic_data.summary") or die;
        while (<ST>) {
            s/^\s+//;
            s/\s+$//;

            next   if m/bytes/;
            next   if m/-----/;

            my ($bytes, $date, $quant, $bases, $reads, $file) = split '\s+', $_;

            if ($bases > 0) {
                $dataBytes{$file} = $bytes;
                $dataQuant{$file} = $quant;
                $dataBases{$file} = $bases;
                $dataReads{$file} = $reads;
            }
        }
        close(ST);
    }

    foreach my $type (sort keys %$seqFiles) {
        my @files = split '\0', $$seqFiles{$type};

        foreach my $sizefile (@files) {
            loadSummary(split '\s+', $sizefile);
        }
    }
}



sub writeSummaries ($) {
    my $name  = shift @_;

    system("mkdir -p species/$name")   if (! -e "species/$name");

    open(ST, "> species/$name/genomic_data.summary") or die "Failed to open 'species/$name/genomic_data.summary' for writing: $!\n";

    print ST "          bytes            date        quantity           bases           reads  file\n";
    print ST "--------------- --------------- --------------- --------------- ---------------  ----------\n";

    foreach my $file (sort keys %dataBytes) {
        my  $date = 0;

        next   if (! exists($dataBytes{$file}));
        next   if (! exists($dataQuant{$file}));
        next   if (! exists($dataBases{$file}));
        next   if (! exists($dataReads{$file}));

        printf(ST "%15d %15d %15s %15d %15d  %s\n", $dataBytes{$file}, $date, $dataQuant{$file}, $dataBases{$file}, $dataReads{$file}, $file);
    }
    close(ST);
}



#
#  If the summary doesn't exist for all the data, download the whole file and summarize.
#
sub downloadFullAndSummarize ($$) {
    my $size  = shift @_;
    my $file  = shift @_;
    my $fdir  = dirname($file);

    #  If the whole summary exists, do not recompute it.

    return   if (exists($dataQuant{$file}) && $dataQuant{$file} eq "all");

    die      if (-e "$file.summary");   #  It should be loaded already!

    #  Fetch the data if it doesn't exist already.

    if ($downloadRAW == 0) {
        printf "FETCH data - size %6.3f GB -- DISABLED\n", $size / 1024 / 1024 / 1024;
        printf "    s3://genomeark/$file\n";
        printf "\n";
        return;
    }

    system("mkdir -p downloads/$fdir");   #  Make a place to download the file
    system("mkdir -p           $fdir");   #  and a place to write the summary.

    if (! -e "downloads/$file") {
        printf "FETCH data - size %6.3f GB\n", $size / 1024 / 1024 / 1024;
        printf "    s3://genomeark/$file \\\n";
        printf "    ->   downloads/$file\n";
        printf "\n";

        system("$aws s3 cp s3://genomeark/$file downloads/$file");
    }

    #  If the download failed....do what?

    if (! -e "downloads/$file") {
        printf "  FAILED.\n";
        return;
    }

    #  If a bam, convert to fastq then summarize, otherwise, summarize directly.

    if ($file =~ m/bam$/) {
        printf "EXTRACT and SUMMARIZE $file.summary\n";

        system("samtools fasta downloads/$file | $seqrequester summarize - > $file.summary");
        #system("gzip -dc downloads/$file | $bamsummarize > $file.summary");
    }
    else {
        printf "SUMMARIZE $file.summary\n";

        system("$seqrequester summarize downloads/$file > $file.summary");
    }

    #  Parse the summary to find the number of bases in the dataset.

    loadSummary($size, $file);
}


#
#  If no estimated number of bases exists, download a part of the file and summarize.
#
sub downloadPartAndSummarize ($$) {
    my $size   = shift @_;
    my $file   = shift @_;
    my $fdir   = dirname($file);

    #my $gigab  = 1;
    #my $gigau  = "GB";
    #my $bytes  = $gigab * 1024 * 1024 * 1024;

    my $gigab  = 9;
    my $gigau  = "MB";
    my $bytes  = $gigab * 1024 * 1024;

    #  If some kind of estimate exists, we'll use that.

    return   if (exists($dataBases{$file}));

    #  If the download exists, summarize it and return the number of bases,

    if (-e "downloads/$file") {
        downloadAndSummarize($size, $file);
        return;
    }

    #  Otherwise, fetch a tiny bit of the data and pass it directly to the summarizer.

    if ($downloadRAW == 0) {
        printf "\n";
        printf "FETCH data - size %6.3f GB -- DISABLED\n", $size / 1024 / 1024 / 1024;
        printf "    s3://genomeark/$file\n";
        printf "\n";
        return;
    }

    system("mkdir -p downloads/$fdir");   #  Make a place to download the file
    system("mkdir -p           $fdir");   #  and a place to write the summary.

    #  While 's3 cp' is happy writing to stdout, 's3api get-object' is not,
    #  and requires a file name.  Rather annoying, since this is where we
    #  could actually get away with piping it to the summarizer.

    if ($file =~ m/bam$/) {
        printf "\n";
        printf "FETCH data - $gigab $gigau out of total size %6.3f GB\n", $size / 1024 / 1024 / 1024;
        printf "    s3://genomeark/$file\n";
        printf "    bamsummarize > $file.summary\n";

        #print("$aws s3api get-object --bucket genomeark --key $file --range bytes=0-$bytes - | samtools fasta | $seqrequester summarize - > $file.${gigab}${gigau}.summary\n");

        my $df = "downloads/$file.${gigab}${gigau}.bam";

        system("$aws s3api get-object --bucket genomeark --key $file --range bytes=0-$bytes $df > $df.err 2>&1")  if (! -e "$df");
        system("samtools fasta $df 2> $df.samtools.err | $seqrequester summarize - > $file.${gigab}${gigau}.summary 2> $df.seqrequester.err");
        #system("rm -f $df");
        printf "\n";
    }
    else {
        printf "\n";
        printf "FETCH data - $gigab $gigau out of total size %6.3f GB\n", $size / 1024 / 1024 / 1024;
        printf "    s3://genomeark/$file\n";
        printf "    seqrequester > $file.summary\n";

        my $df = "downloads/$file.${gigab}${gigau}.gz";

        system("$aws s3api get-object --bucket genomeark --key $file --range bytes=0-$bytes $df > $df.err 2>&1")  if (! -e "$df");
        system("gzip -dc $df 2> $df.gzip.err | $seqrequester summarize - > $file.${gigab}${gigau}.summary 2> $df.seqrequester.err");
        #system("rm -f $df");
        printf "\n";
    }

    loadSummary($size, $file);
}



#
#  Given a species_name, a genomic data type and a list of seqFiles
#  ("filesize datafile" separated by \0 bytes), compute a
#  scalaing factor that will convert a filesize into an estimated number of
#  bases in the (compressed) file.
#
sub estimateRawDataScaling ($$$) {
    my $data     =             shift @_;
    my $name     = $$data{"name_"};
    my $type     =             shift @_;
    my $files    =             shift @_;
    my $scaling;

    return   if (!defined($files) || ($files eq ""));

    my @files    = split '\0', $files;
    my $filesLen = scalar(@files);

    return   if ($filesLen == 0);

    @files = sort { 
        my ($s1, $n1) = split '\s+', $a;
        my ($s2, $n2) = split '\s+', $b;
        return($s1 <=> $s2);
    } @files;

    my $f1 = int(1 * $filesLen / 4);   #  Pick three files representative file (indices)
    my $f2 = int(2 * $filesLen / 4);   #  from the list of input files sorted by size.
    my $f3 = int(3 * $filesLen / 4);

    my ($size1, $file1) = split '\s', $files[$f1];   #  Extract name and size from input list;
    my ($size2, $file2) = split '\s', $files[$f2];   #  bases and seqs are undefined.
    my ($size3, $file3) = split '\s', $files[$f3];

    #  Print a list of the data present, marking which ones we have
    #  summarized already or will use for estimation.

    print " - $type:\n";

    for (my $ii=0; $ii<$filesLen; $ii++) {
        my ($size, $file) = split '\s+', $files[$ii];

        my $se = (exists($dataBases{$file})) ? "S" : " ";    #  Summary exists?
        my $us = (($ii == $f1) ||                            #  Used for estimation?
                  ($ii == $f2) ||
                  ($ii == $f3)) ? "E" : " ";

        printf "    %s%s %12d %s\n", $se, $us, $size, $file;
    }

    #  If needed, download a subset of the data from three files and use that to
    #  estimate the total number of bases in the file.

    downloadPartAndSummarize($size1, $file1);    writeSummaries($name);
    downloadPartAndSummarize($size2, $file2);    writeSummaries($name);
    downloadPartAndSummarize($size3, $file3);    writeSummaries($name);

    #  Get an estimate of (or, rarely, the actual) number of bases in each file.

    if (!exists($dataBases{$file1}) ||        #  Fail if any of the size computations failed.
        !exists($dataBases{$file2}) ||
        !exists($dataBases{$file3})) {
        print "FAILED TO ESTIMATE SIZES FOR:\n";
        print "  File 1 $file1\n"  if (!exists($dataBases{$file1}));
        print "  File 2 $file2\n"  if (!exists($dataBases{$file2}));
        print "  File 3 $file3\n"  if (!exists($dataBases{$file3}));
        exit(1);
    }

    $scaling = (($dataBases{$file1} + $dataBases{$file2} + $dataBases{$file3}) /  #  Compute scaling, then limit
                ($size1             + $size2             + $size3));              #  precision to prevent
    $scaling = int($scaling * 10000) / 10000;                                     #  churn from bad math.

    $$data{"data_${type}_scale"} = $scaling;
}

return(1);
