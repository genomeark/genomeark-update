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


#my $seqrequester = "./seqrequester/build/bin/seqrequester";
my $seqrequester = "seqrequester";
my $bamsummarize = "bamsummarize";

my $downloadRAW = 0;

my %dataBytes;   #  Number of bytes in a genomic_data sequence file
my %dataBases;   #  Number of bases in ....
my %dataReads;   #  Number of reads in ....


#
#  Load existing data summaries from 'species/$name/genomic_data.summary'.
#
#  Then scan seqFiles and add or update summary information for any missing
#  or updated files.
#

sub loadSummary (@) {
    my $size = shift @_;
    my $file = shift @_;

    unlink "$file.summary"   if (-z "$file.summary");   #  Remove incomplete summaries.

    return   if (exists($dataBases{$file}));
    return   if (! -e "$file.summary");

    $dataBytes{$file} = $size;
    $dataBases{$file} = 0;
    $dataReads{$file} = 0;

    open(ST, "< $file.summary");
    while (<ST>) {
        $dataBases{$file} = $1   if (m/^G=(\d+)\s/);
        $dataReads{$file} = $1   if (m/^001.000x\s+(\d+)\s/)
    }
    close(ST);

    return($dataBases{$file});
}

sub loadSummaries ($$) {
    my $name     = shift @_;
    my $seqFiles = shift @_;

    undef %dataBytes;
    undef %dataBases;
    undef %dataReads;

    if (-e "species/$name/genomic_data.summary") {
        open(ST, "< species/$name/genomic_data.summary");
        while (<ST>) {
            s/^\s+//;
            s/\s+$//;

            next   if m/bytes/;
            next   if m/-----/;

            my ($bytes, $date, $bases, $reads, $file) = split '\s+', $_;

            $dataBytes{$file} = $bytes;
            $dataBases{$file} = $bases;
            $dataReads{$file} = $reads;
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

    print ST "          bytes            date           bases           reads  file\n";
    print ST "--------------- --------------- --------------- ---------------  ----------\n";

    foreach my $file (sort keys %dataBytes) {
        my  $date = 0;

        next   if (! exists($dataBytes{$file}));
        next   if (! exists($dataBases{$file}));
        next   if (! exists($dataReads{$file}));

        printf(ST "%15d %15d %15d %15d  %s\n", $dataBytes{$file}, $date, $dataBases{$file}, $dataReads{$file}, $file);
    }
    close(ST);
}



sub downloadAndSummarize ($$$) {
    my $name  = shift @_;
    my $ndir  = dirname($name);
    my $size  = shift @_;
    my $file  = shift @_;
    my $bases = 0;

    if (exists($dataBases{$name})) {   #  Summary exists, just return the
        return($dataBases{$name});     #  number of bases in the file.
    }

    #  Fetch the data if it doesn't exist already.

    if ($downloadRAW == 0) {
        printf "FETCH data - size %6.3f GB -- DISABLED\n", $size / 1024 / 1024 / 1024;
        printf "  aws --no-progress --no-sign-request s3 cp \\\n";
        printf "    s3://genomeark/$name\n";
        printf "\n";

        return($size / 2);
    }

    if (! -e "downloads/$name") {
        printf "FETCH data - size %6.3f GB\n", $size / 1024 / 1024 / 1024;
        printf "  aws --no-progress --no-sign-request s3 cp \\\n";
        printf "    s3://genomeark/$name \\\n";
        printf "         downloads/$name\n";
        printf "\n";

        system("mkdir -p downloads")   if (! -e "downloads");
        system("aws --no-progress --no-sign-request s3 cp s3://genomeark/$name downloads/$name")  if ($downloadRAW);
    }

    #  If the download failed....do what?

    if (! -e "downloads/$name") {
        return(0);
    }

    #  If a bam, convert to fastq then summarize, otherwise, summarize directly.

    if ($name =~ m/bam$/) {
        printf "EXTRACT and SUMMARIZE $name.summary\n";

        system("mkdir -p $ndir")   if (! -e $ndir);
        #system("samtools fasta downloads/$name | $seqrequester summarize - > $name.summary");
        system("gzip -dc downloads/$name | $bamsummarize > $name.summary");
    }
    else {
        printf "SUMMARIZE $name.summary\n";

        system("mkdir -p $ndir")   if (! -e $ndir);
        system("$seqrequester summarize downloads/$name > $name.summary");
    }

    #  Parse the summary to find the number of bases in the dataset.

    return(loadSummary($size, $file));
}



sub downloadPipeSummarize ($$$) {
    my $name  = shift @_;
    my $ndir  = dirname($name);
    my $size  = shift @_;
    my $file  = shift @_;

    #  Summary exists, just return the number of bases in the file.

    if (exists($dataBases{$name})) {
        return($dataBases{$name});
    }

    #  If the download exists, summarize it and return the number of bases,

    if (-e "downloads/$name") {
        return(downloadAndSummarize($name, $size, $file));
    }

    #  Otherwise, fetch the data directly into the summarizer.

    if ($downloadRAW == 0) {
        printf "\n";
        printf "FETCH data - size %6.3f GB -- DISABLED\n", $size / 1024 / 1024 / 1024;
        printf "  aws --no-progress --no-sign-request s3 cp \\\n";
        printf "    s3://genomeark/$name\n";
        printf "\n";

        return($size / 2);
    }
    elsif ($name =~ m/bam$/) {
        printf "\n";
        printf "FETCH data - size %6.3f GB\n", $size / 1024 / 1024 / 1024;
        printf "  aws --no-progress --no-sign-request s3 cp \\\n";
        printf "    s3://genomeark/$name | bamsummarize > $name.summary\n";

        #system("mkdir -p $ndir")   if (! -e $ndir);
        #system("aws --no-progress --no-sign-request s3 cp s3://genomeark/$name - | gzip -dc | $bamsummarize > $name.summary");
        system("mkdir -p $ndir")   if (! -e $ndir);
        system("aws --no-progress --no-sign-request s3 cp s3://genomeark/$name - | samtools fasta - | $seqrequester summarize - > $name.summary");
        printf "\n";
    }
    else {
        printf "\n";
        printf "FETCH data - size %6.3f GB\n", $size / 1024 / 1024 / 1024;
        printf "  aws --no-progress --no-sign-request s3 cp \\\n";
        printf "    s3://genomeark/$name | seqrequester > $name.summary\n";

        system("mkdir -p $ndir")   if (! -e $ndir);
        system("aws --no-progress --no-sign-request s3 cp s3://genomeark/$name - | gzip -dc | $seqrequester summarize - > $name.summary");
        printf "\n";
    }

    return(loadSummary($size, $file));
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

    my ($size1, $name1, $bases1, $seqs1) = split '\s', $files[$f1];   #  Extract name and size from input list;
    my ($size2, $name2, $bases2, $seqs2) = split '\s', $files[$f2];   #  bases and seqs are undefined.
    my ($size3, $name3, $bases3, $seqs3) = split '\s', $files[$f3];

    #  Print a list of the data present, marking which ones we have
    #  summarized already or will use for estimation.

    print " - $type:\n";

    for (my $ii=0; $ii<$filesLen; $ii++) {
        my ($size, $name) = split '\s+', $files[$ii];

        my $se = (exists($dataBases{$name})) ? "S" : " ";    #  Summary exists?
        my $us = (($ii == $f1) ||                            #  Used for estimation?
                  ($ii == $f2) ||
                  ($ii == $f3)) ? "E" : " ";

        printf "    %s%s %12d %s\n", $se, $us, $size, $name;
    }

    #  Retrieve the number of bases in each of the three files, fetching and
    #  summarizing if needed.

    $bases1 = downloadPipeSummarize($name1, $size1, $f1);
    $bases2 = downloadPipeSummarize($name2, $size2, $f2);
    $bases3 = downloadPipeSummarize($name3, $size3, $f3);

    if (($bases1 == 0) ||        #  Fail if any of the size computations failed.
        ($bases2 == 0) ||
        ($bases3 == 0)) {
        print "FAILED TO ESTIMATE SIZES.\n";
        print "  1 - $bases1 - $name1\n";
        print "  2 - $bases2 - $name2\n";
        print "  3 - $bases3 - $name3\n";
        exit(1);
    }

    $scaling = ($bases1 + $bases2 + $bases3) / ($size1 + $size2 + $size3);   #  Compute scaling, limit precision
    $scaling = int($scaling * 10000) / 10000;                                #  to prevent churn from bad math.

    $$data{"data_${type}_scale"} = $scaling;
}

return(1);
