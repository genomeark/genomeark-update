package GenomeArkGenomicData;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(loadSummaries writeSummaries estimateRawDataScaling);

use strict;
use warnings;

use GenomeArkUtility;

#use Time::Local;
use List::Util qw(min max);
use File::Basename;

my $aws          = "aws --no-sign-request";
my $seqrequester = "./seqrequester/build/bin/seqrequester";
my $bamsummarize = "./seqrequester/build/bin/bamsummarize";

my %dataBytes;   #  Number of bytes in a genomic_data sequence file.
my %dataQuant;   #  Amount of raw file summarized (in gigabytes, or "all").
my %dataBases;   #  Number of bases in piece summarized.
my %dataReads;   #  Number of reads in piece summarized.



sub makeSummaryFileName ($$) {
    my $n = shift @_;
    my $s = shift @_;
    my $N = dirname($n) . "/$s-bytes--" . basename($n) . ".summary";

    return($N);
}


#
#  Load a single data summary, which MUST exist.  If the summary is for a
#  partial file, scale the summarized values to approximate the whole file.
#
#  If, after all that work, we end up with no estimate, forget everything.
#
sub loadSummaryFromFile ($$$$) {
    my $filesize  = shift @_;
    my $filename  = shift @_;
    my $sumsize   = shift @_;
    my $summary   = shift @_;

    $dataBytes{$filename} = $filesize;
    $dataQuant{$filename} = $sumsize;
    $dataBases{$filename} = 0;
    $dataReads{$filename} = 0;

    print STDERR "    Import summary from '$summary'\n";

    open(ST, "< $summary") or die;
    while (<ST>) {
        $dataBases{$filename} = $1   if (m/^G=(\d+)\s/);
        $dataReads{$filename} = $1   if (m/^001.000x\s+(\d+)\s/)
    }
    close(ST);

    if ($filesize != $sumsize) {
        $dataBases{$filename} = int($dataBases{$filename} * $filesize / $sumsize);
        $dataReads{$filename} = int($dataReads{$filename} * $filesize / $sumsize);
    }

    if ($dataBases{$filename} == 0) {
        unlink $summary;

        delete $dataBytes{$filename};
        delete $dataQuant{$filename};
        delete $dataBases{$filename};
        delete $dataReads{$filename};
    }
}

#
#  Load a summary, if it exists.
#
sub loadSummary (@) {
    my $filesize = shift @_;
    my $filename = shift @_;

    my $newsize = 0;
    my $newfile = "";

    my $sumnames = makeSummaryFileName($filename, "*");

    #  Pick the summary of the largest amount of data.
    open(LS, "ls $sumnames 2>/dev/null |");
    while (<LS>) {
        chomp;
        if ($_ =~ m!/(\d+)-bytes--.*\.summary$!) {
            if ($1 > $newsize) {
                $newsize = $1;
                $newfile = $_;
            }
        }
    }
    close(LS);

    #  If a full summary exists, always use that.
    if (-e "${filename}.summary") {
        $newsize = $filesize;
        $newfile = "${filename}.summary";
    }

    #  Load the new summary if it summarized more data than we did.

    my $sumsize  = (exists($dataQuant{$filename})) ? $dataQuant{$filename} : 0;

    if ($newsize > $sumsize) {
        loadSummaryFromFile($filesize, $filename, $newsize, $newfile);
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

            $quant = 9 * 1024 * 1024         if ($quant eq "9MB");
            $quant = 1 * 1024 * 1024 * 1024  if ($quant eq "1GB");
            $quant = 2 * 1024 * 1024 * 1024  if ($quant eq "2GB");
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
sub downloadFullAndSummarize ($$$) {
    my $filesize = shift @_;
    my $filename= shift @_;
    my $fdir     = dirname($filename);
    my $download = shift @_;

    #  If the whole summary exists, do not recompute it.

    return   if (exists($dataQuant{$filename}) && $dataQuant{$filename} eq "all");

    die      if (-e "$filename.summary");   #  It should be loaded already!

    #  Fetch the data if it doesn't exist already.

    if ($download == 0) {
        printf "FETCH data - size %6.3f GB -- DISABLED\n", $filesize / 1024 / 1024 / 1024;
        printf "    s3://genomeark/$filename\n";
        printf "\n";
        return;
    }

    system("mkdir -p downloads/$fdir");   #  Make a place to download the file
    system("mkdir -p           $fdir");   #  and a place to write the summary.

    if (! -e "downloads/$filename") {
        printf "FETCH data - size %6.3f GB\n", $filesize / 1024 / 1024 / 1024;
        printf "    s3://genomeark/$filename \\\n";
        printf "    ->   downloads/$filename\n";
        printf "\n";

        system("$aws s3 cp s3://genomeark/$filename downloads/$filename");
    }

    #  If the download failed....do what?

    if (! -e "downloads/$filename") {
        printf "  FAILED.\n";
        return;
    }

    #  If a bam, convert to fastq then summarize, otherwise, summarize directly.

    if ($filename =~ m/bam$/) {
        printf "EXTRACT and SUMMARIZE $filename.summary\n";

        system("samtools fasta downloads/$filename | $seqrequester summarize - > $filename.summary");
        #system("gzip -dc downloads/$filename | $bamsummarize > $filename.summary");
    }
    else {
        printf "SUMMARIZE $filename.summary\n";

        system("$seqrequester summarize downloads/$filename > $filename.summary");
    }

    #  Parse the summary to find the number of bases in the dataset.

    loadSummary($filesize, $filename);
}


#
#  If no estimated number of bases exists, download a part of the file and summarize.
#
sub downloadPartAndSummarize ($$$) {
    my $filesize = shift @_;
    my $filename = shift @_;
    my $fdir     = dirname($filename);
    my $download = shift @_;

    #  How much did we summarize already, and how much do we want to summarize?

    my $sumsize = (exists($dataQuant{$filename})) ? $dataQuant{$filename} : 0;
    my $newsize = 10 * 1024 * 1024;
    my $sumname = makeSummaryFileName($filename, $newsize);


    #  If we've summarized all the data, we can't improve it.
    if    ($sumsize >= $filesize) {
        return;
    }

    #  If there is ALREADY a data file around to summarize, make a full summary.
    elsif (-e "downloads/$filename") {
        downloadFullAndSummarize($filesize, $filename, $download);
        return;
    }

    #  If we've summarized more than requested, don't do more work to get less!
    elsif ($sumsize >= $newsize) {
        return;
    }

    #
    #  Yay!  Let's download and summarize something!
    #

    if ($download == 0) {
        printf "\n";
        printf "  SKIP  - %6.3f GB out of %6.3f GB - s3://genomeark/$filename\n", 0.0, $filesize / 1024 / 1024 / 1024;
        return;
    }

    #  If the file isn't hugely bigger than the size we'd sample, just do the whole thing.
    if ($filesize < 2 * $newsize) {
        downloadFullAndSummarize($filesize, $filename, $download);
        return;
    }

    #  While 's3 cp' is happy writing to stdout, 's3api get-object' is not,
    #  and requires a file name.  Rather annoying, since this is where we
    #  could actually get away with piping it to the summarizer.

    system("mkdir -p downloads/$fdir");   #  Make a place to download the file
    system("mkdir -p           $fdir");   #  and a place to write the summary.

    if ($filename =~ m/bam$/) {
        my $df = "downloads/" . dirname($filename) . "/${newsize}-bytes--" . basename($filename);

        if (! -e "$df") {
            printf "  FETCH     - %6.3f out of %6.3f GB - s3://genomeark/$filename\n", $newsize / 1024 / 1024 / 1024, $filesize / 1024 / 1024 / 1024;
            system("$aws s3api get-object --bucket genomeark --key $filename --range bytes=0-$newsize $df > $df.err 2>&1")  if (! -e "$df");
        } else {
            printf "  CACHED    - %6.3f out of %6.3f GB - s3://genomeark/$filename\n", $newsize / 1024 / 1024 / 1024, $filesize / 1024 / 1024 / 1024;
        }

        printf "  SUMMARIZE - samtools | seqrequester > $sumname\n";
        system("samtools fasta $df 2> $df.samtools.err | $seqrequester summarize - > $sumname 2> $df.seqrequester.err");
        #system("rm -f $df");
    }
    else {
        my $df = "downloads/" . dirname($filename) . "/${newsize}-bytes--" . basename($filename);

        if (! -e "$df") {
            printf "  FETCH     - %6.3f out of %6.3f GB - s3://genomeark/$filename\n", $newsize / 1024 / 1024 / 1024, $filesize / 1024 / 1024 / 1024;
            system("$aws s3api get-object --bucket genomeark --key $filename --range bytes=0-$newsize $df > $df.err 2>&1")  if (! -e "$df");
        } else {
            printf "  CACHED    - %6.3f out of %6.3f GB - s3://genomeark/$filename\n", $newsize / 1024 / 1024 / 1024, $filesize / 1024 / 1024 / 1024;
        }

        printf "  SUMMARIZE - gzip -dc | seqrequester > $sumname\n";
        system("gzip -dc $df 2> $df.gzip.err | $seqrequester summarize - > $sumname 2> $df.seqrequester.err");
        #system("rm -f $df");
    }

    loadSummary($filesize, $filename);
}



#
#  Given a species_name, a genomic data type and a list of seqFiles
#  ("filesize datafile" separated by \0 bytes), compute a
#  scalaing factor that will convert a filesize into an estimated number of
#  bases in the (compressed) file.
#
sub estimateRawDataScaling ($$$$$) {
    my $data     = shift @_;
    my $name     = $$data{"name_"};
    my $type     = shift @_;
    my $files    = shift @_;
    my $errors   = shift @_;
    my $download = shift @_;

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
    printf "\n";

    #  If needed, download a subset of the data from three files and use that to
    #  estimate the total number of bases in the file.

    downloadPartAndSummarize($size1, $file1, $download);    writeSummaries($name);
    downloadPartAndSummarize($size2, $file2, $download);    writeSummaries($name);
    downloadPartAndSummarize($size3, $file3, $download);    writeSummaries($name);

    #  Get an estimate of (or, rarely, the actual) number of bases in each file.

    my $scaling;

    if (!exists($dataBases{$file1}) ||        #  Fail if any of the size computations failed.
        !exists($dataBases{$file2}) ||
        !exists($dataBases{$file3})) {
        if ($download) {
            push @$errors, "  FAILED to estimate sizes for:  (will use approximate scaling estimate)\n";
            push @$errors, "    File 1 $file1\n"   if (!exists($dataBases{$file1}));
            push @$errors, "    File 2 $file2\n"   if (!exists($dataBases{$file2}));
            push @$errors, "    File 3 $file3\n"   if (!exists($dataBases{$file3}));
        }
        else {
            push @$errors, "  WARNING: Size estimates not generated for:  (will use approximate scaling estimate)\n";
            push @$errors, "    File 1 $file1\n"   if (!exists($dataBases{$file1}));
            push @$errors, "    File 2 $file2\n"   if (!exists($dataBases{$file2}));
            push @$errors, "    File 3 $file3\n"   if (!exists($dataBases{$file3}));
        }

        $scaling = 0.0   if ($file1 =~ m!genomic_data/bionano!);
        $scaling = 1.5   if ($file1 =~ m!genomic_data/10x!);           #  Bimodal, ~1.4 and ~1.8
        $scaling = 1.5   if ($file1 =~ m!genomic_data/arima!);         #  Bimodal, ~1.5 and ~1.8
        $scaling = 1.5   if ($file1 =~ m!genomic_data/dovetail!);
        $scaling = 1.8   if ($file1 =~ m!genomic_data/illumina!);
        $scaling = 0.5   if ($file1 =~ m!genomic_data/pacbio!);
        $scaling = 1.1   if ($file1 =~ m!genomic_data/pacbio_hifi!);
        $scaling = 1.6   if ($file1 =~ m!genomic_data/phase!);
    }
    else {
        my $sumbases = $dataBases{$file1} + $dataBases{$file2} + $dataBases{$file3};
        my $sumsizes = $size1             + $size2             + $size3;

        #  Limit precision to avoid differences in rounding
        $scaling = int(10000 * $sumbases / $sumsizes) / 10000.0;
    }

    $$data{"data_${type}_scale"} = $scaling;
}

return(1);
