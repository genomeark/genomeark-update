package GenomeArkGenomicData;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(accumulateData estimateRawDataScaling);

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




sub parseSummary ($$$) {
    my $name  = shift @_;
    my $size  = shift @_;
    my $file  = shift @_;
    my $bases = 0;

    if (-e "$name.summary") {
        open(ST, "< $name.summary");
        while (<ST>) {
            if (m/^G=(\d+)\s/) {
                $bases = $1;
                last;
            }
        }
        close(ST);
    }

    if (($bases == 0) && ($downloadRAW == 0)) {
        $bases = $size;
    }

    return($bases);
}


sub downloadPipeSummarize ($$$) {
    my $name  = shift @_;
    my $ndir  = dirname($name);
    my $size  = shift @_;
    my $file  = shift @_;

    unlink "$name.summary"   if (-z "$name.summary");

    #  If the summary exists, parse it and return the number of bases.

    if (-e "$name.summary") {
        return(parseSummary($name, $size, $file));
    }

    #  If the download exists, summarize it and return the number of bases,

    if (-e "downloads/$name") {
        return(downloadAndSummarize($name, $size, $file));
    }

    #  Otherwise, fetch the data directly into the summarizer.

    if ($name =~ m/bam$/) {
        printf "FETCH data - size %6.3f GB\n", $size / 1024 / 1024 / 1024;
        printf "  aws --no-progress --no-sign-request s3 cp \\\n";
        printf "    s3://genomeark/$name | bamsummarize > \\\n";
        printf "\n";

        system("mkdir -p $ndir")   if (! -e $ndir);
        system("aws --no-progress --no-sign-request s3 cp s3://genomeark/$name - | gzip -dc | $bamsummarize > $name.summary")  if ($downloadRAW);
    }
    else {
        printf "FETCH data - size %6.3f GB\n", $size / 1024 / 1024 / 1024;
        printf "  aws --no-progress --no-sign-request s3 cp \\\n";
        printf "    s3://genomeark/$name | \\\n";
        printf "         downloads/$name\n";
        printf "\n";

        system("mkdir -p $ndir")   if (! -e $ndir);
        system("aws --no-progress --no-sign-request s3 cp s3://genomeark/$name - | $seqrequester summarize - > $name.summary")  if ($downloadRAW);
    }

    return(parseSummary($name, $size, $file));
}


sub downloadAndSummarize ($$$) {
    my $name  = shift @_;
    my $ndir  = dirname($name);
    my $size  = shift @_;
    my $file  = shift @_;
    my $bases = 0;

    unlink "$name.summary"   if (-z "$name.summary");

    #  If the summary exists, parse it and return the number of bases.

    if (-e "$name.summary") {
        return(parseSummary($name, $size, $file));
    }

    #  Fetch the data if it doesn't exist already.

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

    return(parseSummary($name, $size, $file));
}



#
#  Given a species_name, a genomic data type and a list of seqFiles
#  ("filesize species_name/individual" separated by \0 bytes), compute a
#  scalaing factor that will convert a filesize into an estimated number of
#  bases in the (compressed) file.
#
sub estimateRawDataScaling ($$$) {
    my $data     =             shift @_;
    my $name     = $$data{"name_"};
    my $type     =             shift @_;
    my $files    =             shift @_;

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

    print " - $type:\n";

    for (my $ii=0; $ii<$filesLen; $ii++) {
        my ($size, $name) = split '\s+', $files[$ii];

        if (($ii == $f1) ||
            ($ii == $f2) ||
            ($ii == $f3)) {
            printf "    * %12d %s\n", $size, $name;
        } else {
            printf "      %12d %s\n", $size, $name;
        }
    }

    $$data{"data_${type}_scale"} = 0.0;
    return;

    $bases1 = downloadAndSummarize($name1, $size1, $f1);   #  Compute (if needed) and return the number
    $bases2 = downloadAndSummarize($name2, $size2, $f2);   #  of bases in each of these files.
    $bases3 = downloadAndSummarize($name3, $size3, $f3);

    if (($bases1 == 0) ||        #  Fail if any of the size computations failed.
        ($bases2 == 0) ||
        ($bases3 == 0)) {
        print "FAILED TO ESTIMATE SIZES.\n";
        print "  1 - $bases1 - $name1\n";
        print "  2 - $bases2 - $name2\n";
        print "  3 - $bases3 - $name3\n";
        exit(1);
    }

    my $scaling  = 1.0;
    $scaling = ($bases1 + $bases2 + $bases3) / ($size1 + $size2 + $size3);   #  Compute scaling, limit precision
    $scaling = int($scaling * 10000) / 10000;                                #  to prevent churn from bad math.

    $$data{"data_${type}_scale"} = $scaling;
    #return($scaling);
}





#
#  Given input filesize and filename, accumulate the file into
#    seqFiles{datatype} - list of filesizes and filenames for this datatype
#    seqBytes{datatype} - total size of files for this datatype
#    seqIndiv{datatype} - NUL separated list of individuals with data for this datatype
#

sub accumulateData ($$$$$$) {
    my $filesize = shift @_;
    my $filename = shift @_;
    my $seqFiles = shift @_;
    my $seqBytes = shift @_;
    my $seqIndiv = shift @_;
    my $errors   = shift @_;

    my $sName;   #  Name of the species.
    my $iName;   #  Name of the individual.
    my $sf;      #  Output filename, appended to seqFiles
    my $sb;      #  Output filesize, added to seqBytes
    my $si;      #  Output individual, appended to seqIndiv

    if ($filename =~ m!species/(.*)/(.*)/genomic_data!) {
        $sName = $1;
        $iName = $2;
        $sf    = "$filesize $filename\0";
        $sb    = $filesize;
        $si    = "$sName/$iName\0";
    } else {
        die "failed to parse species name and individual from '$filename'\n";
    }


    if ($filename =~ m!/genomic_data/10x/!) {
        return if ($filename =~ m/txt$/);

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"10x"} .= $sf;
            $$seqBytes{"10x"} += $sb;
            $$seqIndiv{"10x"} .= $si;
        }
        else {
            push @$errors, "  Unknown 10x file type in '$filename'\n";
        }
        return;
    }


    if (($filename =~ m!/genomic_data/arima/!) ||
        ($filename =~ m!/genomic_data/baylor-hic/!) ||
        ($filename =~ m!/genomic_data/hic/!)) {
        return if ($filename =~ m/re_bases.txt/);
        return if ($filename =~ m/re_enz.txt/);

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"arima"} .= $sf;
            $$seqBytes{"arima"} += $sb;
            $$seqIndiv{"arima"} .= $si;
        }
        else {
            push @$errors, "  Unknown arima file type in '$filename'\n";
        }
        return;
    }


    if ($filename =~ m!/genomic_data/bionano/!) {
        return if ($filename =~ m/txt$/);

        if      ($filename =~ m/cmap/) {
            $$seqFiles{"bionano"} .= $sf;
            $$seqBytes{"bionano"} += $sb;
            $$seqIndiv{"bionano"} .= $si;
        }
        elsif ($filename =~ m/bnx.gz/) {
            $$seqBytes{"bionano"} += $sb;
            $$seqIndiv{"bionano"} .= $si;
        }
        elsif ($filename =~ m/bnx/) {
            $$seqBytes{"bionano"} += $sb;
            $$seqIndiv{"bionano"} .= $si;
        }
        else {
            push @$errors, "  Unknown bionano file type in '$filename'\n";
        }
        return;
    }

    if ($filename =~ m!/genomic_data/dovetail/!) {
        return if ($filename =~ m/re_bases.txt/);

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"dovetail"} .= $sf;
            $$seqBytes{"dovetail"} += $sb;
            $$seqIndiv{"dovetail"} .= $si;
        }
        else {
            push @$errors, "  Unknown dovetail file type in '$filename'\n";
        }
        return;
    }

    if ($filename =~ m!/genomic_data/illumina/!) {
        return if ($filename =~ m/txt$/);

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"illumina"} .= $sf;
            $$seqBytes{"illumina"} += $sb;
            $$seqIndiv{"illumina"} .= $si;
        }
        else {
            push @$errors, "  Unknown illumina file type in '$filename'\n";
        }
        return;
    }


    #  Completely ignore PacBio scraps.  They shouldn't even exist anymore.
    #
    if ($filename =~ m!/genomic_data/pacbio/!) {
        return if ($filename =~ m/scraps\./);
    }

    #  PacBio CLR has:
    #    bax.h5
    #    subreads.fasta
    #    subreads.bam
    #    subreads.bam.bai
    #    subreads.bam.pbi
    #
    if ($filename =~ m!/genomic_data/pacbio/!) {
        return if ($filename =~ m/txt$/);

        #  Older barcode data
        if    ($filename =~ m/\.ccs\..*\.bam.pbi$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.ccs\..*\.bam$/) {
            $$seqFiles{"pbhifi"} .= $sf;
            $$seqBytes{"pbhifi"} += $sb;
            $$seqIndiv{"pbhifi"} .= $si;
        }

        #  Older non-barcode data.
        elsif ($filename =~ m/\.ccs\.bam\.pbi$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.ccs\.bam\.bai$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.ccs\.bam$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.Q20\.fastq$/) {
            $$seqFiles{"pbhifi"} .= $sf;
            $$seqBytes{"pbhifi"} += $sb;
            $$seqIndiv{"pbhifi"} .= $si;
        }

        #  Newer non-barcode data - this should all be renamed.
        elsif ($filename =~ m/\.reads\.bam\.pbi$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.reads\.bam$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.hifi_reads\.bam$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.hifi_reads\.fastq\.gz$/) {
            $$seqFiles{"pbhifi"} .= $sf;
            $$seqBytes{"pbhifi"} += $sb;
            $$seqIndiv{"pbhifi"} .= $si;
        }

        #  Normal CLR data.
        elsif ($filename =~ m/subreads.bam\.pbi$/) {
            $$seqBytes{"pbclr"} += $sb;
        }
        elsif ($filename =~ m/subreads.bam\.bai$/) {
            $$seqBytes{"pbclr"} += $sb;
        }
        elsif ($filename =~ m/subreads.bam$/) {
            $$seqFiles{"pbclr"} .= $sf;
            $$seqBytes{"pbclr"} += $sb;
            $$seqIndiv{"pbclr"} .= $si;
        }

        else {
            push @$errors, "  Unknown pacbio file type in '$filename'\n";
        }

        return;
    }

    #  PacBio HiFi has:
    #    reads.bam
    #    reads.bam.pbi
    #
    #    hifi_reads.bam
    #    hifi_reads.fastq.gz
    #
    if ($filename =~ m!/genomic_data/pacbio_hifi/!) {
        return if ($filename =~ m/txt$/);

        #  Older data
        if    ($filename =~ m/\.hifi_reads\.bam$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.hifi_reads\.fastq\.gz$/) {
            $$seqFiles{"pbhifi"} .= $sf;
            $$seqBytes{"pbhifi"} += $sb;
            $$seqIndiv{"pbhifi"} .= $si;
        }

        elsif ($filename =~ m/\.subreads\.bam\.pbi$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.subreads\.bam\.bai$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.subreads\.bam$/) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif ($filename =~ m/\.reads\.bam$/) {
            $$seqFiles{"pbhifi"} .= $sf;
            $$seqBytes{"pbhifi"} += $sb;
            $$seqIndiv{"pbhifi"} .= $si;
        }

        else {
            push @$errors, "  Unknown pacbio_hifi file type in '$filename'\n";
        }

        return;
    }


    if ($filename =~ m!/genomic_data/phase/!) {
        return if ($filename =~ m/re_bases.txt/);

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"phase"} .= $sf;
            $$seqBytes{"phase"} += $sb;
            $$seqIndiv{"phase"} .= $si;
        }
        else {
            push @$errors, "  Unknown phase file type in '$filename'\n";
        }
        return;
    }


    push @$errors, "  Unknown data type in '$filename'\n";
}


return(1);
