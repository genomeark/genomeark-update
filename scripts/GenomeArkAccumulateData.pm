package GenomeArkAccumulateData;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(accumulateData);

use strict;
use warnings;

#
#  Given input filesize and filename, accumulate the file into three lists
#  separated by datatype:
#
#    seqFiles{datatype} - list of filesizes and filenames for this datatype
#    seqBytes{datatype} - total size of files for this datatype
#    seqIndiv{datatype} - NUL separated list of individuals with data for this datatype
#

sub accumulateData ($$$$$$) {
    my $filesize = shift @_;   #  Input:  size in bytes of the data file.
    my $filename = shift @_;   #  Input:  name of the data file.
    my $seqFiles = shift @_;   #  In/Out: list of file names by sequence type.  (technically "size name")
    my $seqBytes = shift @_;   #  In/Out: sum  of file sizes by sequence type.
    my $seqIndiv = shift @_;   #  In/Out: list of individual by sequence type.
    my $errors   = shift @_;   #  In/Out: list of potential errors.

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


    #  Nanopore
    #
    #if ($filename =~ m!/genomic_data/ont/!) {
    #}
    #  Nanopore
    #
    if ($filename =~ m!/genomic_data/nanopore/!) {
        #return if ($filename =~ m/txt$/);
        #return if ($filename =~ m/txt$/);
        #return if ($filename =~ m/txt$/);

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"nanopore"} .= $sf;
            $$seqBytes{"nanopore"} += $sb;
            $$seqIndiv{"nanopore"} .= $si;
        } else {
            push @$errors, "  Unknown nanopore file type in '$filename'\n";
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
        #if    ($filename =~ m/\.ccs\..*\.bam.pbi$/) {
        #    $$seqBytes{"pbhifi"} += $sb;
        #}
        #elsif ($filename =~ m/\.ccs\..*\.bam$/) {
        #    $$seqFiles{"pbhifi"} .= $sf;
        #    $$seqBytes{"pbhifi"} += $sb;
        #    $$seqIndiv{"pbhifi"} .= $si;
        #}

        #  Older non-barcode data.
        #elsif ($filename =~ m/\.ccs\.bam\.pbi$/) {
        #    $$seqBytes{"pbhifi"} += $sb;
        #}
        #elsif ($filename =~ m/\.ccs\.bam\.bai$/) {
        #    $$seqBytes{"pbhifi"} += $sb;
        #}
        #elsif ($filename =~ m/\.ccs\.bam$/) {
        #    $$seqBytes{"pbhifi"} += $sb;
        #}
        #elsif ($filename =~ m/\.Q20\.fastq$/) {
        #    $$seqFiles{"pbhifi"} .= $sf;
        #    $$seqBytes{"pbhifi"} += $sb;
        #    $$seqIndiv{"pbhifi"} .= $si;
        #}

        #  Newer non-barcode data - this should all be renamed.
        #elsif ($filename =~ m/\.reads\.bam\.pbi$/) {
        #    $$seqBytes{"pbhifi"} += $sb;
        #}
        #elsif ($filename =~ m/\.reads\.bam$/) {
        #    $$seqBytes{"pbhifi"} += $sb;
        #}
        #elsif ($filename =~ m/\.hifi_reads\.bam$/) {
        #    $$seqBytes{"pbhifi"} += $sb;
        #}
        #elsif ($filename =~ m/\.hifi_reads\.fastq\.gz$/) {
        #    $$seqFiles{"pbhifi"} .= $sf;
        #    $$seqBytes{"pbhifi"} += $sb;
        #    $$seqIndiv{"pbhifi"} .= $si;
        #}

        #  Normal CLR data.
        if    ($filename =~ m/subreads.bam\.pbi$/) {
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

        #  Filtered data

        if    (($filename =~ m/\.hifi_reads\.bam$/) ||
               ($filename =~ m/\.ccs.*\.bam\.pbi$/) ||
               ($filename =~ m/\.ccs.*\.bam\.bai$/) ||
               ($filename =~ m/\.ccs\.bc.*\.bam\.pbi$/) ||
               ($filename =~ m/\.ccs\.bc.*\.bam\.bai$/)) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif (($filename =~ m/\.hifi_reads\.fastq\.gz$/) ||
               ($filename =~ m/\.Q20\.fastq$/) ||
               ($filename =~ m/\.Q20\.fastq\.gz$/) ||
               ($filename =~ m/\.ccs\.bam$/) ||
               ($filename =~ m/\.ccs\.bc.*\.bam$/)) {
            $$seqFiles{"pbhifi"} .= $sf;
            $$seqBytes{"pbhifi"} += $sb;
            $$seqIndiv{"pbhifi"} .= $si;
        }

        #  Ignore unfiltered data.

        elsif (($filename =~ m/\.subreads\.bam\.pbi$/) ||
               ($filename =~ m/\.subreads\.bam\.bai$/) ||
               ($filename =~ m/\.reads\.bam\.pbi$/) ||
               ($filename =~ m/\.reads\.bam\.bai$/)) {
            $$seqBytes{"pbhifi"} += $sb;
        }
        elsif (($filename =~ m/\.subreads\.bam$/) ||
               ($filename =~ m/\.reads\.bam$/)) {
            $$seqBytes{"pbhifi"} += $sb;
        }

        #  Otherwise report confusion.

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
