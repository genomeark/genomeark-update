package GenomeArkAccumulateData;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(scanDataName accumulateData);

use strict;
use warnings;

#
#  Given input filesecs and filename, scanDataName updates the 'modified' date
#  for the species, and builds a private hash of the filenames present.  This
#  is used later to check for bam files without a corresponding fastq file.
#
#  Given input filesize and filename, accumulateData places the filename into
#  three lists separated by datatype:
#    seqFiles{datatype} - list of filesizes and filenames for this datatype
#    seqBytes{datatype} - total size of files for this datatype
#    seqIndiv{datatype} - NUL separated list of individuals with data for this datatype
#

my %filenames;

sub scanDataName ($$$$) {
    my $filesecs = shift @_;   #  Input:  upload time of the data file.
    my $filesize = shift @_;   #  Input:  size in bytes of the data file.
    my $filename = shift @_;   #  Input:  name of the data file.
    my $data     = shift @_;   #  In/Out: data for this species.

    if (!exists($filenames{$$data{"name_"}})) {   #  If this is the first time we've seen
        undef %filenames;                         #  this species, wipe all data.
        $filenames{$$data{"name_"}} = 1;
    }

    #  Save the name of any pacbio_hifi fastq files.
    if ($filename =~ m!pacbio_hifi/(.*)\.f(ast){0,1}q$!) {
        $filenames{$1} = 1;
    }
}


sub saveDataDate ($$) {
    my $filesecs = shift @_;   #  Input:  upload time of the data file.
    my $data     = shift @_;   #  In/Out: data for this species.

    if ($$data{"last_updated"} < $filesecs) {
        $$data{"last_updated"} = $filesecs;
    }
    if ((! exists($$data{"last_raw_data"})) ||
        ($$data{"last_raw_data"} < $filesecs)) {   #  If this isn't set, Raw Data shows
        $$data{"last_raw_data"} = $filesecs;       #  "No data.".
    }
}


sub accumulateData ($$$$$$$$) {
    my $filesecs = shift @_;   #  Input:  upload time of the data file.
    my $filesize = shift @_;   #  Input:  size in bytes of the data file.
    my $filename = shift @_;   #  Input:  name of the data file.
    my $data     = shift @_;   #  In/Out: data for this species.
    my $seqFiles = shift @_;   #  In/Out: list of file names by sequence type.  (technically "size name")
    my $seqBytes = shift @_;   #  In/Out: sum  of file sizes by sequence type.
    my $seqIndiv = shift @_;   #  In/Out: list of individual by sequence type.
    my $errors   = shift @_;   #  In/Out: list of potential errors.
    my $warning  = 0;

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
            saveDataDate($filesecs, $data);
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
            saveDataDate($filesecs, $data);
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
            saveDataDate($filesecs, $data);
        }
        elsif ($filename =~ m/bnx.gz/) {
            $$seqBytes{"bionano"} += $sb;
            $$seqIndiv{"bionano"} .= $si;
            saveDataDate($filesecs, $data);
        }
        elsif ($filename =~ m/bnx/) {
            $$seqBytes{"bionano"} += $sb;
            $$seqIndiv{"bionano"} .= $si;
            saveDataDate($filesecs, $data);
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
            saveDataDate($filesecs, $data);
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
            saveDataDate($filesecs, $data);
        }
        else {
            push @$errors, "  Unknown illumina file type in '$filename'\n";
        }
        return;
    }


    #  Oxford Nanopore
    #
    if ($filename =~ m!/genomic_data/ont/!) {
        #return if ($filename =~ m/txt$/);
        #return if ($filename =~ m/txt$/);
        #return if ($filename =~ m/txt$/);

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"ont"} .= $sf;
            $$seqBytes{"ont"} += $sb;
            $$seqIndiv{"ont"} .= $si;
            saveDataDate($filesecs, $data);
        } else {
            push @$errors, "  Unknown ont file type in '$filename'\n";
        }
        return;
    }


    if ($filename =~ m!/genomic_data/ont_duplex/!) {
        #return if ($filename =~ m/txt$/);
        #return if ($filename =~ m/txt$/);
        #return if ($filename =~ m/txt$/);

        if (($filename =~ m/fastq.gz/) ||
            ($filename =~ m/fq.gz/)) {
            $$seqFiles{"ontduplex"} .= $sf;
            $$seqBytes{"ontduplex"} += $sb;
            $$seqIndiv{"ontduplex"} .= $si;
            saveDataDate($filesecs, $data);
        } else {
            push @$errors, "  Unknown ont_duplex file type in '$filename'\n";
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

        #  Normal CLR data.
        if    ($filename =~ m/subreads.bam\.pbi$/) {
            $$seqBytes{"pacbio"} += $sb;
            saveDataDate($filesecs, $data);
        }
        elsif ($filename =~ m/subreads.bam\.bai$/) {
            $$seqBytes{"pacbio"} += $sb;
            saveDataDate($filesecs, $data);
        }
        elsif ($filename =~ m/subreads.bam$/) {
            $$seqFiles{"pacbio"} .= $sf;
            $$seqBytes{"pacbio"} += $sb;
            $$seqIndiv{"pacbio"} .= $si;
            saveDataDate($filesecs, $data);
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

        if    (($filename =~ m/\.hifi_reads\.f(ast){0,1}q\.gz$/) ||
               ($filename =~ m/\.reads\.f(ast){0,1}q\.gz$/) ||
               ($filename =~ m/\.ccs.bc.*\.f(ast){0,1}q\.gz$/) ||
               ($filename =~ m/\.Q20\.f(ast){0,1}q\.gz$/)) {
            $$seqFiles{"pacbiohifi_fqgz"} .= $sf;
            $$seqBytes{"pacbiohifi_fqgz"} += $sb;
            $$seqIndiv{"pacbiohifi_fqgz"} .= $si;
            saveDataDate($filesecs, $data);
        }

        elsif (($filename =~ m/\.hifi_reads\.f(ast){0,1}q$/) ||
               ($filename =~ m/\.reads\.f(ast){0,1}q$/) ||
               ($filename =~ m/\.ccs.bc.*\.f(ast){0,1}q$/) ||
               ($filename =~ m/\.Q20\.f(ast){0,1}q$/)) {
            if ($warning) {
                push @$errors, "  WARNING: uncompressed pacbio_hifi '$filename'\n";
            }
            #$$seqFiles{"pacbiohifi"} .= $sf;
            #$$seqBytes{"pacbiohifi"} += $sb;
            #$$seqIndiv{"pacbiohifi"} .= $si;
        }

        elsif (($filename =~ m/\.hifi_reads\.bam\.bai$/) ||
               ($filename =~ m/\.hifi_reads\.bam\.pbi$/) ||
               ($filename =~ m/\.ccs.bam\.pbi$/) ||
               ($filename =~ m/\.ccs.bam\.bai$/) ||
               ($filename =~ m/\.ccs.bc.*\.bam\.pbi$/) ||
               ($filename =~ m/\.ccs.bc.*\.bam\.bai$/)) {
            $$seqBytes{"pacbiohifi_bam"} += $sb;
            saveDataDate($filesecs, $data);
        }

        elsif (($filename =~ m/\.hifi_reads\.bam$/) ||
               ($filename =~ m/\.ccs.bam$/) ||
               ($filename =~ m/\.ccs.bc.*\.bam$/)) {
            $$seqFiles{"pacbiohifi_bam"} .= $sf;
            $$seqBytes{"pacbiohifi_bam"} += $sb;
            $$seqIndiv{"pacbiohifi_bam"} .= $si;
            saveDataDate($filesecs, $data);
        }

        #  Ignore unfiltered data, but warn if there isn't a corresponding fastq for it.

        elsif (($filename =~ m/\.subreads\.bam\.pbi$/) ||
               ($filename =~ m/\.reads\.bam\.pbi$/)) {
            $$seqBytes{"pacbiohifi_clr"} += $sb;
            saveDataDate($filesecs, $data);
        }
        elsif (($filename =~ m/\.subreads\.bam\.bai$/) ||
               ($filename =~ m/\.reads\.bam\.bai$/)) {
            $$seqBytes{"pacbiohifi_clr"} += $sb;
            saveDataDate($filesecs, $data);
        }
        elsif (($filename =~ m/\.subreads\.bam$/) ||
               ($filename =~ m/\.reads\.bam$/)) {
            $$seqFiles{"pacbiohifi_clr"} .= $sf;
            $$seqBytes{"pacbiohifi_clr"} += $sb;
            $$seqIndiv{"pacbiohifi_clr"} .= $si;
            saveDataDate($filesecs, $data);

            #  Check that this bam has a corresponding fastq file.
            #    Accipiter_gentilis
            #if ($filename =~ m!pacbio_hifi/(.*)\.bam$!) {
            #    if (! exists($filenames{$1})) {
            #        push @$errors, "  No filtered hifi for '$filename'\n";
            #    }
            #}
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
            saveDataDate($filesecs, $data);
        }
        else {
            push @$errors, "  Unknown phase file type in '$filename'\n";
        }
        return;
    }


    push @$errors, "  Unknown data type in '$filename'\n";
}

return(1);
