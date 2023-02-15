package GenomeArkAccumulateData;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(rememberHiFiName accumulateData);

use strict;
use warnings;

use GenomeArkUpdate;

#
#  rememberHiFiName saves all the hifi fastq data names in a private hash.
#  This is used later to check for bam files without a corresponding fastq
#  file.
#

my %hifinames;   #  Doesn't appear to be used...!

sub rememberHiFiName ($$) {
    my $filename = shift @_;   #  Input:  name of the data file.
    my $data     = shift @_;   #  In/Out: data for this species.

    return  if (! isGenomicDataFile($filename));

    if (!exists($hifinames{$$data{"name_"}})) {   #  If this is the first time we've seen
        undef %hifinames;                         #  this species, wipe all data.
        $hifinames{$$data{"name_"}} = 1;
    }

    #  Save the name of any pacbio_hifi fastq files.
    if ($filename =~ m!pacbio_hifi/(.*)\.f(ast){0,1}q$!) {
        $hifinames{$1} = 1;
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


#
#  Given input filesize and filename, accumulateData places the filename into
#  three lists separated by datatype:
#    seqFiles{datatype} - list of filesizes and filenames for this datatype
#    seqBytes{datatype} - total size of files for this datatype
#    seqIndiv{datatype} - NUL separated list of individuals with data for this datatype
#
sub accumulateData ($$$$$$$$) {
    my $filesecs = shift @_;   #  Input:  upload time of the data file.
    my $filesize = shift @_;   #  Input:  size in bytes of the data file.
    my $filename = shift @_;   #  Input:  name of the data file.

    my $data     = shift @_;   #  In/Out: data for this species.
    my $tiFiles  = shift @_;   #  In/Out: list of file names by sequence type.  (technically "type:indiv size name")
    my $tiBytes  = shift @_;   #  In/Out: sum  of file sizes by sequence type.
    my $tiIndiv  = shift @_;   #  In/Out: list of individual by sequence type.  (technically "type:indiv")
    my $errors   = shift @_;   #  In/Out: list of potential errors.
    my $warning  = 0;

    return  if (! isGenomicDataFile($filename));

    my $sName;   #  Name of the species.
    my $iName;   #  Name of the individual.

    if ($filename =~ m!species/(.*)/(.*)/genomic_data!) {
        $sName = $1;
        $iName = $2;
    } else {
        die "failed to parse species name and individual from '$filename'\n";
    }

    if ($filename =~ m!/genomic_data/10x/!) {
        return if ($filename =~ m/txt$/);

        if (($filename =~ m/fastq.gz$/) ||
            ($filename =~ m/fq.gz$/)) {
            $$tiFiles{"10x:$iName"} .= "10x:$iName $filesize $filename\0";
            $$tiBytes{"10x:$iName"} += $filesize;
            $$tiIndiv{"10x:$iName"} .= "$sName/$iName\0";
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

        if (($filename =~ m/fastq.gz$/) ||
            ($filename =~ m/fq.gz$/)) {
            $$tiFiles{"arima:$iName"} .= "arima:$iName $filesize $filename\0";
            $$tiBytes{"arima:$iName"} += $filesize;
            $$tiIndiv{"arima:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        }
        elsif (($filename =~ m/cram$/) ||
               ($filename =~ m/cram.crai$/)) {
            $$tiFiles{"arima:$iName"} .= "arima:$iName $filesize $filename\0"   if ($filename =~ m/cram$/);
            $$tiBytes{"arima:$iName"} += $filesize;
            $$tiIndiv{"arima:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        }
        else {
            push @$errors, "  Unknown arima file type in '$filename'\n";
        }
        return;
    }


    if ($filename =~ m!/genomic_data/bionano/!) {
        return if ($filename =~ m/txt$/);

        if    ($filename =~ m/cmap.gz$/) {
            #$$tiFiles{"bionano:$iName"} .= "bionano:$iName $filesize $filename\0";
            $$tiBytes{"bionano:$iName"} += $filesize;
            $$tiIndiv{"bionano:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        }
        elsif ($filename =~ m/cmap$/) {
            #$$tiFiles{"bionano:$iName"} .= "bionano:$iName $filesize $filename\0";
            $$tiBytes{"bionano:$iName"} += $filesize;
            $$tiIndiv{"bionano:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        }
        elsif ($filename =~ m/bnx.gz$/) {
            $$tiFiles{"bionano:$iName"} .= "bionano:$iName $filesize $filename\0";
            $$tiBytes{"bionano:$iName"} += $filesize;
            $$tiIndiv{"bionano:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        }
        elsif ($filename =~ m/bnx$/) {
            $$tiFiles{"bionano:$iName"} .= "bionano:$iName $filesize $filename\0";
            $$tiBytes{"bionano:$iName"} += $filesize;
            $$tiIndiv{"bionano:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        }
        else {
            push @$errors, "  Unknown bionano file type in '$filename'\n";
        }
        return;
    }

    if ($filename =~ m!/genomic_data/dovetail/!) {
        return if ($filename =~ m/re_bases.txt/);

        if (($filename =~ m/fastq.gz$/) ||
            ($filename =~ m/fq.gz$/)) {
            $$tiFiles{"dovetail:$iName"} .= "dovetail:$iName $filesize $filename\0";
            $$tiBytes{"dovetail:$iName"} += $filesize;
            $$tiIndiv{"dovetail:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        }
        else {
            push @$errors, "  Unknown dovetail file type in '$filename'\n";
        }
        return;
    }

    if ($filename =~ m!/genomic_data/illumina/!) {
        return if ($filename =~ m/txt$/);

        if (($filename =~ m/fastq.gz$/) ||
            ($filename =~ m/fq.gz$/)) {
            $$tiFiles{"illumina:$iName"} .= "illumina:$iName $filesize $filename\0";
            $$tiBytes{"illumina:$iName"} += $filesize;
            $$tiIndiv{"illumina:$iName"} .= "$sName/$iName\0";
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

        if (($filename =~ m/fastq.gz$/) ||
            ($filename =~ m/fq.gz$/)) {
            $$tiFiles{"ont:$iName"} .= "ont:$iName $filesize $filename\0";
            $$tiBytes{"ont:$iName"} += $filesize;
            $$tiIndiv{"ont:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        } elsif ($filename =~ m/bam$/) {
        } elsif ($filename =~ m/fast5$/) {
            #$$tiFiles{"ont:$iName"} .= {"ont:$iName $filesize $filename\0";
            #$$tiBytes{"ont:$iName"} += $filesize;
            #$$tiIndiv{"ont:$iName"} .= "$sName/$iName\0";
        } else {
            push @$errors, "  Unknown ont file type in '$filename'\n";
        }
        return;
    }


    if ($filename =~ m!/genomic_data/ont_duplex/!) {
        #return if ($filename =~ m/txt$/);
        #return if ($filename =~ m/txt$/);
        #return if ($filename =~ m/txt$/);

        if (($filename =~ m/fastq.gz$/) ||
            ($filename =~ m/fq.gz$/)) {
            $$tiFiles{"ontduplex:$iName"} .= "ontduplex:$iName $filesize $filename\0";
            $$tiBytes{"ontduplex:$iName"} += $filesize;
            $$tiIndiv{"ontduplex:$iName"} .= "$sName/$iName\0";
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
            $$tiBytes{"pacbio:$iName"} += $filesize;
            saveDataDate($filesecs, $data);
        }
        elsif ($filename =~ m/subreads.bam\.bai$/) {
            $$tiBytes{"pacbio:$iName"} += $filesize;
            saveDataDate($filesecs, $data);
        }
        elsif ($filename =~ m/subreads.bam$/) {
            $$tiFiles{"pacbio:$iName"} .= "pacbio:$iName $filesize $filename\0";
            $$tiBytes{"pacbio:$iName"} += $filesize;
            $$tiIndiv{"pacbio:$iName"} .= "$sName/$iName\0";
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
            $$tiFiles{"pacbiohifi_fqgz:$iName"} .= "pacbiohifi_fqgz:$iName $filesize $filename\0";
            $$tiBytes{"pacbiohifi_fqgz:$iName"} += $filesize;
            $$tiIndiv{"pacbiohifi_fqgz:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        }

        elsif (($filename =~ m/\.hifi_reads\.f(ast){0,1}q$/) ||
               ($filename =~ m/\.reads\.f(ast){0,1}q$/) ||
               ($filename =~ m/\.ccs.bc.*\.f(ast){0,1}q$/) ||
               ($filename =~ m/\.Q20\.f(ast){0,1}q$/)) {
            if ($warning) {
                push @$errors, "  WARNING: uncompressed pacbio_hifi '$filename'\n";
            }
            #$$tiFiles{"pacbiohifi:$iName"} .= {"pacbiohifi:$iName $filesize $filename\0";
            #$$tiBytes{"pacbiohifi:$iName"} += $filesize;
            #$$tiIndiv{"pacbiohifi:$iName"} .= "$sName/$iName\0";
        }

        elsif (($filename =~ m/\.hifi_reads\.bam\.bai$/) ||
               ($filename =~ m/\.hifi_reads\.bam\.pbi$/) ||
               ($filename =~ m/\.ccs.bam\.pbi$/) ||
               ($filename =~ m/\.ccs.bam\.bai$/) ||
               ($filename =~ m/\.ccs.bc.*\.bam\.pbi$/) ||
               ($filename =~ m/\.ccs.bc.*\.bam\.bai$/)) {
            $$tiBytes{"pacbiohifi_bam:$iName"} += $filesize;
            saveDataDate($filesecs, $data);
        }

        elsif (($filename =~ m/\.hifi_reads\.bam$/) ||
               ($filename =~ m/\.ccs.bam$/) ||
               ($filename =~ m/\.ccs.bc.*\.bam$/)) {
            $$tiFiles{"pacbiohifi_bam:$iName"} .= "pacbiohifi_bam:$iName $filesize $filename\0";
            $$tiBytes{"pacbiohifi_bam:$iName"} += $filesize;
            $$tiIndiv{"pacbiohifi_bam:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);
        }

        #  Ignore unfiltered data, but warn if there isn't a corresponding fastq for it.

        elsif (($filename =~ m/\.subreads\.bam\.pbi$/) ||
               ($filename =~ m/\.reads\.bam\.pbi$/)) {
            $$tiBytes{"pacbiohifi_clr:$iName"} += $filesize;
            saveDataDate($filesecs, $data);
        }
        elsif (($filename =~ m/\.subreads\.bam\.bai$/) ||
               ($filename =~ m/\.reads\.bam\.bai$/)) {
            $$tiBytes{"pacbiohifi_clr:$iName"} += $filesize;
            saveDataDate($filesecs, $data);
        }
        elsif (($filename =~ m/\.subreads\.bam$/) ||
               ($filename =~ m/\.reads\.bam$/)) {
            $$tiFiles{"pacbiohifi_clr:$iName"} .= "pacbiohifi_clr:$iName $filesize $filename\0";
            $$tiBytes{"pacbiohifi_clr:$iName"} += $filesize;
            $$tiIndiv{"pacbiohifi_clr:$iName"} .= "$sName/$iName\0";
            saveDataDate($filesecs, $data);

            #  Check that this bam has a corresponding fastq file.
            #    Accipiter_gentilis
            #if ($filename =~ m!pacbio_hifi/(.*)\.bam$!) {
            #    if (! exists($hifinames{$1})) {
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

        if (($filename =~ m/fastq.gz$/) ||
            ($filename =~ m/fq.gz$/)) {
            $$tiFiles{"phase:$iName"} .= "phase:$iName $filesize $filename\0";
            $$tiBytes{"phase:$iName"} += $filesize;
            $$tiIndiv{"phase:$iName"} .= "$sName/$iName\0";
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
