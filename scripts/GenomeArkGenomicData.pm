package GenomeArkGenomicData;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(@dataTypes loadSummaries recoverSummaries writeSummaries estimateRawDataScaling);

use strict;
use warnings;

use GenomeArkBionano;
use GenomeArkUtility;

#use Time::Local;
use List::Util qw(min max);
use File::Basename;

my $aws          = "aws --no-sign-request";
my $seqrequester = "./seqrequester/build/bin/seqrequester";

#
#  Summarized data (stored in e.g. 'species/Gallus_gallus/genomic_data.summary') is
#  loaded into five maps from filename to each column in the summary file.
#
#  This data is extended when new data files are summarized (by
#  loadSummaryIfExists() and loadSummaryFromFile()) and is then used to
#  recreate genomic_data.summary at exit.
#

my %fnTypeI;   #  Map filename -> data-type and individual (: separated)
my %fnBytes;   #  Map filename -> number of bytes in a genomic_data sequence file.
my %fnQuant;   #  Map filename -> Amount of raw file summarized (in gigabytes, or "all").
my %fnBases;   #  Map filename -> number of bases in piece summarized.
my %fnReads;   #  Map filename -> number of reads in piece summarized.

#
#  Return the summary file name for some sequence file.  Strips off
#  suffixes, and (optionally) appends the size of the data summarized.
#
sub makeSummaryFileName ($$) {
    my $n = shift @_;
    my $d = dirname($n);
    my $f = basename($n);
    my $s = shift @_;
    my $N;

    $f =~ s/.fasta.gz$//;   #  Remove extensions.  Bases in bam and fastq _should_
    $f =~ s/.fastq.gz$//;   #  be the same, right?
    $f =~ s/.fasta$//;
    $f =~ s/.fastq$//;
    $f =~ s/.bam$//;
    $f =~ s/.fast5$//;

    if (defined($s)) {                    #  In particular, if $s eq "*", this
        $N = "$d/$s-bytes--$f.summary";   #  will match any size - but not
    } else {                              #  the whole data.
        $N = "$d/$f.summary";
    }

    return($N);
}


#
#  Load a single data summary, which MUST exist.  If the summary is for a
#  partial file, scale the summarized values to approximate the whole file.
#
#  If, after all that work, we end up with no estimate, forget everything.
#
sub loadSummaryFromFile ($$$$$$) {
    my $filetype  = shift @_;   #  type:indiv of the file we've summarized.
    my $filesize  = shift @_;   #  Size of the file we've summarized.
    my $filename  = shift @_;   #  Name of the file we've summarized.
    my $sumsize   = shift @_;   #  Amount of the file we summarized.
    my $summary   = shift @_;   #  Name of the summary file.
    my $errors    = shift @_;

    $fnBytes{$filename} = $filesize;
    $fnQuant{$filename} = $sumsize;
    $fnBases{$filename} = 0;
    $fnReads{$filename} = 0;
    $fnTypeI{$filename} = $filetype;

    print "  Import summary from '$summary'\n";

    open(ST, "< $summary") or die;
    while (<ST>) {
        $fnBases{$filename} = $1   if (m/^G=(\d+)\s/);
        $fnReads{$filename} = $1   if (m/^001.000x\s+(\d+)\s/)
    }
    close(ST);

    if ($filesize != $sumsize) {
        $fnBases{$filename} = int($fnBases{$filename} * $filesize / $sumsize);
        $fnReads{$filename} = int($fnReads{$filename} * $filesize / $sumsize);
    }

    if (($fnBases{$filename} == 0) ||
        ($fnReads{$filename} == 0)) {
        push @$errors, "Summary FAILED, no bases or reads found.  Summary saved in '$summary.BAD'.\n";

        unlink "$summary.BAD";
        rename "$summary", "$summary.BAD";

        delete $fnBytes{$filename};
        delete $fnQuant{$filename};
        delete $fnBases{$filename};
        delete $fnReads{$filename};
    }
}

#
#  Load a summary, if any summary exists.
#
sub loadSummaryIfExists ($$$$) {
    my $filetype  = shift @_;
    my $filesize  = shift @_;
    my $filename  = shift @_;
    my $fullname  = makeSummaryFileName($filename, undef);   #  Name of a full summary.
    my $partnames = makeSummaryFileName($filename, "*");     #  Glob to find any partial summaries.
    my $errors    = shift @_;

    my $sumsize   = (exists($fnQuant{$filename})) ? $fnQuant{$filename} : 0;

    my $newsize   = 0;
    my $newname   = "";

    #  If a full summary exists, always use that.
    if (-e $fullname) {
        $newsize = $filesize;
        $newname = $fullname;
    }

    #  Otherwise, search for a partial summary of any size.
    else {
        open(LS, "ls $partnames 2>/dev/null |");
        while (<LS>) {
            chomp;
            if ($_ =~ m!/(\d+)-bytes--.*\.summary$!) {
                if ($1 > $newsize) {
                    $newsize = $1;
                    $newname = $_;
                }
            }
        }
        close(LS);
    }

    #  Load the new summary if it summarized more data than we did.

    if ($newsize > $sumsize) {
        loadSummaryFromFile($filetype, $filesize, $filename, $newsize, $newname, $errors);
    }
}


#
#  Load existing data summaries from 'species/$name/genomic_data.summary'.
#
sub loadSummaries ($$$) {
    my $data    = shift @_;          #  Pointer to global %data.
    my $name    = $$data{"name_"};   #  Name of this species.
    my $tiFiles = shift @_;          #  Map from ti to list of files - all data files are here.
    my $errors  = shift @_;          #  List of output errors.

    #  Clear any data from previous species.

    %fnBytes = ();
    %fnQuant = ();
    %fnBases = ();
    %fnReads = ();
    %fnTypeI = ();

    #
    #  Load precomputed summaries.
    #

    if (-e "species/$name/genomic_data.summary") {
        print "  Load summaries from species/$name/genomic_data.summary\n";

        open(ST, "< species/$name/genomic_data.summary") or die;
        while (<ST>) {
            s/^\s+//;
            s/\s+$//;

            next   if (m/^-/);
            next   if (m/^$/);

            my @v = split '\s+', $_;

            #  If a file description, save the info.

            if (scalar(@v) == 8) {
                my ($type, $indiv, $bytes, $date, $quant, $bases, $reads, $file) = @v;

                $fnBytes{$file} =  $bytes;
                $fnQuant{$file} =  $quant;
                $fnBases{$file} =  $bases;
                $fnReads{$file} =  $reads;
                $fnTypeI{$file} = "$type:$indiv";
            }

            #  If a scaling summary, save the scale.

            if (scalar(@v) == 4) {
                my ($type, $indiv, $scale, $bases) = @v;

                $$data{"data_${type}:${indiv}_scale"} = $scale;   #  Scale factor to convert bytes to bases.
                $$data{"data_${type}:${indiv}_bases"} = $bases;   #  Current estimated bases; updated by scan-bucket.pl
            }
        }
        close(ST);
    }
}


#
#  Scan the list of seqFiles to find any summaries that we don't know about
#  or have been improved since the last run.  This will pull in manually
#  generated summaries and recover a missing (or deleted)
#  genomic_data.summary.
#
sub recoverSummaries ($$$) {
    my $data    = shift @_;          #  Pointer to global %data.
    my $name    = $$data{"name_"};   #  Name of this species.
    my $tiFiles = shift @_;          #  Map from ti to list of files - all data files are here.
    my $errors  = shift @_;          #  List of output errors.

    foreach my $ti (keys %$tiFiles) {
        my @tifs = split '\0', $$tiFiles{$ti};

        foreach my $tif (@tifs) {
            my ($i, $s, $n) = split '\s+', $tif;

            loadSummaryIfExists($i, $s, $n, $errors);
        }
    }
}


#
#  Write all the summary information we know to the precomputed summary file.
#
sub writeSummaries ($) {
    my $data     = shift @_;
    my $name     = $$data{"name_"};

    system("mkdir -p species/$name")   if (! -e "species/$name");

    open(ST, "> species/$name/genomic_data.summary") or die "Failed to open 'species/$name/genomic_data.summary' for writing: $!\n";

    print ST "-  GENOMIC DATA FILE SUMMARIES\n";
    print ST "-\n";
    print ST "-data-type      individual                 bytes            date        quantity           bases           reads   file\n";
    print ST "--------------- ---------------  --------------- --------------- --------------- --------------- ---------------  ----------\n";

    my %types;

    foreach my $file (sort keys %fnBytes) {
        my  $date = 0;   #  For later use, maybe.

        next   if (! exists($fnBytes{$file}));
        next   if (! exists($fnQuant{$file}));
        next   if (! exists($fnBases{$file}));
        next   if (! exists($fnReads{$file}));
        next   if (! exists($fnTypeI{$file}));

        $types{$fnTypeI{$file}}++;

        my ($type, $indiv) = split ':', $fnTypeI{$file};

        printf(ST "%-15s %-15s  %15d %15d %15s %15d %15d  %s\n", $type, $indiv, $fnBytes{$file}, $date, $fnQuant{$file}, $fnBases{$file}, $fnReads{$file}, $file);
    }

    print ST "--------------- ---------------  --------------- --------------- --------------- --------------- ---------------  ----------\n";
    print ST "\n";
    print ST "-  GENOMIC DATA SCALING FACTORS\n";
    print ST "-\n";
    print ST "-datatype       individual       scaling           bases\n";
    print ST "--------------- ---------------  ------- ---------------\n";

    foreach my $ti (sort keys %types) {
        my ($ty, $in) = split ':', $ti;

        my $scale = $$data{"data_${ti}_scale"};
        my $bases = $$data{"data_${ti}_bases"};

        printf(ST "%-15s %-15s  %7.4f %15d\n", $ty, $in, $scale, $bases);
    }

    print ST "--------------- ---------------  ------- ---------------\n";

    close(ST);
}


#
#  If the summary doesn't exist for all the data, download the whole file and summarize.
#
sub downloadFullAndSummarize ($$$$$) {
    my $filetype = shift @_;
    my $filesize = shift @_;
    my $filename = shift @_;
    my $fdir     = dirname($filename);
    my $download = shift @_;
    my $errors   = shift @_;

    my $fullname = makeSummaryFileName($filename, undef);

    #  If the whole summary exists, do not recompute it.

    return   if (exists($fnQuant{$filename}) && $fnQuant{$filename} eq "all");

    die      if (-e $fullname);   #  It should be loaded already!

    #  Fetch the data if it doesn't exist already.

    if ($download == 0) {
        printf "FETCH data - size %6.3f GB -- DISABLED\n", $filesize / 1024 / 1024 / 1024;
        printf "    s3://genomeark/$filename\n";
        printf "\n";
        return;
    }

    if ($filename =~ m/fast5\$/) {
        printf "FETCH data - size %6.3f GB -- CAN'T PROCESS FAST5 FILES\n", $filesize / 1024 / 1024 / 1024;
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
        push @$errors, "DOWNLOAD FAILED for 'downloads/$filename'.\n";
        return;
    }

    #  If a bam or cram, convert to fastq then summarize, otherwise, summarize directly.

    if (($filename =~ m/bam$/) || ($filename =~ m/cram$/)) {
        printf "EXTRACT and SUMMARIZE to $fullname\n";

        system("samtools fasta downloads/$filename | $seqrequester summarize - > $fullname");
    }
    else {
        printf "SUMMARIZE to $fullname\n";

        system("$seqrequester summarize downloads/$filename > $fullname");
    }

    #  Parse the summary to find the number of bases in the dataset.

    loadSummaryFromFile($filetype, $filesize, $filename, $filesize, $fullname, $errors);
}


#
#  If no estimated number of bases exists, download a part of the file and summarize.
#
sub downloadPartAndSummarize ($$$$$$) {
    my $filetype = shift @_;
    my $filesize = shift @_;
    my $filename = shift @_;
    my $fdir     = dirname($filename);
    my $download = shift @_;
    my $downsize = shift @_;
    my $errors   = shift @_;

    #  How much did we summarize already, and how much do we want to summarize?
    #
    #  Generate some summary file names, one for a full summary (to be used
    #  if it exists) and one for the partial summary we're going to make.

    my $oldsize  = (exists($fnQuant{$filename})) ? $fnQuant{$filename} : 0;
    my $sumsize  = $downsize;

    my $fullname = makeSummaryFileName($filename, undef);
    my $partname = makeSummaryFileName($filename, $sumsize);

    #  If we've summarized all the data, we can't improve it.
    if    ($oldsize >= $filesize) {
        return;
    }

    #  If there is a full summary around, use that.
    elsif (-e $fullname) {
        loadSummaryFromFile($filetype, $filesize, $filename, $filesize, $fullname, $errors);
        return;
    }

    #  If there is ALREADY a data file around to summarize, make a full summary.
    elsif (-e "downloads/$filename") {
        downloadFullAndSummarize($filetype, $filesize, $filename, $download, $errors);
        return;
    }

    #  If we've summarized more than requested, don't do more work to get less!
    elsif ($oldsize >= $sumsize) {
        return;
    }

    #
    #  Yay!  Let's download and summarize something!
    #

    if ($download == 0) {
        printf "\n";
        printf "  SKIP  - %6.3f GB out of %6.3f GB - s3://genomeark/$filename\n", 0.0, $filesize / 1024 / 1024 / 1024;

        push @$errors, "  Downloading disabled for file $filename\n";
        return;
    }

    #  If the file isn't hugely bigger than the size we'd sample, just do the whole thing.
    if ($filesize < 2 * $sumsize) {
        downloadFullAndSummarize($filetype, $filesize, $filename, $download, $errors);
        return;
    }

    #  While 's3 cp' is happy writing to stdout, 's3api get-object' is not,
    #  and requires a file name.  Rather annoying, since this is where we
    #  could actually get away with piping it to the summarizer.

    system("mkdir -p downloads/$fdir");   #  Make a place to download the file
    system("mkdir -p           $fdir");   #  and a place to write the summary.

    if (($filename =~ m/bam$/) || ($filename =~ m/cram$/)) {
        my $df = "downloads/" . dirname($filename) . "/${sumsize}-bytes--" . basename($filename);

        if (! -e "$df") {
            printf "  FETCH     - %6.3f out of %6.3f GB - s3://genomeark/$filename\n", $sumsize / 1024 / 1024 / 1024, $filesize / 1024 / 1024 / 1024;
            system("$aws s3api get-object --bucket genomeark --key $filename --range bytes=0-$sumsize $df > $df.err 2>&1")  if (! -e "$df");
        } else {
            printf "  CACHED    - %6.3f out of %6.3f GB - s3://genomeark/$filename\n", $sumsize / 1024 / 1024 / 1024, $filesize / 1024 / 1024 / 1024;
        }

        printf "  SUMMARIZE - samtools | seqrequester > $partname\n";
        system("samtools fasta $df 2> $df.samtools.err | $seqrequester summarize - > $partname 2> $df.seqrequester.err");
        #system("rm -f $df");
    }
    else {
        my $df = "downloads/" . dirname($filename) . "/${sumsize}-bytes--" . basename($filename);

        if (! -e "$df") {
            printf "  FETCH     - %6.3f out of %6.3f GB - s3://genomeark/$filename\n", $sumsize / 1024 / 1024 / 1024, $filesize / 1024 / 1024 / 1024;
            system("$aws s3api get-object --bucket genomeark --key $filename --range bytes=0-$sumsize $df > $df.err 2>&1")  if (! -e "$df");
        } else {
            printf "  CACHED    - %6.3f out of %6.3f GB - s3://genomeark/$filename\n", $sumsize / 1024 / 1024 / 1024, $filesize / 1024 / 1024 / 1024;
        }

        if ($df =~ m/.gz$/) {
            printf "  SUMMARIZE - gzip -dc | seqrequester > $partname\n";
            system("gzip -dc $df 2> $df.gzip.err | $seqrequester summarize - > $partname 2> $df.seqrequester.err");
        } else {
            printf "  SUMMARIZE - seqrequester > $partname\n";
            system("$seqrequester summarize $df > $partname 2> $df.seqrequester.err");
        }
        #system("rm -f $df");
    }

    loadSummaryFromFile($filetype, $filesize, $filename, $sumsize, $partname, $errors);
}



#
#  Given a species_name, a type-x-individual and a list of tiFiles
#  ("type:indiv filesize datafile\0"), compute a scalaing factor that will
#  convert a file size into an estimated number of bases.
#
#  The estimate is two-fold:
#    1) Three files are paritally downloaded and converted to bases.
#    2) Those three estimates are used to compute a global estimate for
#       all files of the same type:individual.
#
#  Outputs:
#    $data{'data_${ti}_scale'}     -- scale factor for this ti
#    $data{'data_${ti}_bases'}     -- estimated bases for all files with the same ti
#

sub estimateRawDataScaling ($$$$$$$) {
    my $data     = shift @_;   #  inout:   Global %data
    my $ti       = shift @_;   #  input:   "type:indiv" we're estimating for
    my $filelist = shift @_;   #  input:   list of "type:indiv filesize datafile\0" to estimate bases for
    my $filesize =       0;    #  compute: total size of all files with this ti
    my $errors   = shift @_;   #  output:  array of errors encountered
    my $missing  = shift @_;   #  output:  array of files with missing data
    my $download = shift @_;   #  input:   true if we're allowed to download missing data files
    my $downsize = shift @_;   #  input:   limit on download size

    return   if (!defined($filelist));
    return   if ($filelist eq "");

    my @filelist = split '\0', $filelist;
    my $filesLen = scalar(@filelist);

    return   if ($filesLen == 0);

    print " - $ti:\n";

    #  Bionano is special since 'bases' is the coverage of the molecules, not
    #  actual bases in the data file.

    if ($ti =~ m/bionano/) {
        my $basestotal = 0;
        my $bytestotal = 0;

        for (my $ii=0; $ii<$filesLen; $ii++) {
            my ($tind, $size, $file) = split '\s+', $filelist[$ii];

            my $name = $file;         #  We key off the basename for some reason.
            $name =~ s/\.bnx\.gz$//;
            $name =~ s/\.bnx$//;

            #  If genomic_data.summary already knows the bases for this file, use that.
            if (exists($fnBases{$name})) {
                printf "    %s%s %12d %s\n", "S", "E", $size, $file;
            }

            #  Otherwise, download and compute, then update the data for genomic_data.summary.
            else {
                printf "    %s%s %12d %s\n", "S", " ", $size, $file;

                my ($bases, $reads) = computeBionanoBases($data, $tind, $size, $file, $missing, $download);

                $fnBytes{$name} = $size;
                $fnQuant{$name} = $size;
                $fnBases{$name} = $bases;  #sLength;
                $fnReads{$name} = $reads;  #nMolecules;
                $fnTypeI{$name} = $tind;
            }

            $basestotal += $fnBases{$name};
            $bytestotal += $fnBytes{$name};
        }

        $$data{"data_${ti}_scale"} = int(10000 * $basestotal / $bytestotal) / 10000.0;
        $$data{"data_${ti}_bases"} = $basestotal;

        print "\n";

        return;
    }

    #  Sort files by size.

    @filelist = sort { 
        my ($i1, $s1, $n1) = split '\s+', $a;
        my ($i2, $s2, $n2) = split '\s+', $b;
        return($s1 <=> $s2);
    } @filelist;

    #  Use the 1/4, 1/2 and 3/4 largest files to estimate a byte:base ratio.

    my $f1 = int(1 * $filesLen / 4);   #  Pick three files representative file (indices)
    my $f2 = int(2 * $filesLen / 4);   #  from the list of input files sorted by size.
    my $f3 = int(3 * $filesLen / 4);

    my ($ind1, $size1, $file1) = split '\s', $filelist[$f1];   #  Extract name and size from input list;
    my ($ind2, $size2, $file2) = split '\s', $filelist[$f2];   #  bases and seqs are undefined.
    my ($ind3, $size3, $file3) = split '\s', $filelist[$f3];

    #  Print a list of the data present, marking which ones we have
    #  summarized already or will use for estimation.

    for (my $ii=0; $ii<$filesLen; $ii++) {
        my ($ind, $size, $file) = split '\s+', $filelist[$ii];

        my $se = (exists($fnBases{$file})) ? "S" : " ";    #  Summary exists?
        my $us = (($ii == $f1) ||                          #  Used for estimation?
                  ($ii == $f2) ||
                  ($ii == $f3)) ? "E" : " ";

        printf "    %s%s %12d %s\n", $se, $us, $size, $file;

        $filesize += $size;
    }
    printf "\n";

    #  If needed, download a subset of the data from those three files and
    #  use it to estimate the total number of bases in the file.

    downloadPartAndSummarize($ind1, $size1, $file1, $download, $downsize, $errors);
    downloadPartAndSummarize($ind2, $size2, $file2, $download, $downsize, $errors);
    downloadPartAndSummarize($ind3, $size3, $file3, $download, $downsize, $errors);

    #  Get an estimate of (or, rarely, the actual) number of bases in each file.
    #
    #  The scaling has precision limited to x.xxxx to avoid differences in rounding.
    #
    #  If any of the estimates do not exist, fall back to a hard coded estimate,
    #  and emit a warning.

    my $scaling;

    if (!exists($fnBases{$file1}) ||
        !exists($fnBases{$file2}) ||
        !exists($fnBases{$file3})) {
        $$missing{ $$data{"name_"} } .= "$ti ";

        my $errstr = "  WARNING: Size estimates not generated for:  (will use approximate scaling estimate)\n";
        $errstr   .= "    File 1 $file1\n"   if (!exists($fnBases{$file1}));
        $errstr   .= "    File 2 $file2\n"   if (!exists($fnBases{$file2}));
        $errstr   .= "    File 3 $file3\n"   if (!exists($fnBases{$file3}));

        push @$errors, $errstr;

        $scaling = 1.5   if ($file1 =~ m!genomic_data/10x!);           #  Bimodal, ~1.4 and ~1.8
        $scaling = 1.5   if ($file1 =~ m!genomic_data/arima!);         #  Bimodal, ~1.5 and ~1.8
        $scaling = 1.5   if ($file1 =~ m!genomic_data/dovetail!);
        $scaling = 1.8   if ($file1 =~ m!genomic_data/illumina!);
        $scaling = 0.5   if ($file1 =~ m!genomic_data/ont!);           #  NOT TRUE!
        $scaling = 1.1   if ($file1 =~ m!genomic_data/ont_duplex!);    #  NOT TRUE!
        $scaling = 0.5   if ($file1 =~ m!genomic_data/pacbio!);
        $scaling = 1.1   if ($file1 =~ m!genomic_data/pacbio_hifi!);
        $scaling = 1.6   if ($file1 =~ m!genomic_data/phase!);
    }
    else {
        my $sumbases = $fnBases{$file1} + $fnBases{$file2} + $fnBases{$file3};
        my $sumsizes = $size1           + $size2           + $size3;

        $scaling = int(10000 * $sumbases / $sumsizes) / 10000.0;
    }

    #  Save the scaling so we can write it in writeSummaries.

    $$data{"data_${ti}_scale"} =     $scaling;
    $$data{"data_${ti}_bases"} = int($scaling * $filesize);
}

return(1);
