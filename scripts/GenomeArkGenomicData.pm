package GenomeArkGenomicData;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(@dataTypes loadSummaries writeSummaries estimateRawDataScaling);

use strict;
use warnings;

use GenomeArkUtility;
use GenomeArkBionano;

#use Time::Local;
use List::Util qw(min max);
use File::Basename;

my $aws          = "aws --no-sign-request";
my $seqrequester = "./seqrequester/build/bin/seqrequester";

my %dataBytes;   #  Number of bytes in a genomic_data sequence file.
my %dataQuant;   #  Amount of raw file summarized (in gigabytes, or "all").
my %dataBases;   #  Number of bases in piece summarized.
my %dataReads;   #  Number of reads in piece summarized.


#
#  Update genomeark.github.io/_layouts/genomeark.html if this changes:
#   - the coverage table
#   - the list of links to download
#
our @dataTypes = qw(10x
                    arima
                    bionano
                    dovetail
                    illumina
                    ont
                    ontduplex
                    pacbio
                    pacbiohifi_fqgz
                    pacbiohifi_bam
                    pacbiohifi_clr
                    phase);


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
sub loadSummaryFromFile ($$$$$) {
    my $filesize  = shift @_;   #  Size of the file we've summarized.
    my $filename  = shift @_;   #  Name of the file we've summarized.
    my $sumsize   = shift @_;   #  Amount of the file we summarized.
    my $summary   = shift @_;   #  Name of the summary file.
    my $errors    = shift @_;

    $dataBytes{$filename} = $filesize;
    $dataQuant{$filename} = $sumsize;
    $dataBases{$filename} = 0;
    $dataReads{$filename} = 0;

    print "  Import summary from '$summary'\n";

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

    if (($dataBases{$filename} == 0) ||
        ($dataReads{$filename} == 0)) {
        push @$errors, "Summary FAILED, no bases or reads found.  Summary saved in '$summary.BAD'.\n";

        unlink "$summary.BAD";
        rename "$summary", "$summary.BAD";

        delete $dataBytes{$filename};
        delete $dataQuant{$filename};
        delete $dataBases{$filename};
        delete $dataReads{$filename};
    }
}

#
#  Load a summary, if any summary exists.
#
sub loadSummaryIfExists ($$$) {
    my $filesize  = shift @_;
    my $filename  = shift @_;
    my $fullname  = makeSummaryFileName($filename, undef);   #  Name of a full summary.
    my $partnames = makeSummaryFileName($filename, "*");     #  Glob to find any partial summaries.
    my $errors    = shift @_;

    my $sumsize   = (exists($dataQuant{$filename})) ? $dataQuant{$filename} : 0;

    my $newsize   = 0;
    my $newname   = "";

    #print "loadSummaryIfExists()- look for '$fullname' and\n";
    #print "                                '$partnames'\n";

    #  If a full summary exists, always use that.
    if (-e $fullname) {
        #print "loadSummaryIfExists()- found    '$fullname'\n";
        $newsize = $filesize;
        $newname = $fullname;
    }

    #  Otherwise, search for a partial summary of any size.
    else {
        open(LS, "ls $partnames 2>/dev/null |");
        while (<LS>) {
            chomp;
            #print "loadSummaryIfExists()- found    '$_'\n";
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
        #print "loadSummaryIfExists()- LOAD     '$fullname'\n";
        loadSummaryFromFile($filesize, $filename, $newsize, $newname, $errors);
    }
}


#
#  Load existing data summaries from 'species/$name/genomic_data.summary',
#  then scan the list of seqFiles to find any new summaries that aren't in
#  the global file - or any summaries that improve over what we know.
#
sub loadSummaries ($$$) {
    my $data     = shift @_;
    my $name     = $$data{"name_"};
    my $seqFiles = shift @_;
    my $errors   = shift @_;

    undef %dataBytes;
    undef %dataQuant;
    undef %dataBases;
    undef %dataReads;

    my $bionanoBases = 0;

    #  Load precomputed summaries.

    if (-e "species/$name/genomic_data.summary") {
        print "  Load summaries from species/$name/genomic_data.summary\n";

        open(ST, "< species/$name/genomic_data.summary") or die;
        while (<ST>) {
            s/^\s+//;
            s/\s+$//;

            next   if (m/^-/);
            next   if (m/^$/);

            my @v = split '\s+', $_;

            #   legacy
            next   if (scalar(@v) == 0);
            next   if ($v[0] eq "bytes");

            if ((scalar(@v) == 6) && ($v[3] > 0)) {
                my ($bytes, $date, $quant, $bases, $reads, $file) = @v;

                $dataBytes{$file} = $bytes;
                $dataQuant{$file} = $quant;
                $dataBases{$file} = $bases;
                $dataReads{$file} = $reads;

                $bionanoBases += $bases   if ($file =~ m!genomic_data/bionano/!);
            }
            if ((scalar(@v == 3)) && ($v[2] > 0)) {
                my ($type, $ext, $scale) = @v;

                $$data{"data_${type}_scale"} = $scale;
            }
        }
        close(ST);

        if ($bionanoBases > 0) {
            $$data{"data_bionano_bases"} = $bionanoBases;
        }
    }

    #  Load any summaries if they are new or better.

    foreach my $type (sort keys %$seqFiles) {
        my @files = split '\0', $$seqFiles{$type};

        foreach my $sizefile (@files) {
            my ($s, $n) = split '\s+', $sizefile;

            loadSummaryIfExists($s, $n, $errors);
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
    print ST "-         bytes            date        quantity           bases           reads  file\n";
    print ST "--------------- --------------- --------------- --------------- ---------------  ----------\n";

    foreach my $file (sort keys %dataBytes) {
        my  $date = 0;   #  For later use, maybe.

        next   if (! exists($dataBytes{$file}));
        next   if (! exists($dataQuant{$file}));
        next   if (! exists($dataBases{$file}));
        next   if (! exists($dataReads{$file}));

        printf(ST "%15d %15d %15s %15d %15d  %s\n", $dataBytes{$file}, $date, $dataQuant{$file}, $dataBases{$file}, $dataReads{$file}, $file);
    }

    print ST "--------------- --------------- --------------- --------------- ---------------  ----------\n";
    print ST "\n";
    print ST "-  GENOMIC DATA SCALING FACTORS\n";
    print ST "-\n";
    print ST "-   datatype suffix scaling\n";
    print ST "------------ ------ -------\n";

    foreach my $type (@dataTypes) {
        my $scale = $$data{"data_${type}_scale"};

        next   if (! defined($scale));
        next   if ($scale <= 0.001);

        printf(ST "%-12s %-6s %7.4f\n", $type, "*", $$data{"data_${type}_scale"});
    }

    print ST "------------ ------ -------\n";

    close(ST);
}


#
#  If the summary doesn't exist for all the data, download the whole file and summarize.
#
sub downloadFullAndSummarize ($$$$) {
    my $filesize = shift @_;
    my $filename = shift @_;
    my $fdir     = dirname($filename);
    my $download = shift @_;
    my $errors   = shift @_;

    my $fullname = makeSummaryFileName($filename, undef);

    #  If the whole summary exists, do not recompute it.

    return   if (exists($dataQuant{$filename}) && $dataQuant{$filename} eq "all");

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

    #  If a bam, convert to fastq then summarize, otherwise, summarize directly.

    if    ($filename =~ m/bam$/) {
        printf "EXTRACT and SUMMARIZE to $fullname\n";

        system("samtools fasta downloads/$filename | $seqrequester summarize - > $fullname");
    }
    else {
        printf "SUMMARIZE to $fullname\n";

        system("$seqrequester summarize downloads/$filename > $fullname");
    }

    #  Parse the summary to find the number of bases in the dataset.

    loadSummaryFromFile($filesize, $filename, $filesize, $fullname, $errors);
}


#
#  If no estimated number of bases exists, download a part of the file and summarize.
#
sub downloadPartAndSummarize ($$$$$) {
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

    my $oldsize  = (exists($dataQuant{$filename})) ? $dataQuant{$filename} : 0;
    my $sumsize  = $downsize;

    my $fullname = makeSummaryFileName($filename, undef);
    my $partname = makeSummaryFileName($filename, $sumsize);

    #  If we've summarized all the data, we can't improve it.
    if    ($oldsize >= $filesize) {
        return;
    }

    #  If there is a full summary around, use that.
    elsif (-e $fullname) {
        loadSummaryFromFile($filesize, $filename, $filesize, $fullname, $errors);
        return;
    }

    #  If there is ALREADY a data file around to summarize, make a full summary.
    elsif (-e "downloads/$filename") {
        downloadFullAndSummarize($filesize, $filename, $download, $errors);
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
        return;
    }

    #  If the file isn't hugely bigger than the size we'd sample, just do the whole thing.
    if ($filesize < 2 * $sumsize) {
        downloadFullAndSummarize($filesize, $filename, $download, $errors);
        return;
    }

    #  While 's3 cp' is happy writing to stdout, 's3api get-object' is not,
    #  and requires a file name.  Rather annoying, since this is where we
    #  could actually get away with piping it to the summarizer.

    system("mkdir -p downloads/$fdir");   #  Make a place to download the file
    system("mkdir -p           $fdir");   #  and a place to write the summary.

    if ($filename =~ m/bam$/) {
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

    loadSummaryFromFile($filesize, $filename, $sumsize, $partname, $errors);
}



#
#  Given a species_name, a genomic data type and a list of seqFiles
#  ("filesize datafile" separated by \0 bytes), compute a
#  scalaing factor that will convert a filesize into an estimated number of
#  bases in the (compressed) file.
#
sub estimateRawDataScaling ($$$$$$$) {
    my $data     = shift @_;
    my $name     = $$data{"name_"};
    my $type     = shift @_;
    my $files    = shift @_;
    my $errors   = shift @_;
    my $missing  = shift @_;
    my $download = shift @_;
    my $downsize = shift @_;
    my $warning  = 0;  #shift @_;

    return   if (!defined($files) || ($files eq ""));

    my @files    = split '\0', $files;
    my $filesLen = scalar(@files);

    return   if ($filesLen == 0);

    @files = sort { 
        my ($s1, $n1) = split '\s+', $a;
        my ($s2, $n2) = split '\s+', $b;
        return($s1 <=> $s2);
    } @files;

    if ($type eq "bionano") {
        $$data{"data_${type}_bases"} = computeBionanoBases($data, \@files, \%dataBytes, \%dataQuant, \%dataBases, \%dataReads, $missing, $download);
        return;
    }

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

    downloadPartAndSummarize($size1, $file1, $download, $downsize, $errors);
    downloadPartAndSummarize($size2, $file2, $download, $downsize, $errors);
    downloadPartAndSummarize($size3, $file3, $download, $downsize, $errors);

    writeSummaries($data);

    #  Get an estimate of (or, rarely, the actual) number of bases in each file.

    my $scaling;

    if (!exists($dataBases{$file1}) ||        #  Fail if any of the size computations failed.
        !exists($dataBases{$file2}) ||
        !exists($dataBases{$file3})) {
        $$missing{ $$data{"name_"} } .= "$type ";

        if ($warning) {
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
        }

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
        my $sumbases = $dataBases{$file1} + $dataBases{$file2} + $dataBases{$file3};
        my $sumsizes = $size1             + $size2             + $size3;

        #  Limit precision to avoid differences in rounding
        $scaling = int(10000 * $sumbases / $sumsizes) / 10000.0;
    }

    #  Save the scaling so we can write it in writeSummaries.
    $$data{"data_${type}_scale"} = $scaling;
}

return(1);
