package GenomeArkAssembly;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(rememberLatestAssembly summarizeAssembly importAssemblySummary);

use strict;
use warnings;

use GenomeArkUtility;

use Time::Local;
#use List::Util;
use File::Basename;

my $aws          = "aws";
my $seqrequester = "./seqrequester/build/bin/seqrequester";

#
#  Decode an assembly name to determine
#    is mito?  is genomic?
#    species name
#    assembly label
#    species tag
#    individual number
#    primary/alternate/maternal/paternal
#    date
#

sub parseAssemblyName ($$$) {
    my $filename = shift @_;
    my $filesecs = shift @_;
    my $verbose  = shift @_;

    #  Initialize the return values.

    my ($sName, $aLabel, $sTag, $sNum, $prialt, $date, $curated, $err, $warn) = (undef, "", "", "", "", "", "uncurated", undef, undef);

    #  Handle mito first, because the second form also matches the generic 'an assembly' regex.

    if    ($filename =~ m!species/(.*)/.*/(.*assembly.+)/(.......)(\d)\.MT\.(\d\d\d\d)(\d\d)(\d\d).fasta.gz$!i) {
        print "\n"                               if ($verbose);
        print "$filename\n"                      if ($verbose);
        print " - A mitochondrial assembly!\n"   if ($verbose);

        ($sName, $aLabel, $sTag, $sNum, $prialt, $date) = ($1, $2, $3, $4, "mito", "$5-$6-$7");
    }

    elsif ($filename =~ m!species/(.*)/.*/(.*assembly.+)/(.......)(\d)[\._](.*)\.(\d\d\d\d)(\d\d)(\d\d)\.MT.fasta.gz$!i) {
        print "\n"                               if ($verbose);
        print "$filename\n"                      if ($verbose);
        print " - A mitochondrial assembly!\n"   if ($verbose);

        ($sName, $aLabel, $sTag, $sNum, $prialt, $date) = ($1, $2, $3, $4, "mito", "$6-$7-$8");
    }

    elsif ($filename =~ m!species/(.*)/.*/(.*assembly.+)/(.......)(\d)[\._]\w+\.[WXYZ]\.\w+\.(\d\d\d\d)(\d\d)(\d\d).fasta.gz$!i) {
        print "\n"                             if ($verbose);
        print "$filename\n"                    if ($verbose);
        print " - A merged trio assembly!\n"   if ($verbose);

        ($sName, $aLabel, $sTag, $sNum, $prialt, $date) = ($1, $2, $3, $4, "mgd", "$5-$6-$7");
    }

    elsif ($filename =~ m!species/(.*)/.*/(.*assembly.+)/(.......)(\d)[\._](.*)\.(\d\d\d\d)(\d\d)(\d\d).fasta.gz$!i) {
        print "\n"                     if ($verbose);
        print "$filename\n"            if ($verbose);
        print " - An assembly!   \n"   if ($verbose);

        #  To handle the correct "mat.cur" and "mat.asm" and incorrect "pri"
        #  we extract the whole 'middle word' and split out the dotted parts.
        #    aDisPic1.pri.20210803.fasta.gz
        #
        #  We're expecting the first word to be 'pri', 'alt', 'mat' or 'pat'.
        #  If we see those in any other word, throw an error.
        #    mBosTau1/assembly_curated/mBosTau1.mat.alt.cur.20200213.fasta.gz
        #    mBosTau1/assembly_curated/mBosTau1.pat.alt.cur.20200213.fasta.gz

        my @p = split '\.', $5;

        #if (!defined($p[0])) { die "No first word in $5 filename '$filename'\n"; }
        #if (!defined($p[1])) { die "No second word in $5 filename '$filename'\n"; }

        if ($p[0] eq "HiC") {   #  Adjust for some intermediate assemblies:
            $p[0] = $p[1];      #    Ara_ararauna/bAraAra1/assembly_vgp_HiC_2.0/bAraAra1.HiC.hap1.20220601.fasta.gz
            $p[1] = "asm";      #    Colius_striatus/bColStr4/assembly_vgp_HiC_2.0/bColStr4.HiC.hap1.20220601.fasta.gz
        }                       #

        if (($p[0] ne "pri")  && ($p[0] ne "alt")  &&
            ($p[0] ne "mat")  && ($p[0] ne "pat")  &&
            ($p[0] ne "hap1") && ($p[0] ne "hap2") &&
            ($p[0] ne "dip")) {
            print " - Filename has '$p[0]' as first word; required 'pri/alt', 'mat/pat', 'hap1/hap2' or 'dip'.\n"   if ($verbose);
            $err = "  Filename '$3$4.$5.$6$7$8.fasta.gz' has '$p[0]' as first word; required 'pri/alt', 'mat/pat', 'hap1/hap2' or 'dip'.\n";
        }

        $p[0] = "hpa"  if ($p[0] eq "hap1");   #  Rename hap1/hap2 to hpa/hpb because everything else is three letters
        $p[0] = "hpb"  if ($p[0] eq "hap2");   #  and hap11 looks weird (where last digit is the individual number).

        if ((defined($p[1])) &&
            ($p[1] ne "asm") &&
            ($p[1] ne "cur") &&
            ($p[1] ne "decon")) {
            print " - Filename has '$p[1]' as second word; required 'cur' or 'asm'.\n"   if ($verbose);
            $err = "  Filename '$3$4.$5.$6$7$8.fasta.gz' has '$p[1]' as second word; required 'cur' or 'asm'.\n";
        }

        ($sName, $aLabel, $sTag, $sNum, $prialt, $date) = ($1, $2, $3, $4, $p[0], "$6-$7-$8");
    }

    else {
        print "\n"                               if ($verbose);
        print "$filename\n"                      if ($verbose);
        print " - Failed to parse filename.\n"   if ($verbose);
        $err = "  Failed to parse '$filename'.\n";
    }

    #  Set the curation/deconamination status.

    if ($filename =~ m/decon/) {
        print " - Decontaminated assembly!\n"   if ($verbose);
        $curated = "decontaminated";
    }
    elsif ($aLabel eq "assembly_curated") {
        print " - Curated assembly!\n"   if ($verbose);
        $curated = "curated";
    }
    else {
        $curated = "uncurated";
    }

    #  Convert the 'date' to seconds (more or less from the epoch).  We used
    #  to return 'filesecs', but some files are uploaded well after they're
    #  generated, even after the file that replaces them.

    my $timesecs = $filesecs;

    if ($date =~ m!(\d\d\d\d)-(\d\d)-(\d\d)$!) {
        my ($y, $m, $d);

        if ($2 < 13) {   #  Normal case, $2 is month,
            ($y, $m, $d) =  ($1, $2-1, $3);
        } else {
            print " - Suspected month/day swap\n"  if ($verbose);
            $err = "  Suspected month/day swap in '$filename'.\n";
            ($y, $m, $d) =  ($1, $3-1, $2);
        }

        $timesecs = timelocal(0, 0, 0, $d, $m, $y);
    }

    #  Return a bunch of stuff.

    return($sName, $aLabel, $sTag, $sNum, $prialt, $date, $timesecs, $curated, $err, $warn);
}



#
#  Reads a sequrequester summary file and returns the n50
#
#

sub loadAssemblySummary ($$$$$$) {
    my $summaryfile = shift @_;
    my $ctgNG       = shift @_;
    my $ctgLG       = shift @_;
    my $ctgLEN      = shift @_;
    my $ctgCOV      = shift @_;
    my $genomeSize  = shift @_;
    my $n50         = 0;
    my $size        = 0;

    if (! -e "$summaryfile") {
        push @$ctgNG,  0;
        push @$ctgLG,  0;
        push @$ctgLEN, 0;
        push @$ctgCOV, 0;

        return(0, 0);
    }

    print "Import '$summaryfile'\n";

    open(SU, "< $summaryfile") or die;
    while (<SU>) {
        my ($ng, $lg, $len, $cov);

        #  Match an NG/LG line:
        #    G=16722                            sum of  ||               length     num
        #    NG         length     index       lengths  ||                range    seqs
        #    ----- ------------ --------- ------------  ||  ------------------- -------
        #    00010        16722         0        16722  ||      16722                 1|-----...
        #    00020        16722         0        16722  ||
        #
        if    ($_ =~ m!0*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\|\|!) {
            $ng  = $1;
            $len = $2;
            $lg  = $3;
            $cov = "-";   #  If not valid, otherwise set below.
            $cov = sprintf("%.2f", $4 / $genomeSize)   if (($4 ne "-") && ($genomeSize > 0));

            $n50 = $len   if ($ng == 50);
        }

        #  Match an empty NG/LG line:
        #    G=16722                            sum of  ||               length     num
        #    NG         length     index       lengths  ||                range    seqs
        #    ----- ------------ --------- ------------  ||  ------------------- -------
        #    ...
        #    00070            -         -            -  || 
        #
        elsif ($_ =~ m!0*(\d+)\s+-\s+-\s+-\s+\|\|!) {
            $ng  = $1;
            $len =  0;
            $lg  =  0;
            $cov = "-";
        }

        #  Match a coverage line:
        #    001.000x                   1        16722  ||
        #
        elsif ($_ =~ m!\d+(\d.\d+x)\s+(\d+)\s+(\d+)\s+\|\|!) {
            $ng   = $1;    #  Actual coverage of genome.
            $lg   = $2;    #  Number of sequences.
            $len  = $3;    #  Actual size of assembly.
            $cov  = $1;    #  Actual size of assembly.
            $size = $3;    #  Actual size of assembly.
        }

        #  Not a line we care about.
        else {
            next;
        }

        push @$ctgNG,                $ng;
        push @$ctgLG,                $lg;
        push @$ctgLEN, prettifyBases($len);
        push @$ctgCOV,               $cov;
    }
    close(SU);

    return($n50, $size);   #  Return the n50 and total size of assembly.
}



sub generateAssemblySummaryHTML ($$$$) {
    my $prialt     = shift @_;
    my $filename   = shift @_;
    my $filesize   = shift @_;
    my $genomeSize = shift @_;
    my $n50table   = "";

    $filename =~ s/.gz//;

    $n50table .= "|\n";
    $n50table .= "  <table class=\"sequence-sizes-table\">\n";
    $n50table .= "  <thead>\n";

    $n50table .= "  <tr>\n";
    $n50table .= "  <th></th>\n";
    $n50table .= "  <th colspan=2 align=center>Contigs</th>\n";
    $n50table .= "  <th colspan=2 align=center>Scaffolds</th>\n";
    $n50table .= "  </tr>\n";

    $n50table .= "  <tr>\n";
    $n50table .= "  <th>NG</th>\n";
    $n50table .= "  <th>LG</th>\n";
    $n50table .= "  <th>Len</th>\n";
    $n50table .= "  <th>LG</th>\n";
    $n50table .= "  <th>Len</th>\n";
    $n50table .= "  </tr>\n";
    $n50table .= "  </thead>\n";
    $n50table .= "  <tbody>\n";

    #  Strip out the N50 chart

    my (@ctgNG, @ctgLG, @ctgLEN, @ctgCOV);
    my (@scfNG, @scfLG, @scfLEN, @scfCOV);

    my ($ctgn50, $ctgsize) = loadAssemblySummary("$filename.ctg.summary", \@ctgNG, \@ctgLG, \@ctgLEN, \@ctgCOV, $genomeSize);
    my ($scfn50, $scfsize) = loadAssemblySummary("$filename.scf.summary", \@scfNG, \@scfLG, \@scfLEN, \@scfCOV, $genomeSize);

    #  Ten rows of actual data.

    for (my $ii=0; $ii<10; $ii++) {
        my $ctgcolor = "";
        my $scfcolor = "";

        #  For non-alternate assemblies (primary, maternal, paternal)
        #  highlight the n-50 line if it is crappy.

        if (($ii == 4) && ($prialt ne "alt")) {
            $ctgcolor = ($ctgn50 <   1000000) ? " style=\"background-color:#ff8888;\"" : " style=\"background-color:#88ff88;\"";
            $scfcolor = ($scfn50 <  10000000) ? " style=\"background-color:#ff8888;\"" : " style=\"background-color:#88ff88;\"";
        }

        if ($ii == 4) {
            $n50table .= "  <tr style=\"background-color:#cccccc;\">";
        } else {
            $n50table .= "  <tr>";
        }

        $n50table .= "<td> $ctgNG[$ii] </td>";
        $n50table .= "<td> $ctgLG[$ii] </td><td$ctgcolor> $ctgLEN[$ii] </td>";  #<td> $ctgCOV[$ii] </td>
        $n50table .= "<td> $scfLG[$ii] </td><td$scfcolor> $scfLEN[$ii] </td>";  #<td> $scfCOV[$ii] </td>
        $n50table .= "</tr>";
    }

    #  And one row of summary.

    $n50table .= "  </tbody>\n";
    $n50table .= "  <tfoot>\n";
    $n50table .= "  <tr>";
    $n50table .= "<th> $ctgNG[10] </th>";
    $n50table .= "<th> $ctgLG[10] </th><th> $ctgLEN[10] </th>";
    $n50table .= "<th> $scfLG[10] </th><th> $scfLEN[10] </th>";
    $n50table .= "</tr>\n";
    $n50table .= "  </tfoot>\n";

    $n50table .= "  </table>\n";

    return($ctgn50, $scfn50, $scfsize, $n50table);
}






#  On the first pass, we're just trying to find the assembly with the most
#  recent date.
#
#  To handle assemblies with multiple curations, we need to track both:
#    Is curated assembly available?
#    Which assembly to use?
#
sub rememberLatestAssembly ($$$$$) {
    my $filesecs = shift @_;
    my $filesize = shift @_;
    my $filename = shift @_;
    my $filebase = $filename;   $filebase =~ s/.gz$//;
    my $data     = shift @_;
    my $errors   = shift @_;

    my ($sName, $aLabel, $sTag, $sNum, $prialt, $date, $secs, $curated, $err, $warn) = parseAssemblyName($filename, $filesecs, 1);

    if (defined($err)) {       #  Oops, something wrong.
        push @$errors, $err;   #  Append the error to the global list of errors,
        return;                #  and bail.
    }

    if (defined($warn)) {      #  Oops, something amss.
        push @$errors, $warn;  #  Append the error to the global list of errors,
    }

    #  If there isn't a date sete for this prialt, set it to whatever we have.

    if (!exists($$data{"${prialt}${sNum}__datesecs"})) {
        print " - Set ${prialt}${sNum} to $date ($curated)\n";
        $$data{"${prialt}${sNum}__datesecs"} = $secs;
        $$data{"${prialt}${sNum}__curated"}  = $curated;
        $$data{"${prialt}${sNum}__filename"} = $filename;
    }

    #  If we've seen a curated assembly, and this one isn't, skip it.

    if (($$data{"${prialt}${sNum}__curated"} eq "curated") && ($curated ne "curated")) {
        print " - Ignore $curated assembly\n";
        return;
    }

    #  If this assembly is curated, and it's the first, save it regardless.

    if (($$data{"${prialt}${sNum}__curated"} ne "curated") && ($curated eq "curated")) {
        print " - Set ${prialt}${sNum} to curated assembly $date\n";
        $$data{"${prialt}${sNum}date"}       = $date;
        $$data{"${prialt}${sNum}__datesecs"} = $secs;
        $$data{"${prialt}${sNum}__curated"}  = $curated;
        $$data{"${prialt}${sNum}__filename"} = $filename;
    }

    #  If the prialt status is curated, but this assembly isn't, something is horribly wrong.

    if (($$data{"${prialt}${sNum}__curated"} eq "curated") &&
        ($curated ne "curated")) {
        die "${prialt}${sNum} status set to curated, but file $filename is $curated.\n";
    }

    #  Pick the later date.

    if ($secs < $$data{"${prialt}${sNum}__datesecs"}) {
        print " - Obsolete\n";
        return;
    }

    #  If the date of this file is newer than the date of our data, reset the date of our data
    #  and wipe out all the old data.

    if ($$data{"${prialt}${sNum}__datesecs"} < $secs) {
        print " - Update ${prialt}${sNum} to $date ", ($curated) ? "($curated)" : "", "\n";
        $$data{"${prialt}${sNum}__datesecs"} = $secs;
        $$data{"${prialt}${sNum}__curated"}  = $curated;
        $$data{"${prialt}${sNum}__filename"} = $filename;
    }
}



#
#  Downloads the assembly and generates a summary.  Does nothing if the
#  summary already exists.
#
#  Reports errors if either the download or summarizing fails.  In this case,
#  no output summary file is generated.
#

sub generateAssemblySummary ($$$$$$) {
    my $filename   = shift @_;
    my $filebase   = $filename;   $filebase =~ s/.gz$//;
    my $dirname    = dirname($filename);
    my $filesize   = shift @_;
    my $type       = shift @_;
    my $genomeSize = shift @_;
    my $errors     = shift @_;
    my $download   = shift @_;

    $genomeSize = 0   if ((!defined($genomeSize)) || ($genomeSize eq ""));

    my $splitopt   = ($type eq "ctg")  ? "-split-n"          : "";
    my $sizeopt    = ($genomeSize > 0) ? "-size $genomeSize" : "";

    my $awsopts    = "--only-show-errors --no-sign-request";  #--quiet --no-progress

    #  If the summary exists and is not empty, we have nothing to do here.

    if ((! -z "$filebase.$type.summary") &&
        (  -e "$filebase.$type.summary")) {
        return;
    }

    #  If we're allowed to download, and we need to download, download.
    #  Generate a warning if the file doesn't exist after 'downloading' it.

    if (($download == 1) && (! -e "downloads/$filebase.gz")) {
        printf "              FETCH asm - size %6.3f GB\n", $filesize / 1024 / 1024 / 1024;

        system("mkdir -p downloads/$dirname");
        system("aws $awsopts s3 cp s3://genomeark/$filebase.gz downloads/$filebase.gz || rm -f downloads/$filebase.gz");
    }
    if (! -e "downloads/$filebase.gz") {
        #ush @$errors, "  Downloading disabled for file $filebase.gz\n"   if ($download == 0);
        push @$errors, "  Failed to download $filebase.gz\n"              if ($download == 1);
        return;
    }

    #  Summarize it!

    print "              SUMMARIZING as $type\n"                                if ($genomeSize == 0);
    print "              SUMMARIZING as $type with genome size $genomeSize\n"   if ($genomeSize  > 0);

    system("mkdir -p $dirname");
    system("$seqrequester summarize -1x $splitopt $sizeopt downloads/$filebase.gz > $filebase.$type.summary || rm -f $filebase.$type.summary");

    #  Remove incomplete summaries, and generate a warning.

    unlink "$filebase.scf.summary"                          if (  -z "$filebase.$type.summary");
    push @$errors, "  Failed to summarize $filebase.gz\n"   if (! -e "$filebase.$type.summary");
}



#
#  On the second pass, we process the assembly with the most recent date.
#
#  But if there are multiple assemblies with the same date, we get a little confused
#  and show only the last one encountered.
#

sub summarizeAssembly ($$$$$$$) {
    my $filesecs = shift @_;
    my $filesize = shift @_;
    my $filename = shift @_;
    my $filebase = $filename;   $filebase =~ s/.gz$//;
    my $data     = shift @_;
    my $errors   = shift @_;
    my $missing  = shift @_;
    my $download = shift @_;

    my ($sName, $aLabel, $sTag, $sNum, $prialt, $date, $secs, $curated, $err, $warn) = parseAssemblyName($filename, $filesecs, 0);

    if (defined($err)) {       #  Oops, something wrong.
        return;                #  But we've already saved the error, so just bail.
    }

    if (defined($warn)) {      #  Oops, something amiss.
    }                          #  But we've already saved the error, so do nothing.

    #return  if ($secs     ne $$data{"${prialt}${sNum}__datesecs"});
    return  if (!exists($$data{"${prialt}${sNum}__filename"}));       #  Already processed the correct assembly.
    return  if ($filename ne $$data{"${prialt}${sNum}__filename"});   #  Not the correct assembly.

    delete $$data{"${prialt}${sNum}__datesecs"};   #  We're done with these; they were
    delete $$data{"${prialt}${sNum}__curated"};    #  only used to decide which $filename
    delete $$data{"${prialt}${sNum}__filename"};   #  should be processed here.

    #  If more than 9 individuals:
    #    _layouts/genomeark.html needs to have more assembly entries added.
    #    utilities/make-*-table.pl needs to have more individuals added.
    #
    if ($sNum > 9) {
        die "Too many assemblies; update _layouts/genomeark.html and utilities/make-*-table.pl.\n";
    }

    #  Log that we're using this assembly.

    printf "  %8s <- %s\n", "$prialt$sNum", $filename;

    #  And generate summaries of the contigs and scaffolds.

    generateAssemblySummary($filename, $filesize, "ctg", $$data{"genome_size"}, $errors, $download)   if (! -e "$filebase.ctg.summary");
    generateAssemblySummary($filename, $filesize, "scf", $$data{"genome_size"}, $errors, $download)   if (! -e "$filebase.scf.summary");

    my (@ctgNG, @ctgLG, @ctgLEN, @ctgCOV);
    my (@scfNG, @scfLG, @scfLEN, @scfCOV);

    my ($ctgn50, $ctgsize) = loadAssemblySummary("$filename.ctg.summary", \@ctgNG, \@ctgLG, \@ctgLEN, \@ctgCOV, $$data{"genome_size"});
    my ($scfn50, $scfsize) = loadAssemblySummary("$filename.scf.summary", \@scfNG, \@scfLG, \@scfLEN, \@scfCOV, $$data{"genome_size"});

    #unlink "$filebase.ctg.summary"   if (-z "$filebase.$type.summary");

    if ((! -e "$filebase.ctg.summary") ||
        (! -e "$filebase.scf.summary")) {
        $$missing{ $$data{"name_"} } = 1;

        if ($download == 0) {
            push @$errors, "  Size stats not generated for '$filebase'\n";
        }
        else {
            push @$errors, "  FAILED to generate size stats for '$filebase'\n"   
        }
    }
}



sub importAssemblySummary ($$$$$$) {
    my $filesecs = shift @_;
    my $filesize = shift @_;
    my $filename = shift @_;
    my $filebase = $filename;   $filebase =~ s/.gz$//;
    my $data     = shift @_;
    my $errors   = shift @_;
    my $missing  = shift @_;

    my ($sName, $aLabel, $sTag, $sNum, $prialt, $date, $secs, $curated, $err, $warn) = parseAssemblyName($filename, $filesecs, 0);

    if (defined($err)) {       #  Oops, something wrong.
        return;                #  But we've already saved the error, so just bail.
    }

    if (defined($warn)) {      #  Oops, something amiss.
    }                          #  But we've already saved the error, so do nothing.

    #return  if ($secs     ne $$data{"${prialt}${sNum}__datesecs"});
    return  if (!exists($$data{"${prialt}${sNum}__filename"}));       #  Already processed the correct assembly.
    return  if ($filename ne $$data{"${prialt}${sNum}__filename"});   #  Not the correct assembly.

    delete $$data{"${prialt}${sNum}__datesecs"};   #  We're done with these; they were
    delete $$data{"${prialt}${sNum}__curated"};    #  only used to decide which $filename
    delete $$data{"${prialt}${sNum}__filename"};   #  should be processed here.

    # Set date, assembly label, pretty-print file size, path to file, etc, etc.

    $$data{"${prialt}${sNum}date"}     = $date;
    $$data{"${prialt}${sNum}version"}  = $aLabel;

    $$data{"${prialt}${sNum}filesize"} = prettifySize($filesize);

    $$data{"${prialt}${sNum}seq"}      = "https://s3.amazonaws.com/genomeark/$filename";

    ($$data{"${prialt}${sNum}n50ctg"},
     $$data{"${prialt}${sNum}n50scf"},
     $$data{"${prialt}${sNum}length"},
     $$data{"${prialt}${sNum}sizes"})  = generateAssemblySummaryHTML($prialt, $filename, $filesize, $$data{"genome_size"});

    #  Update the assembly status.  This used to be based on n50 of the
    #  primary, but that isn't useful anymore (there were no low-quality
    #  draft assemblies with n50 < 1 Mbp) and so was removed.  See e.g.,
    #  19cbac3d3ff8dc3564162d6bf24e64ed645600f3.

    if (($prialt eq "pri") ||
        ($prialt eq "mat") || ($prialt eq "pat") ||
        ($prialt eq "hpa") || ($prialt eq "hpb") ||
        ($prialt eq "mgd") ||
        ($prialt eq "dip")) {
        my $st;

        $st = "draft";                                        #  Decide if this assembly is
        $st = "curated"  if ($curated eq "curated");          #  a draft (by default) or
        $st = "curated"  if ($curated eq "decontaminated");   #  curated.
        $st = "curated"  if ($prialt eq "mgd");

        if ($$data{"assembly_status"} ne "curated") {   #  If the assembly status isn't curated,
            $$data{"assembly_status"} = $st;            #  set it to whatever.
        }

        printf "  %8s (%s) -> species (%s)\n", "$prialt$sNum", $st, $$data{"assembly_status"};
    }
}


return(1);
