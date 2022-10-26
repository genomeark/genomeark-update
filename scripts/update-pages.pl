#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use lib "$FindBin::RealBin";

use Time::Local;
#use List::Util;

use GenomeArkUtility;
use GenomeArkUpdate;
use GenomeArkGenbank;
use GenomeArkProject;
use GenomeArkMetadata;
use GenomeArkAssembly;
use GenomeArkAccumulateData;
use GenomeArkGenomicData;

#use Slack;
#use Update;
#use Schedule;

#
#  Parse options.
#

my @speciesList;
my $updatePages = 1;

foreach my $arg (@ARGV) {
    if    ($arg eq "--no-update") {
        $updatePages = 0;
    }
    else {
        if ($arg =~ m!^species/(.*)/?$!) {   #  Expect to find a list of species to process
            push @speciesList, $1;           #  on the command line.  Allow either
        } else {                             #  'species/NAME' (from e.g., scan-bucket.pl
            push @speciesList, $arg;         #  species/T*) or the actual name.
        }
    }
}

#
#  Main
#

die "ERROR: 'downloads/genomeark.ls' doesn't exist, can't update.\n"   if (! -e "downloads/genomeark.ls");
die "ERROR: 'downloads/species-list' doesn't exist, can't update.\n"   if (! -e "downloads/species-list");

my $ga = "../genomeark.github.io";

@speciesList = discoverSpecies(@speciesList);   #  List of all species in genomeark-metadata or supplied on command line.

my $lastupdate  = loadBucketFiles();        #  Private list of all files in the bucket.
my $nGenBank    = loadGenbankMap();         #  Private map from ToLID to GenBank accession.
my $nProject    = loadProjectMap();         #  Private map from species name_ to project.

my %missingData;                            #  List of species that need to download data.
my @potentialErrors;                        #  List of fatal errors we encountered.
my @unknownFiles;                           #  List of files we don't know how to process.

#  Create output directories, now that we know what projects exist.

system("mkdir -p $ga/_genomeark");

foreach my $proj ("genomeark", getAllProjectNames()) {
    system("mkdir -p $ga/_$proj-all");
    system("mkdir -p $ga/_$proj-curated-assembly");
    system("mkdir -p $ga/_$proj-hqdraft-assembly");
    system("mkdir -p $ga/_$proj-lqdraft-assembly");
    system("mkdir -p $ga/_$proj-raw-data-only");
}

#  Make index pages for each of the categories.

foreach my $proj ("genomeark", getAllProjectNames()) {
    makeIndexPage($proj, "$ga", "$proj-all",              "All Species");
    makeIndexPage($proj, "$ga", "$proj-curated-assembly", "Species with Curated Assemblies");
    makeIndexPage($proj, "$ga", "$proj-hqdraft-assembly", "Species with High-Quality Draft Assembles");
    makeIndexPage($proj, "$ga", "$proj-lqdraft-assembly", "Species with Low-Quality Draft Assembles");
    makeIndexPage($proj, "$ga", "$proj-raw-data-only",    "Species without an Assembly");
}

#  Iterate over all species and do work!

foreach my $species (@speciesList) {
    my @speciesFiles = ();
    my @speciesSizes = ();
    my @speciesEpoch = ();

    #  Rebuild each species.md file.

    my %seqFiles = ();
    my %seqBytes = ();
    my %seqIndiv = ();
    my %meta     = ();
    my %data     = ();

    print "\n";
    print "----------------------------------------\n";
    print "Processing species  $species\n";

    my $name = loadSpeciesMetadata(\%data, $species, \%meta, \@potentialErrors);

    setGenbankIDs(\%data);

    $data{"last_updated"} = getBucketFiles($name, \@speciesFiles, \@speciesSizes, \@speciesEpoch);

    print "  with ", scalar(@speciesFiles), " files.\n";

    print "\n";
    print "----------\n";
    print "Assemblies, pass 1: remember the latest assembly\n";

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        if (isAssemblyFile($filename, 0)) {
            if ($data{"last_updated"} < $filesecs) {
                $data{"last_updated"} = $filesecs;
            }

            rememberLatestAssembly($filesecs, $filesize, $filename, \%data, \@potentialErrors);
        }
    }


    print "\n";
    print "----------\n";
    print "Assemblies, pass 2: import summaries\n";

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        next   if (! isAssemblyFile($filename, 0));

        importAssemblySummary($filesecs, $filesize, $filename, \%data, \@potentialErrors, \%missingData);
    }

    #
    #  Update genome size to one of the assembly sizes, if needed.
    #

    $data{"genome_size"} = $data{"pri1length"}   if ($data{"genome_size"} == 0) && defined($data{"pri1length"});
    $data{"genome_size"} = $data{"pri2length"}   if ($data{"genome_size"} == 0) && defined($data{"pri2length"});
    $data{"genome_size"} = $data{"pri3length"}   if ($data{"genome_size"} == 0) && defined($data{"pri3length"});
    $data{"genome_size"} = $data{"pri4length"}   if ($data{"genome_size"} == 0) && defined($data{"pri4length"});

    $data{"genome_size_display"} = prettifyBases($data{"genome_size"});



    print "\n";
    print "----------\n";
    print "Genomic Data, accumulation\n";

    foreach my $type (@dataTypes) {
        $seqFiles{$type} = "";
        $seqBytes{$type} = 0;
        $seqIndiv{$type} = "";
    }

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        next   if (! isGenomicDataFile($filename));

        accumulateData($filesecs, $filesize, $filename, \%data, \%seqFiles, \%seqBytes, \%seqIndiv, \@potentialErrors);
    }


    print "\n";
    print "----------\n";
    print "Page Creation\n";

    loadSummaries(\%data, \%seqFiles, \@potentialErrors);   #  Load existing summaries.

    #
    #  Figure out how much and what types of data exist.
    #

    sub uniquifyStringArray ($) {    #  Split a \0 delimited string into an
        my @a = split '\0', $_[0];   #  array of unique elements.
        my %a;

        foreach my $a (@a) {
            $a{$a}++;
        }

        return(sort keys %a);
    }

    sub setCoverage ($$) {
        my $b = shift @_;
        my $g = shift @_;

        if ($g > 0) {
            return(sprintf("%.2fx", $b / $g));
        } else {
            return("N/A");
        }
    }


    foreach my $type (@dataTypes) {
        next  if ($seqBytes{$type} == 0);

        #  Argh, the location of the reads is not the same as 'type'.

        my $path = $type;
        $path = "pacbio_hifi"   if ($type eq "pacbiohifi_fqgz");
        $path = "pacbio_hifi"   if ($type eq "pacbiohifi_bam");
        $path = "pacbio_hifi"   if ($type eq "pacbiohifi_clr");
        $path = "ont_duplex"    if ($type eq "ontduplex");

        my $tnam = $type;
        $type = "10x"              if ($tnam eq "10x        ");
        $type = "arima"            if ($tnam eq "Arima      ");
        $type = "bionano"          if ($tnam eq "Bionano    ");
        $type = "dovetail"         if ($tnam eq "Dovetail   ");
        $type = "illumina"         if ($tnam eq "Illumina   ");
        $type = "ont"              if ($tnam eq "ONT Simplex");
        $type = "ontduplex"        if ($tnam eq "ONT Duplex ");
        $type = "pacbio"           if ($tnam eq "PB CLR     ");
        $type = "pacbiohifi_fqgz"  if ($tnam eq "PB HiFi Q20");
        $type = "pacbiohifi_bam"   if ($tnam eq "PB HfFi All");
        $type = "pacbiohifi_clr"   if ($tnam eq "PB HiFi Sub");
        $type = "phase"            if ($tnam eq "Phase      ");

#
#  DON'T list hifi subreads
#

        foreach my $k (uniquifyStringArray($seqIndiv{$type})) {
            $data{"data_${type}_links"} .= "s3://genomeark/species/$k/genomic_data/${path}/<br>";
        }

        if ($type eq "bionano") {                 #  BioNano has no scaling factor,
            my $b = $data{"data_${type}_bases"};  #  but it knows the bases.

            if (defined($b)) {
                $data{"data_${type}_bytes"}    = sprintf("%.3f GB", $seqBytes{$type} / 1024 / 1024 / 1024);
                $data{"data_${type}_coverage"} = setCoverage($b, $data{"genome_size"});
                $data{"data_${type}_bases"}    = prettifyBases($b);
                $data{"data_${type}_files"}    = 666;
            }
        }
        else {
            next  if (!defined($data{"data_${type}_scale"}));
            next  if ($data{"data_${type}_scale"} < 0.001);
            next  if ($data{"genome_size"} == 0);

            $data{"data_${type}_bytes"}    = sprintf("%.3f GB", $seqBytes{$type} / 1024 / 1024 / 1024);
            $data{"data_${type}_coverage"} = setCoverage($seqBytes{$type} * $data{"data_${type}_scale"}, $data{"genome_size"});
            $data{"data_${type}_bases"}    = prettifyBases($seqBytes{$type} * $data{"data_${type}_scale"});
            $data{"data_${type}_files"}    = 666; #$seqFiles{$type};
        }
    }

    #
    #  Set the classification.
    #    If no assembly, put it under either "some data" or "all data".
    #    Otherwise, put it under the assembly classification.
    #
    #  Add classifications for:
    #    HiFi + ONT
    #    Pac + HIC/BioNano
    #
    #  HIC = arima or dovetail or phase
    #
    #  DIFFERENT from what it was:
    #    all == pac and 10x and HIC and bionano
    #

    sub hasCoverage ($@) {
        my $data = shift @_;
        foreach my $type (@_) {
            next   if (!defined($$data{"data_${type}_coverage"}));
            next   if ($$data{"data_${type}_coverage"} eq "N/A");
            return(1);
        }
        return(0);
    }

    my $hasHQ = hasCoverage(\%data, "pacbiohifi_fqgz", "pacbiohifi_bam", "ontduplex");             #  Has High Quality Long Reads
    my $hasLR = hasCoverage(\%data, "pacbio", "ont");                                              #  Has Long Reads
    my $hasSR = hasCoverage(\%data, "10x", "arima", "bionano", "dovetail", "illumina", "phase");   #  Has Short Reads
    my $hasPH = hasCoverage(\%data, "10x", "bionano");                                             #  Has Phasing Reads
    my $hasSC = hasCoverage(\%data, "arima", "dovetail", "phase");                                 #  Has Scaffolding Reads

    #
    #  Create symlinks to the categories.
    #
    #  Name changes should change:
    #  1)  the label used in this script ("low-quality-draft")
    #  2)  the directory name in the symlink
    #  3)  the actual directories (both "_genomeark-low-quality-draft" and "genomeark-low-quality-draft")
    #  4)  the index.html in genomeark-low-quality-draft/  (both the title and the collection name)
    #  5)  _config.yml (in two places) and index.html in the root directory.
    #

    if ($updatePages) {
        foreach my $proj ("genomeark", getAllProjectNames()) {
            unlink("$ga/_$proj-all/$name.md");
            unlink("$ga/_$proj-curated-assembly/$name.md");
            unlink("$ga/_$proj-hqdraft-assembly/$name.md");
            unlink("$ga/_$proj-lqdraft-assembly/$name.md");
            unlink("$ga/_$proj-raw-data-only/$name.md");
        }

        sub makeLink ($$$) {       #  Create a link to the actual data file, which is ALWAYS in "../_genomeark/".
            my $name = shift @_;   #  Species name.
            my $proj = shift @_;   #  Project: 'genomeark', 'vg', 't2t', etc.
            my $catg = shift @_;   #  Category: 'all', 'curated-assembly', etc.
            symlink("../_genomeark/$name.md", "$ga/_$proj-$catg/$name.md") or die "Failed to make symlink for '$ga/_$proj-$catg/$name.md': $!.\n";
        }

        foreach my $proj ("genomeark", split '\s+', getProjectID($data{"name_"})) {
            makeLink($name, $proj, "all");

            if    ($data{"assembly_status"} eq "curated") { makeLink($name, $proj, "curated-assembly") }
            elsif ($data{"assembly_status"} eq "hqdraft") { makeLink($name, $proj, "hqdraft-assembly") }
            elsif ($data{"assembly_status"} eq "lqdraft") { makeLink($name, $proj, "lqdraft-assembly") }
            else                                          { makeLink($name, $proj, "raw-data-only")    }
        }
    }

    #  And reset the classifications to strings we can use in the display.

    $data{"data_status"}     = "<em style=\"color:red\">no data</em>"                          if ($data{"data_status"}) eq "none";
    $data{"data_status"}     = "<em style=\"color:orange\">some data</em>"                     if ($data{"data_status"}) eq "some";
    $data{"data_status"}     = "<em style=\"color:green\">all data</em>"                       if ($data{"data_status"}) eq "all";

    $data{"data_status"}  = "'";
    $data{"data_status"} .= "<em style=\"color:" . (($hasHQ) ? "green" : "red") . "\">HQ Long</em> ::: ";
    $data{"data_status"} .= "<em style=\"color:" . (($hasLR) ? "green" : "red") . "\">Long</em> ::: ";
    $data{"data_status"} .= "<em style=\"color:" . (($hasSR) ? "green" : "red") . "\">Short</em> ::: ";
    $data{"data_status"} .= "<em style=\"color:" . (($hasPH) ? "green" : "red") . "\">Phasing</em> ::: ";
    $data{"data_status"} .= "<em style=\"color:" . (($hasSC) ? "green" : "red") . "\">Scaffolding</em>";
    $data{"data_status"} .= "'";

    $data{"assembly_status"} = "<em style=\"color:red\">none</em>"                    if ($data{"assembly_status"} eq "none");
    $data{"assembly_status"} = "<em style=\"color:red\">low-quality draft</em>"       if ($data{"assembly_status"} eq "lqdraft");
    $data{"assembly_status"} = "<em style=\"color:orange\">high-quality draft</em>"   if ($data{"assembly_status"} eq "hqdraft");
    $data{"assembly_status"} = "<em style=\"color:green\">curated</em>"               if ($data{"assembly_status"} eq "curated");

    #  If no date set -- no raw data, no assemblies, no anything -- set it to the
    #  date of this update (genomeark.ls's date).

    if ($data{"last_updated"} == 0) {
        print "WARNING: no files; last_updated set to 'now'.\n";
        #$data{"last_updated"} = 0;
    }

    #  Done.  Write the output to _genomeark.  Note that the _proj-whatever/
    #  directories are just symlinks; this is the only copy of the data
    #  itself.

    if ($updatePages) {
        saveData(\%data, "$ga/_genomeark/$name.md");
    }

    print "\n";
}




#  And whatever errors we saved.

if (scalar(@potentialErrors > 0)) {
    print "\n";
    print "----------------------------------------\n";
    print "Potential errors found:\n";
    foreach my $l (@potentialErrors) {
        print "  $l";
    }
    print "\n";
}


if (scalar(@unknownFiles > 0)) {
    print "\n";
    print "----------------------------------------\n";
    print "Unknown files:\n";
    foreach my $l (sort @unknownFiles) {
        print "  $l";
    }
    print "\n";
}


if (scalar(keys %missingData > 0)) {
    print "\n";
    print "----------------------------------------\n";
    print "Species with missing data files:\n";
    foreach my $l (sort keys %missingData) {
        printf "  %3d - %s\n", $missingData{$l}, $l;
    }
    print "\n";
}
