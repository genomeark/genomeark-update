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

#use YAML::XS;

#use Slack;
#use Update;
#use Schedule;

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

    die "setCoverage input bases not defined.\n"  if (!defined($b));
    die "setCoverage genome size not defined.\n"  if (!defined($g));

    if ($g > 0) {
        return(sprintf("%.2fx", $b / $g));
    } else {
        return("N/A");
    }
}


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

my $ga;

$ga = "../genomeark.github.io"   if (-e "../genomeark.github.io");   #  Legacy location.
$ga =  "./genomeark.github.io"   if (-e  "./genomeark.github.io");   #  Submodule location.

die "Didn't find genomeark.github.io.\n"  if (!defined($ga));

@speciesList = discoverSpecies(@speciesList);   #  List of all species in genomeark-metadata or supplied on command line.

my $lastupdate  = loadBucketFiles();        #  Private list of all files in the bucket.
my $nGenBank    = loadGenbankMap();         #  Private map from ToLID to GenBank accession.
my $nProject    = loadProjectMetadata();    #  Private list of known genome projects.

my %missingData;                            #  List of species that need to download data.
my %templateMeta;                           #  List of species that have only template metadata
my %noProjectID;                            #  List of species that have no project ID
my @potentialErrors;                        #  List of fatal errors we encountered.
my @potentialDataErrors;                    #  List of fatal errors we encountered.
my @unknownFiles;                           #  List of files we don't know how to process.

#  Create output directories, now that we know what projects exist.

system("mkdir -p $ga/_genomeark");

foreach my $proj ("genomeark", getAllProjectNames()) {
    system("mkdir -p $ga/_$proj-all");
    system("mkdir -p $ga/_$proj-curated-assembly");
    system("mkdir -p $ga/_$proj-draft-assembly");
    system("mkdir -p $ga/_$proj-raw-data-only");
}

#  Make index pages for each of the categories.

foreach my $proj ("genomeark", getAllProjectNames()) {
    makeIndexPage($proj, "$ga", "$proj-all",              "All Species");
    makeIndexPage($proj, "$ga", "$proj-curated-assembly", "Species with Curated Assemblies");
    makeIndexPage($proj, "$ga", "$proj-draft-assembly",   "Species with Draft Assembles");
    makeIndexPage($proj, "$ga", "$proj-raw-data-only",    "Species without an Assembly");
}

#  Iterate over all species and do work!

foreach my $species (@speciesList) {

    print "\n";
    print "----------------------------------------\n";
    print "Processing species  $species\n";

    #  Load metadata, copy some of it into the .md data output.

    my %data = ();
    my $name = loadSpeciesMetadata(\%data, $species, \%templateMeta, \@potentialErrors);

    setGenbankIDs(\%data);

    #  Extract the files associated with this species, returns the latest file-date.

    my @speciesFiles = ();  #  A one-to-one list of file-names, file-sizes and 
    my @speciesSizes = ();  #  file-dates for all files associated with this
    my @speciesEpoch = ();  #  species.

    $data{"last_updated"} = getBucketFiles($name, \@speciesFiles, \@speciesSizes, \@speciesEpoch);

    print "  with ", scalar(@speciesFiles), " files.\n";

    ########################################

    print "\n";
    print "----------\n";
    print "Assemblies, pass 1: remember the latest assembly\n";

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        if (isAssemblyFile($filename, "sequence", \@potentialErrors)) {
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

        next   if (! isAssemblyFile($filename, "sequence", \@potentialErrors));

        importAssemblySummary($filesecs, $filesize, $filename, \%data, \@potentialErrors, \%missingData);
    }

    print "\n";
    print "----------\n";
    print "Assemblies, pass 2: import metadata\n";

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];
        my $dir;

        next   if (! isAssemblyFile($filename, "metadata", \@potentialErrors));

        print STDERR "METADATA '$filename'\n";

        {
            my $md = loadYAMLasString("downloads/$filename");
            my @fc = split '/', $filename;

            my $mdh = {};                              #  Make an entry for this metadata file.
            $$mdh{"ident"} = "md$ii";                  #
            $$mdh{"title"} = join '/', @fc[3..$#fc];   #  Remove stuff we already know.
            $$mdh{"data"}  = $md;                      #

            if (! exists($data{"mds"})) {              #  Add an empty array to the data.
                $data{"mds"} = [];
            }

            my $arr = $data{"mds"};                    #  Push entry onto list.
            push @$arr, $mdh;
        }
    }

    #  Update genome size to one of the assembly sizes, if needed.

    $data{"genome_size"} = $data{"pri1length"}   if ($data{"genome_size"} == 0) && defined($data{"pri1length"});
    $data{"genome_size"} = $data{"pri2length"}   if ($data{"genome_size"} == 0) && defined($data{"pri2length"});
    $data{"genome_size"} = $data{"pri3length"}   if ($data{"genome_size"} == 0) && defined($data{"pri3length"});
    $data{"genome_size"} = $data{"pri4length"}   if ($data{"genome_size"} == 0) && defined($data{"pri4length"});
    $data{"genome_size"} = $data{"pri5length"}   if ($data{"genome_size"} == 0) && defined($data{"pri5length"});
    $data{"genome_size"} = $data{"pri6length"}   if ($data{"genome_size"} == 0) && defined($data{"pri6length"});

    $data{"genome_size_display"} = prettifyBases($data{"genome_size"});

    ########################################

    print "\n";
    print "----------\n";
    print "Genomic Data, accumulation\n";

    my %tiFiles = ();   #  Map from "type:indiv" to list of "type:indiv filesize filename\0"
    my %tiBytes = ();   #  Map from "type:indiv" to sum of data size
    my %tiIndiv = ();   #  Map from "type:indiv" to list of "species-name/individual-name\0"

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        accumulateData($speciesEpoch[$ii], $speciesSizes[$ii], $speciesFiles[$ii],
                       \%data,
                       \%tiFiles, \%tiBytes, \%tiIndiv,
                       \@potentialDataErrors);
    }


    print "\n";
    print "----------\n";
    print "Page Creation\n";

    #  Load summaries of data files into %data.
    #    $data{'data_${ti}_scale'}     -- scale factor for this ti
    #    $data{'data_${ti}_bases'}     -- estimated bases for all files with the same ti
    #
    #  In update-pages.pl, the rest of the summary data is private.  DO NOT CALL writeSummaries()!
    #
    loadSummaries(\%data, \%tiFiles, \@potentialErrors);

    #  Clear data for the bytes/bases in each ti and type.
    #
    foreach my $ti (keys %tiFiles) {
        my ($ty, $in) = split ':', $ti;

        $data{"data_${ti}_bytes_v"} = 0;
        $data{"data_${ti}_bases_v"} = 0;
        $data{"data_${ty}_bytes_v"} = 0;
        $data{"data_${ty}_bases_v"} = 0;
    }

    #  Sum bytes and estimated bases in each ti and type.
    #
    foreach my $ti (keys %tiFiles) {
        my @tiflist = split '\0', $tiFiles{$ti};

        foreach my $tif (@tiflist) {
            my ($ti, $bytes, $filename) = split '\s+', $tif;
            my ($ty, $id) = split ':', $ti;

            $data{"data_${ti}_bytes_v"} += $bytes;
            $data{"data_${ty}_bytes_v"} += $bytes;

            $data{"data_${ti}_bases_v"} += $bytes * $data{"data_${ti}_scale"};
            $data{"data_${ty}_bases_v"} += $bytes * $data{"data_${ti}_scale"};

            #print "$ti ", $data{"data_${ti}_bytes_v"}, " ", $data{"data_${ti}_bases_v"}, " $ty ", $data{"data_${ty}_bytes_v"}, " ", $data{"data_${ti}_bases_v"}, "\n";
        }
    }

    #  Convert the summed bytes/bases to display values.
    #    $ti is unique, but $ty is not.
    #    
    foreach my $ti (keys %tiFiles) {
        my ($ty, $id) = split ':', $ti;

        die  if (! exists($data{"data_${ti}_bytes_v"}));

        $data{"data_${ti}_bytes"}    = prettifySize($data{"data_${ti}_bytes_v"});
        $data{"data_${ti}_bases"}    = prettifyBases($data{"data_${ti}_bases_v"});
        $data{"data_${ti}_bases"}    = "N/A" if ($ty eq "bionano");
        $data{"data_${ti}_coverage"} = setCoverage($data{"data_${ti}_bases_v"}, $data{"genome_size"});

        delete $data{"data_${ti}_bytes_v"};
        delete $data{"data_${ti}_bases_v"};

        next if (! exists($data{"data_${ty}_bytes_v"}));

        $data{"data_${ty}_bytes"}    = prettifySize($data{"data_${ty}_bytes_v"});
        $data{"data_${ty}_bases"}    = prettifyBases($data{"data_${ty}_bases_v"});
        $data{"data_${ty}_bases"}    = "N/A" if ($ty eq "bionano");
        $data{"data_${ty}_coverage"} = setCoverage($data{"data_${ty}_bases_v"}, $data{"genome_size"});

        delete $data{"data_${ty}_bytes_v"};
        delete $data{"data_${ty}_bases_v"};
    }


    #  Build a list of links to data for each ti.
    #
    #  https://genomeark.s3.amazonaws.com/index.html?prefix=species/Gallus_gallus
    #  https://42basepairs.com/browse/s3/genomeark/species/Gallus_gallus
    #
    foreach my $ti (keys %tiFiles) {
        my ($ty, $in) = split ':', $ti;

        #  Argh, the location of the reads is not the same as 'type'.
        my $path = $ty;
        $path = "pacbio_hifi"   if ($ty eq "pacbiohifi_dcfqgz");
        $path = "pacbio_hifi"   if ($ty eq "pacbiohifi_fqgz");
        $path = "pacbio_hifi"   if ($ty eq "pacbiohifi_bam");
        $path = "pacbio_hifi"   if ($ty eq "pacbiohifi_clr");
        $path = "ont_duplex"    if ($ty eq "ontduplex");

        my $sp = $data{"name_"};

        $data{"data_${ti}_links"} .= "s3://genomeark/species/$sp/$in/genomic_data/$path/<br>";

        $data{"data_${ti}_s3url"} .= "https://genomeark.s3.amazonaws.com/index.html?prefix=species/$sp/$in/genomic_data/$path/";
        $data{"data_${ti}_s3gui"} .= "https://42basepairs.com/browse/s3/genomeark/species/$sp/$in/genomic_data/$path/";
    }

    #
    #  Build a list of the data types that exist for display in the species list pages.
    #

    sub hasBases ($@) {
        my $data = shift @_;
        foreach my $type (@_) {
            next   if (!defined($$data{"data_${type}_bases"}));
            next   if ($$data{"data_${type}_bases"} eq "N/A");
            return(1);
        }
        return(0);
    }

    my @have;
    push @have, "<em style=\"color:forestgreen\">PacBio CLR</em>"  if hasBases(\%data, "pacbio");
    push @have, "<em style=\"color:forestgreen\">PacBio HiFi</em>" if hasBases(\%data, "pacbiohifi_fqgz", "pacbiohifi_bam");
    push @have, "<em style=\"color:forestgreen\">ONT Simplex</em>" if hasBases(\%data, "ont");
    push @have, "<em style=\"color:forestgreen\">ONT Duplex</em>"  if hasBases(\%data, "ontduplex");
    push @have, "<em style=\"color:forestgreen\">10x</em>"         if hasBases(\%data, "10x");
    push @have, "<em style=\"color:forestgreen\">Bionano</em>"     if hasBases(\%data, "bionano");
    push @have, "<em style=\"color:forestgreen\">Arima</em>"       if hasBases(\%data, "arima");
    push @have, "<em style=\"color:forestgreen\">Dovetail</em>"    if hasBases(\%data, "dovetail");
    push @have, "<em style=\"color:forestgreen\">Phase</em>"       if hasBases(\%data, "phase");
    push @have, "<em style=\"color:forestgreen\">Illumina</em>"    if hasBases(\%data, "illumina");

    if (scalar(@have) > 0) {
        $data{"data_status"}  = "'" . join(" ::: ", @have) . "'";
    } else {
        $data{"data_status"}  = "'<em style=\"color:maroon\">No data</em>'";
    }

    #
    #  Remove individual name from keys.
    #

    my $sn = $data{"short_name"};

    foreach my $key (keys %data) {
        if ($key =~ m/^(data_.*):$sn(\d_.*)$/) {
            $data{"$1-$2"} = $data{$key};
            delete $data{$key};
        }
    }


    #
    #  Create symlinks to the categories.
    #
    #  Name changes should change:
    #  1)  the label used in this script ("draft-assembly")
    #  2)  the directory name in the symlink
    #  3)  the actual directories (both "_genomeark-draft-assembly" and "genomeark-draft-assembly")
    #  4)  the index.html in genomeark-draft-assembly/  (both the title and the collection name)
    #  5)  _config.yml (in two places) and index.html in the root directory.
    #

    if ($updatePages) {
        foreach my $proj ("genomeark", getAllProjectNames()) {
            unlink("$ga/_$proj-all/$name.md");
            unlink("$ga/_$proj-curated-assembly/$name.md");
            unlink("$ga/_$proj-draft-assembly/$name.md");
            unlink("$ga/_$proj-raw-data-only/$name.md");
        }

        sub makeLink ($$$) {       #  Create a link to the actual data file, which is ALWAYS in "../_genomeark/".
            my $name = shift @_;   #  Species name.
            my $proj = shift @_;   #  Project: 'genomeark', 'vg', 't2t', etc.
            my $catg = shift @_;   #  Category: 'all', 'curated-assembly', etc.
            symlink("../_genomeark/$name.md", "$ga/_$proj-$catg/$name.md") or die "Failed to make symlink for '$ga/_$proj-$catg/$name.md': $!.\n";
        }

        my $pl = getProjectID($data{"name_"});
        my @pl = ( "genomeark" );

        if ($pl) {
            push @pl, split '\s+', $pl;
        } else {
            $noProjectID{$data{'name_'}}++;
            #push @potentialErrors, "  Species $data{'name_'} has no projectID.\n";
        }

        foreach my $proj (@pl) {
            makeLink($name, $proj, "all");

            if    ($data{"assembly_status"} eq "curated") { makeLink($name, $proj, "curated-assembly") }
            elsif ($data{"assembly_status"} eq "draft")   { makeLink($name, $proj, "draft-assembly")   }
            else                                          { makeLink($name, $proj, "raw-data-only")    }
        }
    }

    #
    #  Convert the assembly status to colored html (we use the raw label just above here).
    #

    $data{"assembly_status"} = "<em style=\"color:maroon\">No assembly</em>"    if ($data{"assembly_status"} eq "none");
    $data{"assembly_status"} = "<em style=\"color:orangered\">Draft</em>"       if ($data{"assembly_status"} eq "draft");
    $data{"assembly_status"} = "<em style=\"color:forestgreen\">Curated</em>"   if ($data{"assembly_status"} eq "curated");

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

#
#  Output whatever errors we saved.
#

if (scalar(keys %noProjectID > 0)) {
    print "\n";
    print "----------------------------------------\n";
    print "Found ", scalar(keys %noProjectID), " species with no project ID assigned.\n";
}

if (scalar(@potentialErrors > 0)) {
    print "\n";
    print "----------------------------------------\n";
    print "Potential errors found (not showing ", scalar(@potentialDataErrors), " data errors; see scan-bucket output):\n";
    foreach my $l (sort @potentialErrors) {
        print "  $l";
    }
    print "\n";
}

exit(0);
