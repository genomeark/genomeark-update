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
my $downloadData       = 1;
my $downloadAssemblies = 1;
my $downloadSize       = 2 * 1024 * 1024 * 1024;

while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if    ($arg eq "--no-download") {       #  Don't do any downloads.
        $downloadData = 0;
        $downloadAssemblies = 0;
    }
    elsif ($arg eq "--max-download") {       #  Limit data downloads to X MB.
        $downloadSize  = shift @ARGV;
        $downloadSize *= 1024 * 1024;
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

@speciesList = discoverSpecies(@speciesList);   #  List of all species in genomeark-metadata or supplied on command line.

my $lastupdate  = loadBucketFiles();        #  Private list of all files in the bucket.
my $nGenBank    = loadGenbankMap();         #  Private map from ToLID to GenBank accession.

my %missingData;                            #  List of species that need to download data.
my %templateMeta;                           #  List of species that have only template metadata.
my @potentialErrors;                        #  List of fatal errors we encountered.
my @unknownFiles;                           #  List of files we don't know how to process.

foreach my $species (@speciesList) {
    print "\n";
    print "----------------------------------------\n";
    print "Processing species  $species\n";

    #  Load metadata, copy some of it into the .md data output.

    my %meta = ();
    my %data = ();
    my $name = loadSpeciesMetadata(\%data, $species, \%meta, \%templateMeta, \@potentialErrors);

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

        next   if (! isAssemblyFile($filename, 0));

        if ($data{"last_updated"} < $filesecs) {
            $data{"last_updated"} = $filesecs;
        }

        rememberLatestAssembly($filesecs, $filesize, $filename, \%data, \@potentialErrors);
    }


    print "\n";
    print "----------\n";
    print "Assemblies, pass 2: summarize assembly\n";

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        next   if (! isAssemblyFile($filename, 0));

        summarizeAssembly($filesecs, $filesize, $filename, \%data, \@potentialErrors, \%missingData, $downloadAssemblies);
    }


    print "\n";
    print "----------\n";
    print "Genomic Data, accumulation\n";

    my %tiFiles = ();   #  Map from "type:indiv" to list of "type:indiv filesize filename\0"
    my %tiBytes = ();   #  Map from "type:indiv" to sum of data size
    my %tiIndiv = ();   #  Map from "type:indiv" to list of "species-name/individual-name\0"

    #  Save the filenames of any hifi fastq data.
    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        rememberHiFiName($speciesFiles[$ii], \%data);
    }

    #  Place filesize and filename into lists keyed by 'datatype:individual'.
    #  Only seqFiles is used.
    #    "type:indiv filesize filename\0"
    #
    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        accumulateData($speciesEpoch[$ii], $speciesSizes[$ii], $speciesFiles[$ii],
                       \%data,
                       \%tiFiles, \%tiBytes, \%tiIndiv,
                       \@potentialErrors);
    }


    print "\n";                                  #  Estimate the number of bases in all raw data files by
    print "----------\n";                        #  examining a handful and scaling the rest.
    print "Genomic Data, bases estimation\n";    #  BIONANO is different and not reported.
    print "   (S: summary exists)\n";
    print "   (E: used for byte:base estimation)\n";
    print "\n";

    #  Load existing summaries from 'species/$name/genomic_data.summary', then scan
    #  the list of files for any additional summary files.
    #
    loadSummaries(\%data, \%tiFiles, \@potentialErrors);

    #  Load any additional summaries, or update existing with new/better summaries.
    #
    recoverSummaries(\%data, \%tiFiles, \@potentialErrors);

    #  Compute scale factors to convert 'bytes' to 'bases' for each
    #  type-x-individual.  Populates:
    #    $data{'data_${ti}_scale'}     -- scale factor for this ti
    #    $data{'data_${ti}_bases'}     -- estimated bases for all files with the same ti
    #
    foreach my $ti (sort keys %tiFiles) {
        estimateRawDataScaling(\%data, $ti, $tiFiles{$ti},
                               \@potentialErrors, \%missingData, $downloadData, $downloadSize);
    }

    #  Update data on disk.
    #
    writeSummaries(\%data);

    print "\n";
}

#
#  Output whatever errors we saved.
#

if (scalar(@potentialErrors > 0)) {
    print "\n";
    print "----------------------------------------\n";
    print "Potential errors found:\n";
    foreach my $l (sort @potentialErrors) {
        print "  $l";
    }
    print "\n";
}


if (scalar(keys %templateMeta > 0)) {
    print "\n";
    print "----------------------------------------\n";
    print "Species with template metadata:\n";
    foreach my $l (sort keys %templateMeta) {
        printf "  %s\n", $l;
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
        my @m = split '\s+', $missingData{$l};
        my %c;

        foreach my $m (@m) {
            $c{$m}++;
        }

        foreach my $k (keys %c) {
            printf "  %40s  %16s  %d files\n", $l, $k, $c{$k};
        }
    }
    print "\n";
}

exit(0);
