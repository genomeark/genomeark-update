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

foreach my $arg (@ARGV) {
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
my @potentialErrors;                        #  List of fatal errors we encountered.
my @unknownFiles;                           #  List of files we don't know how to process.

foreach my $species (@speciesList) {
    my @speciesFiles;
    my @speciesSizes;
    my @speciesEpoch;

    #  Rebuild each species.md file.

    my %seqFiles;
    my %seqBytes;
    my %seqIndiv;
    my %meta;
    my %data;


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

        scanDataName($filesecs, $filesize, $filename, \%data);
    }

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        next   if (! isGenomicDataFile($filename));

        accumulateData($filesecs, $filesize, $filename, \%data, \%seqFiles, \%seqBytes, \%seqIndiv, \@potentialErrors);
    }


    print "\n";                                  #  Estimate the number of bases in all raw data files by
    print "----------\n";                        #  examining a handful and scaling the rest.
    print "Genomic Data, bases estimation\n";    #  BIONANO is different and not reported.

    loadSummaries(\%data, \%seqFiles, \@potentialErrors);   #  Load existing summaries.

    foreach my $type (@dataTypes) {
        estimateRawDataScaling(\%data, $type, $seqFiles{$type}, \@potentialErrors, \%missingData, $downloadData, $downloadSize);
    }

    writeSummaries(\%data);

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
        my @m = split '\s+', $missingData{$l};
        my %c;

        foreach my $m (@m) {
            $c{$m}++;
        }

        foreach my $k (keys %c) {
            printf "  %40s  %12s  %d files\n", $l, $k, $c{$k};
        }
    }
    print "\n";
}
