#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use lib "$FindBin::RealBin";

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

use Time::Local;
#use List::Util;

#
#  Parse options.
#

my @speciesList;
my $downloadData       = 1;
my $downloadAssemblies = 1;

foreach my $arg (@ARGV) {
    if    ($arg eq "--no-download") {
        $downloadData = 0;
        $downloadAssemblies = 0;
    }
    elsif ($arg eq "--no-download") {
    }
    elsif ($arg eq "--no-download") {
    }
    else {
        push @speciesList, $arg;
    }
}


#
#  Main
#

die "ERROR: 'downloads/genomeark.ls' doesn't exist, can't update.\n"            if (! -e "downloads/genomeark.ls");
die "ERROR: 'genomeark-metadata/species-list' doesn't exist, can't update.\n"   if (! -e "genomeark-metadata/species-list");

@speciesList = discoverSpecies(@speciesList);   #  List of all species in genomeark-metadata or supplied on command line.

my $lastupdate  = loadBucketFiles();        #  Private list of all files in the bucket.
my $nGenBank    = loadGenbankMap();         #  Private map from ToLID to GenBank accession.

my @potentialErrors;                        #  List of fatal errors we encountered.
my @unknownFiles;                           #  List of files we don't know how to process.

foreach my $species (@speciesList) {
    my @speciesFiles;
    my @speciesSizes;
    my @speciesEpoch;

    #  Rebuild each species.md file.

    my %seqFiles;   #  Number of files for each type of data
    my %seqBytes;   #  Total size in bytes for each type
    my %seqIndiv;   #  List of individuals with data
    my %meta;
    my %data;


    print "\n";
    print "----------------------------------------\n";
    print "Processing species  $species\n";

    loadSpeciesMetadata(\%data, $species, \%meta);    #  Load genomeark-metadata.
    setGenbankIDs(\%data);                            #  Find GenBank accessions.

    $data{"last_updated"} = getBucketFiles($data{"name_"}, \@speciesFiles, \@speciesSizes, \@speciesEpoch);

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
    print "Assemblies, pass 2: summarize assembly\n";

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        if (isAssemblyFile($filename, 0)) {
            summarizeAssembly($filesecs, $filesize, $filename, \%data, \@potentialErrors, $downloadAssemblies);
        }
    }


    print "\n";
    print "----------\n";
    print "Genomic Data, accumulation\n";

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        if (isGenomicDataFile($filename)) {
            if ($data{"last_updated"} < $filesecs) {
                $data{"last_updated"} = $filesecs;
            }
            if ((! exists($data{"last_raw_data"})) ||
                ($data{"last_raw_data"} < $filesecs)) {   #  If this isn't set, Raw Data shows
                $data{"last_raw_data"} = $filesecs;       #  "No data.".
            }

            accumulateData($filesize, $filename, \%seqFiles, \%seqBytes, \%seqIndiv, \@potentialErrors);
        }
    }


    print "\n";                                  #  Estimate the number of bases in all raw data files by
    print "----------\n";                        #  examining a handful and scaling the rest.
    print "Genomic Data, bases estimation\n";    #  BIONANO is different and not reported.

    loadSummaries($data{"name_"}, \%seqFiles);   #  Load any summaries we have, so we can skip it in estimateRawDataScaling().

    estimateRawDataScaling(\%data, "10x",      $seqFiles{"10x"},      \@potentialErrors, $downloadData);  #  seqFiles is a \0 separated list of
    estimateRawDataScaling(\%data, "arima",    $seqFiles{"arima"},    \@potentialErrors, $downloadData);  #    'filesize \s datafile
    estimateRawDataScaling(\%data, "dovetail", $seqFiles{"dovetail"}, \@potentialErrors, $downloadData);  #
    estimateRawDataScaling(\%data, "illumina", $seqFiles{"illumina"}, \@potentialErrors, $downloadData);  #  datafile is 'species/NAME/INDIVIDUAL/genomic_data/TYPE/FILE
    estimateRawDataScaling(\%data, "nanopore", $seqFiles{"nanopore"}, \@potentialErrors, $downloadData);
    estimateRawDataScaling(\%data, "pbclr",    $seqFiles{"pbclr"},    \@potentialErrors, $downloadData);
    estimateRawDataScaling(\%data, "pbhifi",   $seqFiles{"pbhifi"},   \@potentialErrors, $downloadData);
    estimateRawDataScaling(\%data, "phase",    $seqFiles{"phase"},    \@potentialErrors, $downloadData);

    writeSummaries($data{"name_"});

    #print "\n";                  #  Report any unprocessed files.
    #print "----------\n";
    #print "Unknown files:\n";

    for (my $ii=0; $ii<scalar(@speciesFiles); $ii++) {
        my $filesecs = $speciesEpoch[$ii];
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        next  if (isGenomicDataFile($filename));
        next  if (isAssemblyFile($filename, 0));

        #print "  UNKNOWN $filename\n";

        push @unknownFiles, "  $filename\n";
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
