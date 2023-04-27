package GenomeArkUpdate;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(discoverSpecies loadBucketFiles getBucketFiles isGenomicDataFile isAssemblyFile);

use strict;
use warnings;

use GenomeArkUtility;

use Time::Local;
use List::Util qw(min max);

#  In memory cache of files in genomeark.
my @genomeArkEpoch;
my @genomeArkSizes;
my @genomeArkFiles;

#  In memory cache of files in genomeark, filtered to the current species.
#my @speciesEpoch;
#my @speciesSizes;
#my @speciesFiles;



#
#  Discover species names from genomeark-metadata/species-list.
#
#  If input names are supplied, they are checked against species-list and
#  errors are reported if any are not found.
#  

sub discoverSpecies (@) {
    my %species;   #  
    my @species;   #  List of valid species, return value.
 
    print "\n";
    print "Discovering species names.\n";

    foreach my $s (@_) {
        $species{$s}++   if ($s ne "");
    }

    open(SL, "< downloads/species-list");
    if (scalar(keys %species) == 0) {        #  If no input names, just slurp in the
        while (<SL>) {                       #  whole list.
            s/^\s+//;                        #
            s/\s+$//;
            push @species, $_;
        }
    }
    else {
        while (<SL>) {                       #  Otherwise, check the name against
            s/^\s+//;                        #  the hash and add it if found.
            s/\s+$//;                        #
            if (exists($species{$_})) {
                push @species, $_;
                delete $species{$_};
            }
        }
    }
    close(SL);

    print "  Found ", scalar(@species), " species.\n";

    if (scalar(keys %species) > 0) {         #  Any remaining names are errors.
        print "\n";
        print "  ERROR: Did NOT find ", scalar(keys %species), " species:\n";
        foreach my $s (keys %species) {
            print "  ERROR:   $s\n";
        }
        exit(1);
    }

    return(@species);
}



#
#  Loads the (filtered) list of files in the bucket into global arrays:
#    genomeArkFiles
#    genomeArkSizes
#    genomeArkEpoch
#
#  Returns the time stamp on 'downloads/genomeark.ls' which we interpret as
#  the time of the update.
#

sub loadBucketFiles () {
    my $lastupdate = (stat("downloads/genomeark.ls"))[9];

    open(GLS, "< downloads/genomeark.ls");
    while (<GLS>) {
        my ($filedate, $filetime, $filesize, $filename, $filesecs) = split '\s+', $_;

        my @fileComps   = split '/', $filename;
        my $speciesName = $fileComps[1];
        my $asmName     = $fileComps[3];
        my $seconds     = 0;

        if ("$filedate $filetime" =~ m/(\d\d\d\d)-(\d\d)-(\d\d)\s+(\d\d):(\d\d):(\d\d)/) {
            my ($yr, $mo, $da, $hr, $mn, $sc) = ($1, $2, $3, $4, $5, $6);

            $filesecs = timelocal($sc, $mn, $hr, $da, $mo-1, $yr); 
        } else {
            die "Failed to parse date ('$filedate') and time ('$filetime') for file '$filename'\n";
        }

        push @genomeArkFiles, $filename;
        push @genomeArkSizes, $filesize;
        push @genomeArkEpoch, $filesecs;
    }
    close(GLS);

    return($lastupdate);
}



#
#  Extracts files for a single speices from genomeArkFiles et alia.
#

sub getBucketFiles ($$$$) {
    my $sname  = shift @_;
    my $sfiles = shift @_;
    my $ssizes = shift @_;
    my $sepoch = shift @_;
    my $lastupdate  = 0;

    undef @$sfiles;
    undef @$ssizes;
    undef @$sepoch;

    for (my $ii=0; $ii<scalar(@genomeArkFiles); $ii++) {
        next if ($genomeArkFiles[$ii] !~ m!$sname!);

        push @$sfiles, $genomeArkFiles[$ii];
        push @$ssizes, $genomeArkSizes[$ii];
        push @$sepoch, $genomeArkEpoch[$ii];

        $lastupdate = max($lastupdate, $genomeArkEpoch[$ii]);
    }

    return($lastupdate);
}



#
#  Return true/false if the filename looks like it's a data file
#  or assembly we should be reporting.
#

sub isGenomicDataFile ($) {
    my $filename   = shift @_;

    return($filename =~ m!genomic_data!)
}


#
#  Returns true if type == "sequence" and file is either mito or genomic.
#  Returns true if type == "metadata" and file is sample or assembly metadata.
#

sub isAssemblyFile ($$$) {
    my $filename = shift @_;
    my $type     = shift @_;
    my $errors   = shift @_;

    return(0)   if ($filename =~ m!genomic_data!);
    return(0)   if ($filename =~ m!additional_haplotigs!);

    return(0)   if ($filename =~ m!assembly_v0!);
    return(0)   if ($filename =~ m!bwgenome_!);
    return(0)   if ($filename =~ m!concatenated-!);


    my $isMeta = (($filename =~ m!species/.*/${ToLIDregex}/.*assembly.*/${ToLIDregex}.*\.y[a]{0,1}ml$!i) ||
                  ($filename =~ m!species/.*/${ToLIDregex}/${ToLIDregex}.*\.y[a]{0,1}ml$!i));

    my $isMito = (($filename =~ m!species/.*/${ToLIDregex}/.*assembly.*/${ToLIDregex}.MT.\d\d\d\d\d\d\d\d\.fasta\.gz$!i) ||
                  ($filename =~ m!species/.*/${ToLIDregex}/.*assembly.*/${ToLIDregex}.\w+.\w+\.\d\d\d\d\d\d\d\d\.MT\.fasta\.gz$!i));

    my $isGeno =  ($filename =~ m!species/.*/${ToLIDregex}/.*assembly.*/${ToLIDregex}.\w+.\w+.\d\d\d\d\d\d\d\d\.fasta\.gz$!i);


    if    ($type eq "metadata") {
        return(1)   if ($isMeta);
        return(0)   if ($isMito || $isGeno);
    }
    elsif ($type eq "sequence") {
        return(0)   if ($isMeta);
        return(1)   if ($isMito || $isGeno);
    }
    else {
        die "Expecting 'metadata' or 'sequence' in isAssemblyFile().\n";
    }


    push @$errors, "  Unknown assembly file '$filename' for '$type'; meta='$isMeta' mito='$isMito' geno='$isGeno'\n";
    return(0);
}


return(1);
