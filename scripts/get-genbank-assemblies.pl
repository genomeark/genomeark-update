#!/usr/bin/env perl

#  Use NCBI datasets command to create the file that Bri used to create
#  with esearch.

use strict;
use JSON;

my $Usage = 'get-genbank-assemblies.pl <BioProject accession>';

$#ARGV==0
    or die $Usage;

my $bioproject = $ARGV[0];

my $jsonassemblies = `datasets summary genome accession $bioproject`;

#print $jsonassemblies;

my $decodedassemblies = decode_json($jsonassemblies);

if (ref $decodedassemblies eq 'HASH') {
    for my $hashkey (keys %{$decodedassemblies}) {
	    my $keyref = ref $decodedassemblies->{$hashkey} || $decodedassemblies->{$hashkey};
	    print "$hashkey\t$keyref\n";
	    if ($keyref eq 'ARRAY') {
                my @values = @{$decodedassemblies->{$hashkey}};
		foreach my $value (@values) {
		    my $valueref = ref $value;
		    if ($valueref eq 'HASH') {
			    my $accession = $value->{'accession'};
			    my $assemblyinfo = $value->{'assembly_info'};
			    my $assemblyinforef = ref $assemblyinfo;

			    my $assemblytype = $assemblyinfo->{'assembly_type'};
			    my $assemblyname = $assemblyinfo->{'assembly_name'};
			    print "$accession\t$assemblyname\t$assemblytype\n";
	            }
		}
            }
    }
}

