package GenomeArkGenbank;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(loadGenbankMap setGenbankIDs);

use strict;
use warnings;

use Time::Local;
use List::Util;

my %genbankMap;

#
#  Load a map from TolID to GenBankID.
#

sub loadGenbankMap () {

    if (! -e "species-data/genbank.map") {
        die "No species-data/genbank.map found.\n";
    }

    open(GB, "< species-data/genbank.map");
    while (<GB>) {
        chomp;

        if ($_ =~ m/^(GCA_\d+.\d+)\s+(.[A-Z].....)(\d)\s+(...)$/) {
            my $acc = $1;
            my $nam = $2;
            my $ind = $3;
            my $typ = $4;

            if (exists($genbankMap{"$nam$typ"})) {
                $genbankMap{"$nam$typ"} .= " $acc";
            } else {
                $genbankMap{"$nam$typ"}  =   $acc;
            }

            if (exists($genbankMap{"$nam$ind$typ"})) {
                die "Multiple genbank ACCs for '$nam$ind$typ'.\n";
            } else {
                $genbankMap{"$nam$ind$typ"} = $acc;
            }
        } else {
            die "Failed to match genbank line '$_'\n";
        }
    }
    close(GB);

    return(scalar(keys %genbankMap));
}



#
#  Set GenBank IDs for the short_name found in %data.
#

sub setGenbankIDs ($) {
    my $data       = shift @_;

    my $n = $$data{"short_name"};

    if (exists($genbankMap{"${n}pri"}))    { $$data{"genbank_pri"} = $genbankMap{"${n}pri"}; }
    if (exists($genbankMap{"${n}alt"}))    { $$data{"genbank_alt"} = $genbankMap{"${n}alt"}; }
    if (exists($genbankMap{"${n}mat"}))    { $$data{"genbank_mat"} = $genbankMap{"${n}mat"}; }
    if (exists($genbankMap{"${n}pat"}))    { $$data{"genbank_pat"} = $genbankMap{"${n}pat"}; }
}


return(1);
