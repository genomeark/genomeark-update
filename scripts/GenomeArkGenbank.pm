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
    my @projects;

    open(F, "ls projects/genbank.*.map 2>/dev/null |");
    while (<F>) {
        chomp;
        push @projects, $_;
    }
    close(F);

    if (scalar(@projects) == 0) {
        die "No projects/genbank.*.map found.\n";
    }

    print "Loading individual to genbank mappings.\n";

    foreach my $project (@projects) {
        my ($success, $fail) = (0, 0);

        print "  Reading '$project'.\n";

        open(GB, "< $project");
        while (<GB>) {
            chomp;

            #  GCA_009859025.1        bTaeGut2        alt
            if ($_ =~ m/^(GCA_\d+.\d+)\s+([a-z]+[A-Z][a-z]+[A-Z][a-z]+)(\d+)\s+(...)$/) {
                my $acc = $1;   #  GCA_009859025.1
                my $nam = $2;   #  bTaeGut
                my $ind = $3;   #  2
                my $typ = $4;   #  alt

                if (exists($genbankMap{"$nam$typ"})) {
                    $genbankMap{"$nam$typ"} .= " ${nam}${ind}:${acc}";
                } else {
                    $genbankMap{"$nam$typ"}  =  "${nam}${ind}:${acc}";
                }

                $success++;
            } else {
                $fail++;
                #print STDERR "Failed to match genbank line '$_'\n";
            }
        }
        close(GB);

        printf "    Found %5d genbank mappings.\n", $success;
        printf "          %5d genbank mappings were invalid.\n", $fail;
    }

    print "  Removing duplicates.\n";
    
    my $accsBefore = 0;
    my $accsAfter  = 0;

    foreach my $k (keys %genbankMap) {
        my %vals;
        my @vals = split '\s+', $genbankMap{$k};

        $accsBefore += scalar(@vals);

        foreach my $v (@vals) {
            $vals{$v} = 1;
        }

        $accsAfter += scalar(keys %vals);

        $genbankMap{$k} = join ' ', sort keys %vals;
    }

    print "    Removed ", $accsBefore - $accsAfter, " duplicates; $accsAfter accessions remain.\n";

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
