package GenomeArkMetadata;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(loadSpeciesMetadata saveData loadData);

use strict;
use warnings;

use GenomeArkUtility;

#use Time::Local;
#use List::Util;


#
#  Load metadata from our cache, or fetch it from the bucket if it is newer.
#

sub loadSpeciesMetadata ($$$$) {
    my $data    = shift @_;
    my $species = shift @_;
    my $meta    = shift @_;
    my $errors  = shift @_;

    my  @keys;
    my  @lvls;

    my  $lvl = 0;

    undef %$data;
    undef %$meta;

    if (! -e "downloads/species/$species/metadata.yaml") {
        system("mkdir -p downloads/species/$species");
        system("aws s3 cp s3://genomeark/species/$species/metadata.yaml downloads/species/$species/metadata.yaml");
    }

    open(MD, "< downloads/species/$species/metadata.yaml") or die;
    while (<MD>) {
        chomp;

        if (m/^(\s*)(\S*):\s*(.*)$/) {
            my $indent = $1;
            my $key    = $2;   $key   =~ s/^\s+//;  $key   =~ s/\s+$//;
            my $value  = $3;   $value =~ s/^\s+//;  $value =~ s/\s+$//;

            my $len    = length($indent);

            if      ($len  < $lvl) {
                while ($len < $lvl) {
                    $lvl -= $lvls[-1];
                    pop @keys;
                    pop @lvls;
                }
            }

            if ($len == $lvl) {
                pop @keys;
                push @keys, $key;

            } elsif ($len >  $lvl) {
                push @keys, $key;
                push @lvls, $len - $lvl;

                $lvl = $len;
            }

            $key = join '.', @keys;

            $$meta{$key} = $value;
        }
    }
    close(MD);

    die "No meta{species.name} found?\n"  if ($$meta{"species.name"} eq "");

    my @n = split '\s+', $$meta{"species.name"};

    die "species.name '", $$meta{"species.name"}, "' has ", scalar(@n), " components, expected 2.\n" if (scalar(@n) < 2);

    $$data{"name"}                = $$meta{"species.name"};         #  Name with a space:     'Species name'
    $$data{"name_"}               = "$n[0]_$n[1]";                  #  Name with underscore:  'Species_name'
    $$data{"short_name"}          = $$meta{"species.short_name"};

    if ($$data{"name_"} ne $species) {
        push @$errors, "  Species '$species' is not the same as name '" . $$data{"name_"} . "'\n";
    }

    $$data{"common_name"}         = $$meta{"species.common_name"};
    $$data{"common_name"}         = ""   if (!defined($$data{"common_name"}));

    $$data{"taxon_id"}            = $$meta{"species.taxon_id"};
    $$data{"taxon_id"}            = ""   if (!defined($$data{"taxon_id"}));

    $$data{"genome_size"}         = $$meta{"species.genome_size"};
    $$data{"genome_size"}         = 0   if (!defined($$data{"genome_size"}));
    $$data{"genome_size"}         = 0   if ($$data{"genome_size"} eq "");

    if (($$data{"genome_size"} > 0) && ($$data{"genome_size"} < 1000)) {
        push @$errors, "  Genome size $$data{'genome_size'} for $$data{'name'} ($$data{'common_name'}) is suspiciously low; assuming Gbp is implied.\n";
        $$data{"genome_size"} *= 1000 * 1000 * 1000;
    }

    $$data{"genome_size_method"}  = $$meta{"species.genome_size_method"};
    $$data{"genome_size_method"}  = ""  if (!defined($$data{"genome_size_method"}));

    $$data{"data_status"}         = "none";  #<em style=\"color:red\">no data</em>";
    $$data{"assembly_status"}     = "none";  #<em style=\"color:red\">no assembly</em>";

    #$data{"last_raw_data"}       = Do not set here; should only be present if data exists.
    $$data{"last_updated"}        = 0;

    return($$data{"name_"});      #  Return the proper Species_name for this species.
}



#
#  Write and read data from our markdown page.
#

sub saveData ($$) {
    my $data = shift @_;
    my $file = shift @_;

    print "  Write to '$file'\n";
    open(MD, "> $file") or die "Failed to open '$file' for write: $!\n";

    print MD "---\n";
    foreach my $key (sort keys %$data) {
        next if ($key eq ":");

        die "undef data{$key}\n"  if (!defined($$data{$key}));

        if    (ref $$data{$key} eq "ARRAY") {
            print MD "$key:\n";
            foreach my $val (@$data{$key}) {
                print MD " - $val\n";
            }
        }
        elsif ($$data{$key} ne "") {
            chomp $$data{$key};
            print MD "$key: $$data{$key}\n";
        }
    }

    print MD "---\n";
    print MD $$data{":"}   if (exists($$data{":"}));

    close(MD);
}


sub loadData ($$) {
    my $species = shift @_;
    my $data    = shift @_;
    my $key     = undef;
    my $lvl     = 0;
    my @keys;
    my @lvls;

    undef %$data;   #  Forget everything we know about some other species.

    open(MD, "< ../_genomeark/$species.md") or die;
    while (<MD>) {
        chomp;

        if    (m!^---$!) {                        #  Match the YAML delimiters.
            $key    = undef;
        }
        elsif (m!(\w+):\s*\|$!) {                 #  Match any multi-line pre-formatted parameters.
            $key = $1;
            $$data{$key} = "|\n";
        }
        elsif (defined($key) && (m!^\s\s\S!)) {   #  Add data to multi-line parameters.
            $$data{$key} .= "$_\n";
        }
        elsif (m!(\w+):\s*(.*)$!) {               #  Match any key-value pairs.
            $key       = undef;
            $$data{$1} = $2;
        }
        else {                                    #  Match any extra stuff in the file.
            $$data{":"} .= "$_\n";
        }
    }
    close(MD);
}

