#
#  Load the existing .md page for a species.
#

sub loadMeta ($$) {
    my $species = shift @_;
    my $meta    = shift @_;

    my  @keys;
    my  @lvls;

    my  $lvl = 0;

    undef %$meta;

    open(MD, "< vgp-metadata/species/$species.yaml");
    while (<MD>) {
        chomp;

        if (m/^(\s*)(\S*):\s*(.*)$/) {
            my $indent = $1;
            my $key    = $2;   $key   =~ s/^\s+//;  $key   =~ s/\s+$//;
            my $value  = $3;   $value =~ s/^\s+//;  $value =~ s/\s+$//;

            my $len    = length($indent);

            #print "\n";
            #print "WORK $len $lvl $key -> $value\n";

            if      ($len  < $lvl) {
                while ($len < $lvl) {
                    #print "pop     ", $keys[-1], " len=$len lvl=$lvl\n";
                    $lvl -= $lvls[-1];
                    pop @keys;
                    pop @lvls;
                }
            }

            if ($len == $lvl) {
                #print "replace ", $keys[-1], "\n";
                pop @keys;
                push @keys, $key;

            } elsif ($len >  $lvl) {
                #print "append  ", $keys[-1], "\n";
                push @keys, $key;
                push @lvls, $len - $lvl;

                $lvl = $len;
            }

            $key = join '.', @keys;

            #print "$key: $value\n";

            $$meta{$key} = $value;
        }
    }

    die "No meta{species.name} found?\n"  if ($$meta{"species.name"} eq "");

    my @n = split '\s+', $$meta{"species.name"};

    die "species.name '", $$meta{"species.name"}, "' has ", scalar(@n), " components, expected 2.\n" if (scalar(@n) < 2);

    return("$n[0]_$n[1]");      #  Return the proper Species_name for this species.
}

