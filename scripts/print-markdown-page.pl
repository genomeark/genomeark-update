

sub printData ($$$) {
    my $species = shift @_;
    my $data    = shift @_;

    my @keys = keys %$data;

    @keys = sort @keys;

    foreach my $key (@keys) {
        next if ($key eq ":");

        chomp $$data{$key};   #  Remove the newline on multi-line data.

        if ($$data{$key} ne "") {
            print "$key: $$data{$key}\n";
        } else {
            print "$key:\n";
        }
    }
    print $$data{":"};
}

