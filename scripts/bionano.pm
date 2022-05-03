sub computeBionanoBases ($) {
    my $name    = shift @_;
    my @files;

    for (my $ii=0; $ii<scalar(@speciesSizes); $ii++) {
        my $filesize = $speciesSizes[$ii];
        my $filename = $speciesFiles[$ii];

        push @files, "$filesize\0$filename"   if ($filename =~ m!$name.*genomic_data/bionano/.*bnx.gz$!);
    }

    #foreach my $f (@files) {
    #    printf "BIONANO %12d - '%s'\n", split '\0', $f;
    #}

    #  Download all the cmap files.

    foreach my $f (@files) {
        my ($size, $name) = split '\0', $f;

        if ((! -e "downloads/$name") &&
            (! -e "$name.summary")) {
            printf "FETCH bionano - size %6.3f GB\n", $size / 1024 / 1024 / 1024;
            printf "  aws --no-progress --no-sign-request s3 cp \\\n";
            printf "    s3://genomeark/$name \\\n";
            printf "         downloads/$name\n";
            printf "\n";

            system("mkdir -p downloads/$name");
            system("rmdir    downloads/$name");
            system("aws --no-progress --no-sign-request s3 cp s3://genomeark/$name downloads/$name")  if ($downloadRAW);
        }
    }

    #  Parse each bnx to find the molecule sizes

    foreach my $f (@files) {
        my ($size, $name) = split '\0', $f;

        next  if (! -e "downloads/$name");
        next  if (  -e "$name.summary");

        print "PARSING $name\n";

        my $nMolecules = 0;
        my $sLength    = 0;

        system("mkdir -p $name");
        system("rmdir    $name");

        open(B, "gzip -dc downloads/$name |");
        while (<B>) {
            if (m/^0\s/) {
                my @v = split '\s+', $_;

                $nMolecules += 1;
                $sLength    += $v[2];
            }
        }
        close(B);

        open(S, "> $name.summary") or die "Failed to open '$name.summary' for writing: $!\n";;
        print S "Molecules: $nMolecules\n";
        print S "Length:    $sLength\n";
        close(S);
    }

    #  Read the summaries to find the bases.

    my $bases = 0;

    foreach my $f (@files) {
        my ($size, $name) = split '\0', $f;

        if (-e "$name.summary") {
            open(S, "< $name.summary") or die "Failed to open '$name.summary' for reading: $!\n";
            while (<S>) {
                $bases += $1   if (m/Length:\s+(.*)$/);
            }
            close(S);
        }
    }

    return($bases);
}
