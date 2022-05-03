






    #
    #  Finalize the genomic_data by adding to %data.
    #

    $seqFiles{"dovetail"} /= 2;
    $seqFiles{"illumina"} /= 2;
    $seqFiles{"phase"}    /= 2;

    foreach my $type (qw(10x arima bionano dovetail illumina pbsubreads pbhifi pbhifisub phase)) {
        if ($seqBytes{$type} > 0) {
            $data{"data_${type}_bytes"}    = sprintf("%.3f GB", $seqBytes{$type} / 1024 / 1024 / 1024);
            $data{"data_${type}_coverage"} = "N/A";
            $data{"data_${type}_bases"}    = "N/A";
            $data{"data_${type}_files"}    = $seqFiles{$type};
        }
    }

    #  Update genome size if not set.

    if ($data{"genome_size"} == 0)   { $data{"genome_size"} = $data{"pri1length"}; }
    if ($data{"genome_size"} == 0)   { $data{"genome_size"} = $data{"pri2length"}; }
    if ($data{"genome_size"} == 0)   { $data{"genome_size"} = $data{"pri3length"}; }
    if ($data{"genome_size"} == 0)   { $data{"genome_size"} = $data{"pri4length"}; }

    $data{"genome_size_display"} = prettifyBases($data{"genome_size"});

    #
    #  Estimate the number of bases in all the raw data files.
    #    BIONANO is totally different and not supported.
    #    HiFi subreads aren't computed either, not useful and big.
    #

    $data{"data_10x_scale"}        = estimateRawDataScaling($data{"name_"}, "10x")          if ($data{"data_10x_scale"}        == 0);
    $data{"data_arima_scale"}      = estimateRawDataScaling($data{"name_"}, "arima")        if ($data{"data_arima_scale"}      == 0);
    $data{"bionano_scale"}         = 0;
    $data{"data_dovetail_scale"}   = estimateRawDataScaling($data{"name_"}, "dovetail")     if ($data{"data_dovetail_scale"}   == 0);
    $data{"data_illumina_scale"}   = estimateRawDataScaling($data{"name_"}, "illumina")     if ($data{"data_illumina_scale"}   == 0);
    $data{"data_pbsubreads_scale"} = estimateRawDataScaling($data{"name_"}, "pbsubreads")   if ($data{"data_pbsubreads_scale"} == 0);
    $data{"data_pbhifi_scale"}     = estimateRawDataScaling($data{"name_"}, "pbhifi")       if ($data{"data_pbhifi_scale"}     == 0);
    #data{"data_pbhifisub_scale"}  = estimateRawDataScaling($data{"name_"}, "pbhifisub")    if ($data{"data_pbhifisub_scale"}  == 0);
    $data{"data_pbhifisub_scale"}  = 0;
    $data{"data_phase_scale"}      = estimateRawDataScaling($data{"name_"}, "phase")        if ($data{"data_phase_scale"}      == 0);

    #  Figure out how much and what types of data exist.

    my $dataPac = 0;
    my $data10x = 0;
    my $dataHIC = 0;
    my $dataBio = 0;

    sub uniquifyStringArray ($) {    #  Split a \0 delimited string into an
        my @a = split '\0', $_[0];   #  array of unique elements.
        my %a;

        foreach my $a (@a) {
            $a{$a}++;
        }

        return(sort keys %a);
    }


    if (($seqBytes{"10x"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"10x"})) {
            $data{"data_10x_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/10x/ .<br>";
        }

        if (($data{"data_10x_scale"} > 0) && ($data{"genome_size"} > 0)) {
            $data{"data_10x_coverage"}        = sprintf("%.2fx", $seqBytes{"10x"} * $data{"data_10x_scale"} / $data{"genome_size"});
            $data{"data_10x_bases"}           = prettifyBases($seqBytes{"10x"} * $data{"data_10x_scale"});
        }
        $data10x++;
    }

    if (($seqBytes{"arima"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"arima"})) {
            $data{"data_arima_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/arima/ .<br>";
        }

        if (($data{"data_arima_scale"} > 0) && ($data{"genome_size"} > 0)) {
            $data{"data_arima_coverage"}      = sprintf("%.2fx", $seqBytes{"arima"} * $data{"data_arima_scale"} / $data{"genome_size"});
            $data{"data_arima_bases"}         = prettifyBases($seqBytes{"arima"} * $data{"data_arima_scale"});
        }
        $dataHIC++;
    }

    if (($seqBytes{"bionano"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"bionano"})) {
            $data{"data_bionano_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/bionano/ .<br>";
        }

        if ($data{"genome_size"} > 0) {    #  No scaling for bionano!
            my $b = computeBionanoBases($data{"name_"});
            $data{"data_bionano_coverage"}    = sprintf("%.2fx", $b / $data{"genome_size"});
            $data{"data_bionano_bases"}       = prettifyBases($b);
        }
        $dataBio++;
    }

    if (($seqBytes{"dovetail"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"dovetail"})) {
            $data{"data_dovetail_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/dovetail/ .<br>";
        }

        if (($data{"data_dovetail_scale"} > 0) && ($data{"genome_size"} > 0)) {
            $data{"data_dovetail_coverage"}   = sprintf("%.2fx", $seqBytes{"dovetail"} * $data{"data_dovetail_scale"} / $data{"genome_size"});
            $data{"data_dovetail_bases"}      = prettifyBases($seqBytes{"dovetail"} * $data{"data_dovetail_scale"});
        }
        $dataHIC++;
    }

    if (($seqBytes{"illumina"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"illumina"})) {
            $data{"data_illumina_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/illumina/ .<br>";
        }

        if (($data{"data_illumina_scale"} > 0) && ($data{"genome_size"} > 0)) {
            $data{"data_illumina_coverage"}   = sprintf("%.2fx", $seqBytes{"illumina"} * $data{"data_illumina_scale"} / $data{"genome_size"});
            $data{"data_illumina_bases"}      = prettifyBases($seqBytes{"illumina"} * $data{"data_illumina_scale"});
        }
    }

    if (($seqBytes{"pbsubreads"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"pbsubreads"})) {
            $data{"data_pbsubreads_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/pacbio/ . --exclude \"*ccs*bam*\"<br>";
        }

        if (($data{"data_pbsubreads_scale"} > 0) && ($data{"genome_size"} > 0)) {
            $data{"data_pbsubreads_coverage"} = sprintf("%.2fx", $seqBytes{"pbsubreads"} * $data{"data_pbsubreads_scale"} / $data{"genome_size"});
            $data{"data_pbsubreads_bases"}    = prettifyBases($seqBytes{"pbsubreads"} * $data{"data_pbsubreads_scale"});
        }
        $dataPac++;
    }

    if (($seqBytes{"pbhifi"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"pbhifi"})) {
            $data{"data_pbhifi_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/pacbio/ . --exclude \"*subreads.bam*\"<br>";
        }

        if (($data{"data_pbhifi_scale"}) && ($data{"genome_size"} > 0)) {
            $data{"data_pbhifi_coverage"} = sprintf("%.2fx", $seqBytes{"pbhifi"} * $data{"data_pbhifi_scale"} / $data{"genome_size"});
            $data{"data_pbhifi_bases"}    = prettifyBases($seqBytes{"pbhifi"} * $data{"data_pbhifi_scale"});
        }
        $dataPac++;
    }

    if (($seqBytes{"pbhifisub"} > 0)) {
        #foreach my $k (uniquifyStringArray($seqIndiv{"pbhifisub"})) {
        #    $data{"data_pbhifisub_links"} .= "aws s3 ......<br>";   #  NOT listed on the page.
        #}

        if (($data{"data_pbhifisub_scale"} > 0) && ($data{"genome_size"} > 0)) {
            $data{"data_pbhifisub_coverage"} = sprintf("%.2fx", $seqBytes{"pbhifisub"} * $data{"data_pbhifisub_scale"} / $data{"genome_size"});
            $data{"data_pbhifisub_bases"}    = prettifyBases($seqBytes{"pbhifisub"} * $data{"data_pbhifisub_scale"});
        }
        $dataPac++;
    }

    if (($seqBytes{"phase"} > 0)) {
        foreach my $k (uniquifyStringArray($seqIndiv{"phase"})) {
            $data{"data_phase_links"} .= "aws s3 --no-sign-request sync s3://genomeark/species/$k/genomic_data/phase/ .<br>";
        }

        if (($data{"data_phase_scale"} > 0) && ($data{"genome_size"} > 0)) {
            $data{"data_phase_coverage"}      = sprintf("%.2fx", $seqBytes{"phase"} * $data{"data_phase_scale"} / $data{"genome_size"});
            $data{"data_phase_bases"}         = prettifyBases($seqBytes{"phase"} * $data{"data_phase_scale"});
        }
        $dataHIC++;
    }

    #  Set the classification.
    #  If no assembly, put it under either "some data" or "all data".
    #  Otherwise, put it under the assembly classification.

    if (($dataPac > 0) ||
        ($data10x > 0) ||
        ($dataHIC > 0) ||
        ($dataBio > 0)) {
        $data{"data_status"} = "some";
    }

    if (($dataPac > 0) &&
        ($data10x > 0) &&
        ($dataHIC > 0) &&
        ($dataBio > 0)) {
        $data{"data_status"} = "all";
    }

    #  Create symlinks to the categories.
    #
    #  Name changes should change:
    #  1)  the label used in this script ("low-quality-draft")
    #  2)  the directory name in the symlink
    #  3)  the actual directories (both "_genomeark-low-quality-draft" and "genomeark-low-quality-draft")
    #  4)  the index.html in genomeark-low-quality-draft/  (both the title and the collection name)
    #  5)  _config.yml (in two places) and index.html in the root directory.

    my $name = $data{"name_"};

    unlink("../_genomeark-curated-assembly/$name.md");
    unlink("../_genomeark-high-quality-draft-assembly/$name.md");
    unlink("../_genomeark-low-quality-draft-assembly/$name.md");
    unlink("../_genomeark-complete-data/$name.md");
    unlink("../_genomeark-incomplete-data/$name.md");

    if    ($data{"assembly_status"} eq "curated") {
        system("mkdir -p ../_genomeark-curated-assembly");
        symlink("../_genomeark/$name.md", "../_genomeark-curated-assembly/$name.md") or die "Failed to make symlink for curated assembly: $!.\n";
    }

    elsif ($data{"assembly_status"} eq "high-quality-draft") {
        system("mkdir -p ../_genomeark-high-quality-draft-assembly");
        symlink("../_genomeark/$name.md", "../_genomeark-high-quality-draft-assembly/$name.md") or die "Failed to make symlink for high-quality-draft assembly: $!.\n";
    }

    elsif ($data{"assembly_status"} eq "low-quality-draft-assembly") {
        system("mkdir -p ../_genomeark-low-quality-draft");
        symlink("../_genomeark/$name.md", "../_genomeark-low-quality-draft-assembly/$name.md") or die "Failed to make symlink for low-quaity-draft assembly: $!.\n";
    }

    elsif ($data{"data_status"} eq "all") {
        system("mkdir -p ../_genomeark-complete-data");
        symlink("../_genomeark/$name.md", "../_genomeark-complete-data/$name.md") or die "Failed to make symlink for complete data: $!.\n";
    }

    elsif ($data{"data_status"} eq "some") {
        system("mkdir -p ../_genomeark-incomplete-data");
        symlink("../_genomeark/$name.md", "../_genomeark-incomplete-data/$name.md") or die "Failed to make symlink for incomplete data: $!.\n";
    }

    else {
        system("mkdir -p ../_genomeark-incomplete-data");
        symlink("../_genomeark/$name.md", "../_genomeark-incomplete-data/$name.md") or die "Failed to make symlink for incomplete data (catch all): $!.\n";
    }

    #  And reset the classifications to strings we can use in the display.

    $data{"data_status"}     = "<em style=\"color:red\">no data</em>"                          if ($data{"data_status"}) eq "none";
    $data{"data_status"}     = "<em style=\"color:orange\">some data</em>"                     if ($data{"data_status"}) eq "some";
    $data{"data_status"}     = "<em style=\"color:green\">all data</em>"                       if ($data{"data_status"}) eq "all";

    $data{"assembly_status"} = "<em style=\"color:red\">no assembly</em>"                      if ($data{"assembly_status"} eq "none");
    $data{"assembly_status"} = "<em style=\"color:red\">low-quality draft assembly</em>"       if ($data{"assembly_status"} eq "low-quality-draft");
    $data{"assembly_status"} = "<em style=\"color:orange\">high-quality draft assembly</em>"   if ($data{"assembly_status"} eq "high-quality-draft");
    $data{"assembly_status"} = "<em style=\"color:green\">curated assembly</em>"               if ($data{"assembly_status"} eq "curated");

    #  If no date set -- no raw data, no assemblies, no anything -- set it to the
    #  date of this update (genomeark.ls's date).

    if ($data{"last_updated"} == 0) {
        print "WARNING: no files; last_updated set to 'now'.\n";
        $data{"last_updated"} = $lastUpdate;
    }

    #  Done.  Write the output.

    if ($name ne $species) {
        print "WARNING: species '$species' is not the same as name '$name' (using name instead)\n";
    }

    #printData($species, \%data);
    saveData($name, \%data);

