package GenomeArkMetadata;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(loadSpeciesMetadata saveData loadYAMLasString);

use strict;
use warnings;

use GenomeArkProject;

use YAML::XS;

#
#  Load metadata for a species.
#

sub loadSpeciesMetadata ($$$$) {
    my $data     = shift @_;
    my $species  = shift @_;
    my $template = shift @_;
    my $errors   = shift @_;
    my $mdf;

    undef %$data;

    #  Remove template if real exists.

    if   ((-e "genomeark-metadata/species/$species.yaml") &&
          (-e "genomeark-metadata/species/$species.yaml.template")) {
        unlink "genomeark-metadata/species/$species.yaml.template";
    }

    #  Find metadata file.

    if    (-e "genomeark-metadata/species/$species.yaml") {
        $mdf  = "genomeark-metadata/species/$species.yaml";
    }
    elsif (-e "genomeark-metadata/species/$species.yaml.template") {
        $mdf  = "genomeark-metadata/species/$species.yaml.template";
        $$template{$species}++;
        #push @$errors, "  Species '$species' has only template metadata.\n";
    }
    else {
        die "No metadata found for species '$species'.\n";
    }

    #  Read the metadata.

    #print "Using libYAML ", YAML::XS::LibYAML::libyaml_version(), "\n";
    my $meta = YAML::XS::LoadFile($mdf);

    #  Do a bunch of sanity checks and copy the relevant bits from metadata to our data.

    die "No meta{species.name} found?\n"  if ($meta->{species}->{name} eq "");

    my @n = split '\s+', $meta->{species}->{name};

    die "species.name '", $meta->{species}->{name}, "' has ", scalar(@n), " components, expected at least 2.\n" if (scalar(@n) < 2);

    $$data{"name"}                = $meta->{species}->{name};      #  Name with a space:     'Species name'
    $$data{"name_"}               = join "_", @n;                  #  Name with underscore:  'Species_name'
    $$data{"short_name"}          = $meta->{species}->{short_name};

    if ($$data{"name_"} ne $species) {
        push @$errors, "  Species '$species' is not the same as name '" . $$data{"name_"} . "'\n";
    }

    $$data{"common_name"}         = $meta->{species}->{common_name};
    $$data{"common_name"}         = ""   if (!defined($$data{"common_name"}));

    $$data{"taxon_id"}            = $meta->{species}->{taxon_id};
    $$data{"taxon_id"}            = ""   if (!defined($$data{"taxon_id"}));

    $$data{"genome_size"}         = $meta->{species}->{genome_size};
    $$data{"genome_size"}         = 0   if (!defined($$data{"genome_size"}));
    $$data{"genome_size"}         = 0   if ($$data{"genome_size"} eq "");

    if (($$data{"genome_size"} > 0) && ($$data{"genome_size"} < 1000)) {
        push @$errors, "  Genome size $$data{'genome_size'} for $$data{'name'} ($$data{'common_name'}) is suspiciously low; assuming Gbp is implied.\n";
        $$data{"genome_size"} *= 1000 * 1000 * 1000;
    }

    $$data{"genome_size_method"}  = $meta->{species}->{genome_size_method};
    $$data{"genome_size_method"}  = ""  if (!defined($$data{"genome_size_method"}));

    #  Pull in an exact copy of the metadata file.  We could use YAML::XS::Dump() to
    #  regenerate the yaml, but comments are lost and formatting changes.
    if ($mdf =~ m/y[a]{0,1}ml$/) {
        $$data{"metadata"} = loadYAMLasString($mdf);
    }

    $$data{"data_status"}         = "none";  #<em style=\"color:red\">no data</em>";
    $$data{"assembly_status"}     = "none";  #<em style=\"color:red\">no assembly</em>";

    #$data{"last_raw_data"}       = Do not set here; should only be present if data exists.
    $$data{"last_updated"}        = 0;

    $$data{"data_use_text"}          = $meta->{species}->{data_use_text};
    $$data{"data_use_url"}           = $meta->{species}->{data_use_url};

    $$data{"data_use_project"}       = $meta->{species}->{data_use_project};

    $$data{"point_of_contact"}       = $meta->{species}->{point_of_contact};
    $$data{"point_of_contact_url"}   = $meta->{species}->{point_of_contact_url};
    $$data{"point_of_contact_email"} = $meta->{species}->{point_of_contact_email};

    $$data{"project"}                = $meta->{species}->{project};

    foreach my $p (@{$meta->{species}->{project}}) {
        addProjectToSpecies($p, $$data{name_});
    }

    return($$data{"name_"});      #  Return the proper Species_name for this species.
}


#
#  Read a yaml file into a single string, omitting the '---' lines.
#
sub loadYAMLasString ($) {
    my $mdf = shift @_;
    my $str;

    open(Y, "< $mdf") or die "Failed to open yaml file '$mdf' for reading: $!\n";
    while (<Y>) {
        $str .= "$_"   if ($_ !~ "^---");
    }
    close(Y);

    return($str);
}


#
#  Write and read data from our markdown page.
#

sub saveData ($$) {
    my $data = shift @_;
    my $file = shift @_;

    print "  Write to '$file'\n";

    open(MD, "> $file") or die "Failed to open '$file' for write: $!\n";
    print MD YAML::XS::Dump($data);
    print MD "---\n";
    close(MD);
}

