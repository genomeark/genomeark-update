package GenomeArkBionano;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(computeBionanoBases);

use strict;
use warnings;

use GenomeArkUtility;

#use Time::Local;
use List::Util qw(min max);
use File::Basename;

my $aws          = "aws --no-sign-request";


sub computeBionanoBases ($$$$$$$$) {
    my $data      = shift @_;
    my $files     = shift @_;   #  An array!
    my $dataBytes = shift @_;
    my $dataQuant = shift @_;
    my $dataBases = shift @_;
    my $dataReads = shift @_;
    my $missing   = shift @_;
    my $download  = shift @_;

    my $bases     = 0;

    foreach my $sizefile (@$files) {
        my ($size, $file, $name);

        #  Process only the cmap and cmap.gz files, ignore the bnx files.

        ($size, $file, $name) = ($1, "$2$3", $2)   if ($sizefile =~ m/^(\d+)\s(.*)(\.cmap\.gz)$/);
        ($size, $file, $name) = ($1, "$2$3", $2)   if ($sizefile =~ m/^(\d+)\s(.*)(\.cmap)$/);

        next   if (! defined($file));

        if ((! -e "downloads/$file") &&
            (! -e "$name.summary")) {

            if ($download == 0) {
                printf "  SKIP      -    all out of %6.3f MB - s3://genomeark/$file\n", $size / 1024 / 1024;

                $$missing{ $$data{"name_"} } .= "bionano ";
            }
            else {
                printf "  FETCH     -    all out of %6.3f MB - s3://genomeark/$file\n", $size / 1024 / 1024;

                system("mkdir -p downloads/" . dirname($file));
                system("aws --no-progress --no-sign-request s3 cp s3://genomeark/$file downloads/$file > downloads/$file.err 2>&1");
            }
        }

        #  Parse to find the molecule sizes.

        if ((-e "downloads/$file") &&
            (! exists($$dataBases{$name}))) {
            print "  PARSE\n";

            my $bgn = 9999999999;
            my $end = 0;

            my $nMolecules = 0;
            my $sLength    = 0;

            if ($file =~ m/gz$/) {
                open(B, "gzip -dc downloads/$file |") or die;
            } else {
                open(B, "< downloads/$file") or die;
            }
            while (<B>) {
                my @v = split '\s+', $_;

                next   if ($v[0] =~ m/^#/);

                if ($v[4] == 1) {
                    $bgn = min($v[5], $bgn);
                    $end = max($v[5], $end);
                }
                else {
                    $nMolecules += 1;
                    $sLength    += $end - $bgn;

                    #print "$v[0] $v[1] $end $bgn ", $end - $bgn, "\n";

                    $bgn = 9999999999;
                    $end = 0;
                }
            }
            close(B);

            $$dataBytes{$name} = 0;
            $$dataQuant{$name} = 0;
            $$dataBases{$name} = $sLength;
            $$dataReads{$name} = $nMolecules;

            #system("mkdir -p " . dirname($name));

            #open(S, "> $name.summary") or die "Failed to open '$name.summary' for writing: $!\n";;
            #print S "Molecules: $nMolecules\n";
            #print S "Length:    $sLength\n";
            #close(S);
        }

        #if (-e "$name.summary") {
        #    open(S, "< $name.summary") or die "Failed to open '$name.summary' for reading: $!\n";
        #    while (<S>) {
        #        $bases += $1   if (m/Length:\s+(.*)$/);
        #    }
        #    close(S);
        #}
    }

    #return($bases);
}
