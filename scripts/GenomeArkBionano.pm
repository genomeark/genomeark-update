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


sub computeBionanoBases ($$$$$$) {
    my $data      = shift @_;
    my $tind      = shift @_;
    my $size      = shift @_;
    my $file      = shift @_;
    my $missing   = shift @_;
    my $download  = shift @_;

    my $bases     = 0;

    #  Download data.

    if (! -e "downloads/$file") {
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

    if (! -e "downloads/$file") {
        die "Failed to download bionano data for '$file'.\n";
    }

    #  Parse to find the molecule sizes.

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

        if    ($v[0] =~ m/^#/) {
        }
        elsif ($v[0] =~ m/^0/) {
            $nMolecules += 1;
            $sLength    += $v[2];
        }
        elsif ($v[0] =~ m/^1/) {
        }
        elsif ($v[0] =~ m/^Q/) {
        }
        else {
        }
    }
    close(B);

    #  Return bases and reads.

    print STDERR "nM $nMolecules sL $sLength\n";

    return($sLength, $nMolecules);
}
