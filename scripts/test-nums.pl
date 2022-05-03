#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use lib "$FindBin::RealBin";

use GenomeArkUtility;

sub test ($) {
    my $s = shift @_;

    my $a = "'" . prettifySize($s) . "'";
    my $b = "'" . prettifyBases($s) . "'";

    printf "%14s - %12s %12s\n", $s, $a, $b;
}


test(0);
test(511);
test(512);
test(513);
test(994);
test(995);
test(996);






foreach my $s (0,
               511, 512, 513,
               994, 995, 996, 999, 1000, 1002, 1004, 1005, 1023, 1024, 1025,
               2047, 2048, 4095, 4096,
               499990, 499994, 499995, 499996,
               512000, 524287, 524288, 542289,
               1047551, 1047552, 1047553, 1048570, 1048572, 1048574, 1048575, 1048576
    ) {
    my $a = "'" . prettifySize($s) . "'";
    my $b = "'" . prettifyBases($s) . "'";#

    printf "%14s - %12s %12s\n", $s, $a, $b;
}
