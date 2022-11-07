package GenomeArkUtility;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(prettifySize prettifyBases);

use strict;
use warnings;

#use Time::Local;
#use List::Util;



#  Format a size in bytes nicely.
#
sub prettifySize ($) {
    my $n = int(shift @_);
    my @d;
    my @p;

    die "prettifySize input size not defined.\n"  if (!defined($n));

    $d[0] = 1024;           $p[0] = 896 * $d[0];
    $d[1] = 1024 * $d[0];   $p[1] = 896 * $d[1];
    $d[2] = 1024 * $d[1];   $p[2] = 896 * $d[2];
    $d[3] = 1024 * $d[2];   $p[3] = 896 * $d[3];

    if    ($n <  128 ) { return(sprintf( "0.1 KiB")); }
    elsif ($n < $p[0]) { return(sprintf("%.1f KiB", ($n / $d[0]))); }
    elsif ($n < $p[1]) { return(sprintf("%.1f MiB", ($n / $d[1]))); }
    elsif ($n < $p[2]) { return(sprintf("%.1f GiB", ($n / $d[2]))); }
    else               { return(sprintf("%.1f TiB", ($n / $d[3]))); }
}


#  Format a size in basepars nicely.
#
#  The four cases will show:
#    up to "### bp"
#    up to "500.## Kbp"
#    up to "500.## Mbp"
#          "###.## Gbp"  (from 0.50 Gbp on)
#
#  Explicitly round and truncate the numbers, since printf() seems to
#  do it differently on different platforms leading to a bit of churn - for example,
#  one platform will show 5.00 and another will show 5.01.
#
sub prettifyBases ($) {
    my $n = shift @_;

    die "prettifyBases input bases not defined.\n"  if (!defined($n));

    if (!defined($n)) { return(undef); }   #  Unset genome size.
    if ($n eq "-")    { return("$n");  }   #  Empty row in N50 table.
    if ($n eq "")     { return(""); }      #  Supplied but unset genome size in yaml.
    if ($n == 0)      { return(""); }      #  Zero genome size.

    if    ($n <      1000) { return("$n  bp");                                              }
    elsif ($n <    499995) { return(sprintf("%.2f Kbp", int($n /       10 + 0.5) / 100.0)); }
    elsif ($n < 499999995) { return(sprintf("%.2f Mbp", int($n /    10000 + 0.5) / 100.0)); }
    else                   { return(sprintf("%.2f Gbp", int($n / 10000000 + 0.5) / 100.0)); }
}



