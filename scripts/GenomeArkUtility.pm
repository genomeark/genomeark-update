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
    my @p = ( 0, 1<<10, 1<<20, 1<<30, 1<<40 );

    die "prettifySize undefined input.\n"  if (!defined($n));

    if    ($n < $p[1] - $p[0]) { return(sprintf("1 KiB")) }
    elsif ($n < $p[2] - $p[1]) { return(sprintf("%.0f KiB", ($n/$p[1]))) }  # + (($n % $p[1]) >= ($p[1]/2)))); }
    elsif ($n < $p[3] - $p[2]) { return(sprintf("%.0f MiB", ($n/$p[2]))) }  # + (($n % $p[2]) >= ($p[2]/2)))); }
    elsif ($n < $p[4] - $p[3]) { return(sprintf("%.0f GiB", ($n/$p[3]))) }  # + (($n % $p[3]) >= ($p[3]/2)))); }
    else                       { return(sprintf("%.0f TiB", ($n/$p[4]))) }  # + (($n % $p[4]) >= ($p[4]/2)))); }
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

    die "prettifyBases undefined input.\n"  if (!defined($n));

    if (!defined($n)) { return(undef); }   #  Unset genome size.
    if ($n eq "-")    { return("$n");  }   #  Empty row in N50 table.
    if ($n eq "")     { return(""); }      #  Supplied but unset genome size in yaml.
    if ($n == 0)      { return(""); }      #  Zero genome size.

    if    ($n <      1000) { return("$n  bp");                                              }
    elsif ($n <    499995) { return(sprintf("%.2f Kbp", int($n /       10 + 0.5) / 100.0)); }
    elsif ($n < 499999995) { return(sprintf("%.2f Mbp", int($n /    10000 + 0.5) / 100.0)); }
    else                   { return(sprintf("%.2f Gbp", int($n / 10000000 + 0.5) / 100.0)); }
}



