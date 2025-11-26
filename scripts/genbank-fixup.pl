#!/usr/bin/env perl

#  Use genbank-fetch.sh to pull the genbank XML and parse out a few fields.
#  Use genbank-fixup.pl to clean up the result and write genbank.map.

use strict;

my %found;

while (<STDIN>) {
    chomp;

    if ($_ =~ m/^(GC[AF]_\d+.\d+)\s+(\S+)\s+(.*)$/) {
        my $acc = $1;  #  Correct as is.
        my $lab = $2;  #  Minor mods needed.
        my $dsc = $3;  #  A mess.

        my $nam = "XXX";  #  Properly formatted VGP individual name
        my $suf = "XXX";  #  Rest of the name label
        my $typ = "XXX";  #  Type of this assembly (pri alt mat pat)

        #  Fix mistakes.

        if (($acc eq "GCA_022578405.1") && ($lab eq "fEmbJac1.0.p"))   { $lab = "fEmbJac1.0.a"; }

        #  Fix up obvious mistakes in the name.

        if ($lab =~ m/^g(adMor.*)$/)       { $lab = "fG$1"; }
        if ($lab =~ m/^mCalJa(1.2.*)$/)    { $lab = "mCalJac$1"; }

        #  Ignore bizarre stuff.

        if ($lab eq "DD_fGirMul_XY1")      { next; }
        if ($lab eq "UOA_Angus_1")         { next; }
        if ($lab eq "UOA_Brahman_1")       { next; }
        if ($lab eq "ZJU1.0")              { next; }
        if ($lab eq "mRatBN7.1")           { next; }
        if ($lab eq "mRatBN7.2")           { next; }

        #  Find the name and suffix, or emit a dire warning.

        if      ($lab =~ m/^([a-z]+[A-Z][a-z]+[A-Z][a-z]+[0-9]+)(.*)$/) {
            $nam = $1;
            $suf = $2;
        } elsif ($lab =~ m/^ASM/) {
        } else {
            printf STDERR "WARN: %-15s %-10s %-6s lab %-18s dsc %s\n", $acc, $nam, $typ, $lab, $dsc;
            next;
        }

        #  Use the suffix of the name to decide the type of this assembly.

        if    ($suf =~ m/mat\.[A-Z]/)   { $typ = "mgd"; }
        elsif ($suf =~ m/pat\.[A-Z]/)   { $typ = "mgd"; }

        elsif ($suf =~ m/pri/)          { $typ = "pri"; }
        elsif ($suf =~ m/alt/)          { $typ = "alt"; }

        elsif ($suf =~ m/mat/)          { $typ = "mat"; }
        elsif ($suf =~ m/pat/)          { $typ = "pat"; }

        elsif ($suf =~ m/\.p$/)         { $typ = "pri"; }
        elsif ($suf =~ m/\.a$/)         { $typ = "alt"; }

        #elsif ($suf =~ m/mat/)          { $typ = "mat"; }
        #elsif ($suf =~ m/pat/)          { $typ = "pat"; }

        #  If not decided, use the description.
        #
        #  There seem to only be three forms here:
        #    % cut -f 3 genbank.map.raw | sort | uniq -c
        #    258 alternate-pseudohaplotype
        #     57 haploid
        #    229 haploid (principal pseudohaplotype of diploid)

        if ($typ eq "XXX") {
            if    ($dsc =~ m/principal/)  { $typ = "pri"; }
            elsif ($dsc =~ m/alter/)      { $typ = "alt"; }

            elsif ($dsc =~ m/maternal/)   { $typ = "mat"; }
            elsif ($dsc =~ m/paternal/)   { $typ = "pat"; }

            elsif ($dsc =~ m/haploid/)    { $typ = "pri"; }
        }

        if ($typ eq "XXX") {
            printf STDERR "UNKN: %-15s %-10s %-6s lab %-18s dsc %s\n", $acc, $nam, $typ, $lab, $dsc;
        }

        #  Hardcode some ugly ones.

        #if ($lab eq "bTaeGut2.mat.v3")    { $typ = "mat"; }
        #if ($lab eq "bTaeGut2pat")        { $typ = "pat"; }

        #if ($lab eq "bTaeGut2.pri.v2")    { $typ = "pri"; }
        #if ($lab eq "bTaeGut2.p.v1.alt")  { $typ = "alt"; }

        #if ($lab eq "eAstRub1.3")   { $typ = "pri"; }
        #if ($lab eq "fAstCal1.2")   { $typ = "pri"; }
        #if ($lab eq "fBetSpl5.2")   { $typ = "pri"; }
        #if ($lab eq "fSalTru1.1")   { $typ = "pri"; }

        next   if ($lab eq "mRhiFer1_v1.h");
        next   if ($lab eq "mRhiFer1_v1.p");

        #  Emit output or fail.
        if ($typ eq "XXX") {
            printf STDERR "FAIL: %-15s %-10s %-6s lab %-18s dsc %s\n", $acc, $nam, $typ, $lab, $dsc;
        } else {
            #printf STDERR "PASS: %-15s %-10s %-6s lab %-18s dsc %s\n", $acc, $nam, $typ, $lab, $dsc;
            print STDOUT "$acc\t$nam\t$typ\n";
        }

        if (exists($found{"$nam$typ"})) {
            print STDERR "DUPL: acc '$acc' lab '$lab'\n";
        }
        $found{"$nam$typ"} = 1;
    }

    else {
        die "Failed to match '$_'\n";
    }
}

exit(0);
