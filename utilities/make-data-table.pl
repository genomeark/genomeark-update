#!/usr/bin/env perl
use strict;

my %tyn;

$tyn{"pacbio"}          = "PacBio CLR";
$tyn{"pacbiohifi_fqgz"} = "PacBio HiFi (.fq, Q&ge;20)";
$tyn{"pacbiohifi_bam"}  = "PacBio HiFi (.bam, Q&ge;0)";
$tyn{"pacbiohifi_clr"}  = "PacBio HiFi (subreads)";
$tyn{"ont"}             = "Oxford Nanopore Simplex";
$tyn{"ontduplex"}       = "Oxford Nanopore Duplex";
$tyn{"10x"}             = "10x";
$tyn{"bionano"}         = "Bionano";
$tyn{"arima"}           = "Arima";
$tyn{"dovetail"}        = "Dovetail Genomics";
$tyn{"phase"}           = "Phase";
$tyn{"illumina"}        = "Illumina";

print "{% comment %}IMPORT FROM MAKE-DATA-TABLE BEGINS{% endcomment %}\n";
print "\n";
print "<table class=\"raw-data-table\">\n";
print "<thead><tr><th class=\"label\">Datatype</th>\n";
print "           <th class=\"label\">Individual</th>\n";
print "           <th>Bases</th>\n";
print "           <th>Coverage</th>\n";
print "           <th>Bytes</th>\n";
print "           <th>Access</th>\n";
print "       </tr></thead>\n";
print "<tbody>\n";
foreach my $ty (qw(pacbio pacbiohifi_fqgz pacbiohifi_bam ont ontduplex 10x bionano arima dovetail phase illumina)) {
    print "\n";
    print "{% assign rs = 0 %}\n";
    print "\n";
    foreach my $id (qw(1 2 3 4 5 6)) {
        print "{% if page.data_${ty}-${id}_bytes %}{% assign rs = rs | plus: 1 %}{% endif %}\n";
    }
    print "\n";
    print "{% if rs > 0 %}\n";
    print "  {% assign first = true %}\n";
    print "\n";
    foreach my $id (qw(1 2 3 4 5 6)) {
        print "  {% if page.data_${ty}-${id}_bytes %}\n";
        print "    {% if first %}\n";
        print "      {% assign first = false %}\n";
        print "      <tr class=\"raw-data-toprow\">\n";
        print "          <td rowspan=\"{{ rs }}\" class=\"label\">$tyn{$ty}</td>\n";
        print "          <td class=\"label\">{{ page.short_name }}${id}</td>\n";
        print "          <td>{{ page.data_${ty}-${id}_bases   }}</td>\n";
        print "          <td>{{ page.data_${ty}-${id}_coverage}}</td>\n";
        print "          <td>{{ page.data_${ty}-${id}_bytes   }}</td>\n";
        print "          <td><a title=\"Select files to download from S3.\" href=\"{{ page.data_${ty}-${id}_s3url }}\"><i class=\"svg-icon2 download\"></i></a>\n";
        print "              <a title=\"Preview files using 42basepairs.com.\" href=\"{{ page.data_${ty}-${id}_s3gui }}\"><i class=\"svg-icon2 view\"></i></a></td></tr>\n";
        print "    {% else %}\n";
        print "      <tr>\n";
        print "          <td class=\"label\">{{ page.short_name }}${id}</td>\n";
        print "          <td>{{ page.data_${ty}-${id}_bases   }}</td>\n";
        print "          <td>{{ page.data_${ty}-${id}_coverage}}</td>\n";
        print "          <td>{{ page.data_${ty}-${id}_bytes   }}</td>\n";
        print "          <td><a title=\"Select files to download from S3.\" href=\"{{ page.data_${ty}-${id}_s3url }}\"><i class=\"svg-icon2 download\"></i></a>\n";
        print "              <a title=\"Preview files using 42basepairs.com.\" href=\"{{ page.data_${ty}-${id}_s3gui }}\"><i class=\"svg-icon2 view\"></i></a></td></tr>\n";
        print "    {% endif %}\n";
        print "  {% endif %}\n";
        print "\n";
    }
    print "{% endif %}\n";
    print "\n";
}
print "</tbody>\n";
print "</table>\n";
print "\n";
print "<div class=\"genome-size\">\n";
print "Bases and Coverage are approximate.<br>\n";
print "{% if page.genome_size_display %}\n";
print "Coverage based on genome size {{ page.genome_size_display }}.<br>\n";
print "{% else %}\n";
print "No genome size estimate available.<br>\n";
print "{% endif %}\n";
print "Last upload on {{ page.last_raw_data | date: \"%-d %B %Y\" }}.\n";
print "</div>\n";
print "\n";
print "{% comment %}IMPORT FROM MAKE-DATA-TABLE ENDS{% endcomment %}\n";
