#!/usr/bin/env perl
#
#############################################################################
 #
 #  Part of GenomeArk, a collection of genomic raw data and assemblies.
 #
 #  Except as indicated otherwise, this is a 'United States Government
 #  Work', and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 #
 ##

use strict;

#
#  Expects an S3 'ls --recursive' list on input and outputs a list of files
#  that are useful for displaying statistics on genomeark.github.io.
#
#  Takes about half a minute on the bucket as of June 2022.
#
#  Also generates a list of the species with a metadata file,
#  and a report of the species without a metadata file.
#

#  time perl ../scripts/filter-aws-bucket-list.pl --raw genomeark.ls.raw.20220601 --filtered genomeark.ls --species-list species-list

my %speciesList;   #  Count of the number of filtered files for each species.
my %speciesMeta;   #  Defined if 'Genus_species/metadata.yaml' is found.

my $lsRaw = $ARGV[1];
my $lsFilt = $ARGV[2];
my $spList = undef;

while (scalar(@ARGV) > 0) {
    my $opt = shift @ARGV;

    if    ($opt eq "--raw") {
        $lsRaw = shift @ARGV;
    }
    elsif ($opt eq "--filtered") {
        $lsFilt = shift @ARGV;
    }
    elsif ($opt eq "--species-list") {
        $spList = shift @ARGV;
    }
    else {
        die "Unknown option '$opt'\n";
    }
}

open(RAW,  "< $lsRaw")  or die "Failed to open --raw '$lsRaw' for reading: $!\n";
open(FILT, "> $lsFilt") or die "Failed to open --filtered '$lsFilt' for writing: $!\n";
open(SPLI, "> $spList") or die "Failed to open --species-list '$spList' for writing: $!\n";

while (<RAW>) {
    chomp;

    my ($filedate, $filetime, $filesize, $filename, $filesecs) = split '\s+', $_;

    my @fileComps   = split '/', $filename;
    my $speciesName = $fileComps[1];
    my $asmName     = $fileComps[3];
    my $seconds     = 0;

    if ($filename =~ m!species/\w+_\w+/metadata.yaml$!) {
        my @v = split '/', $_;
        $speciesMeta{$v[1]}++;
    }

    next if ($filename =~ m!/$!);          #  Why are you giving me directories?

    next if ($filename =~ m!^working!);
    next if ($filename =~ m!^galaxy!);
    next if ($filename =~ m!^species/insects!);

    next if ($filename =~ m!/intermediate!i);
    next if ($filename =~ m!/Intermidiates!i);           #  One guy has this.

    next if ($filename =~ m!.DS_Store$!);

    next if ($filename =~ m!.sh$!);

    next if ($filename =~ m!nohup.out$!);
    next if ($filename =~ m!nohup.err$!);

    next if ($filename =~ m!.log$!);
    next if ($filename =~ m!.log.gz$!);
    next if ($filename =~ m!.log.xz$!);

    next if ($filename =~ m!mother$!);
    next if ($filename =~ m!father$!);

    next if ($filename =~ m!.md5$!);
    next if ($filename =~ m!md5sum!);

    next if ($filename =~ m!readme!i);
    next if ($filename =~ m!readme.txt!i);
    next if ($filename =~ m!report.txt!i);
    next if ($filename =~ m!version.txt!i);

    next if ($filename =~ m!summary!i);
    next if ($filename =~ m!summary.txt!i);

    next if ($filename =~ m!manifest.txt!i);

    next if ($filename =~ m!upload_report!i);

    next if ($filename =~ m!.pdf$!);

    next if ($filename =~ m!.xml$!);

    next if ($filename =~ m!.csv$!); 
    next if ($filename =~ m!.csv.gz$!);
    next if ($filename =~ m!.csv.xz$!);

    next if ($filename =~ m!.tsv$!);
    next if ($filename =~ m!.tsv.gz$!);
    next if ($filename =~ m!.tsv.xz$!);

    next if ($filename =~ m!.zip$!);

    next if ($filename =~ m!.gff3$!);
    next if ($filename =~ m!.gff3.gz$!);
    next if ($filename =~ m!.gff3.xz$!);

    next if ($filename =~ m!.yaml$!);
    next if ($filename =~ m!.yml$!);

    next if ($filename =~ m!.tar$!);
    next if ($filename =~ m!.tar.gz$!);
    next if ($filename =~ m!.tar.xz$!);

    next if ($filename =~ m!/Stats!);

    next if ($filename =~ m!.gfastats$!);

    next if ($filename =~ m!/troubleshooting!i);
    next if ($filename =~ m!/annotation!i);
    next if ($filename =~ m!/test!i);
    next if ($filename =~ m!/analyses!i);
    next if ($filename =~ m!/comparative_analyses!i);

    next if ($filename =~ m!/transcriptomic_data!i);
    next if ($filename =~ m!/evaluation!i);
    next if ($filename =~ m!/bam_to_fasta!i);            #  fAngAng1/assembly_vgp_standard_1.6/bam_to_fasta (and others)
    next if ($filename =~ m!additional_haplotigs!);

    next if ($filename =~ m!aligned.genomecov!);

    next if ($filename =~ m!/qc/logs!i);
    next if ($filename =~ m!/qc/mash!i);
    next if ($filename =~ m!/qc/meryl!i);
    next if ($filename =~ m!/qc/stats!i);
    next if ($filename =~ m!/qc/busco!i);
    next if ($filename =~ m!/qc/asm-stats!i);
    next if ($filename =~ m!/qc/genomescope!i);
    next if ($filename =~ m!/genomescope!i);
    next if ($filename =~ m!/meryl_genomescope!i);
    next if ($filename =~ m!/katplot!i);
    next if ($filename =~ m!/FASTK!i);
    next if ($filename =~ m!/busco!i);
    next if ($filename =~ m!/merfin!i);
    next if ($filename =~ m!/mercury!i);

    next if ($filename =~ m!/assembly_v0/!);

    next if ($filename =~ m!bwgenome!i);

    next if ($filename =~ m!\.pretext$!i);
    next if ($filename =~ m!\.pretext\.gz$!i);
    next if ($filename =~ m!\.pretext\.png\.gz$!i);
    next if ($filename eq "species/Mesoplodon_densirostris/mMesDen1/assembly_curated/mMesDen1.pri.cur.20220330.png.gz");

    #  Some specific crud that is either large or useless.

    next if ($filename =~ m!/genomic_data/.*/scripts/!);

    next if ($filename =~ m!/genomic_data/nanopore/.*fast5$!);
    next if ($filename =~ m!/genomic_data/nanopore/.*ont_run_stats.txt$!);
    next if ($filename =~ m!/genomic_data/nanopore/.*telemetry.js.xz$!);

    next if ($filename =~ m!/aBomBom1/genomic_data/bionano/exp_refineFinal1/!);
    next if ($filename =~ m!/aBomBom1/genomic_data/pacbio/fasta!);
    next if ($filename =~ m!/bBucAby1/Test!);
    next if ($filename =~ m!/bCalAnn1/genomic_data/nanopore.*clean.fastq.gz$!);
    next if ($filename =~ m!/bGeoTri1/genomic_data/pacbio/old/!);
    next if ($filename =~ m!/fAloSap1/vgp_assembly_2.0/evaluation!);
    next if ($filename =~ m!/fAngAng1/assembly_vgp_standard_1.6/Scaffolding!);
    next if ($filename =~ m!/mCalJac1/SDA/!);
    next if ($filename =~ m!/mZalCal1/assembly_rockefeller_1.6/longranger!);
    next if ($filename =~ m!/rCheMyd1/assembly_vgp_standard_1.6/evaluation!);    #  LOTS of BUSCO intermediates
    next if ($filename =~ m!/sCarCar1/rawData/!);
    next if ($filename =~ m!/bHirRus1/Hirundo/!);
    next if ($filename =~ m!/mCalJac1/chrY/!);
    next if ($filename =~ m!/bSylAtr1/bSylAtr1_s4.arrow!);
    next if ($filename =~ m!/bTaeGut1/assembly_mt_milan/temp!);
    next if ($filename =~ m!/sCarCar1/assemblies/white_shark_09Feb2016_V8sgr!);

    #  Random single files (some included above)

    next if ($filename eq "species/Alosa_sapidissima/fAloSap1/vgp_assembly/fAloSap1.pri.v0.fasta.gz");
    next if ($filename eq "species/Bombina_bombina/aBomBom1/genomic_data/bionano/auto_noise/autoNoise1.errbin");
    next if ($filename eq "species/Bos_taurus/mBosTau1/assembly_adelaide/bionano_joyce/mBosTau1_mat_Saphyr_DLE1.bnx.gz");
    next if ($filename eq "species/Bos_taurus/mBosTau1/assembly_adelaide/bionano_joyce/mBosTau1_mat_Saphyr_DLE1_nonhap_ES_sdb_m3.cmap.gz");
    next if ($filename eq "species/Callithrix_jacchus/mCalJac1/assembly_nhgri_rockefeller_trio_1.6/pretext/mCalJac1.mat.bam");
    next if ($filename eq "species/Callithrix_jacchus/mCalJac1/assembly_nhgri_rockefeller_trio_1.6/pretext/mCalJac1.pat.bam");
    next if ($filename eq "species/Callithrix_jacchus/mCalJac1/qc/mCalJac1.k21.hist");
    next if ($filename eq "species/Calypte_anna/bCalAnn1/genomic_data/bionano/hybridScaffold_two_enzymes.xml.gz");
    next if ($filename eq "species/Carcharodon_carcharias/sCarCar1/assemblies/bams_chicago/chicagoLibrary1.final.scaffolds.snap.md.sorted.bam");
    next if ($filename eq "species/Carcharodon_carcharias/sCarCar1/assemblies/gapfilled_assembly/jelly.out.fasta");
    next if ($filename eq "species/Carcharodon_carcharias/sCarCar2/assembly_vgp_standard_1.6/sCarCar2_pri.asm.20200727.fasta.gz");
    next if ($filename eq "species/Cariama_cristata/bCarCri1/purge_haplotigs/curated.artefacts.fasta.gz");
    next if ($filename eq "species/Chiroxiphia_lanceolata/bChiLan1/assembly_ru_canu_2.1/bChiLan1.contigs.fasta.gz");
    next if ($filename eq "species/Chiroxiphia_lanceolata/bChiLan1/assembly_ru_canu_2.1/bChiLan1.correctedReads.fasta.gz");
    next if ($filename eq "species/Chiroxiphia_lanceolata/bChiLan1/assembly_ru_canu_2.1/bChiLan1.unassembled.fasta.gz");
    next if ($filename eq "species/Choloepus_didactylus/mChoDid1/assembly_berlinSanger_vgp_1.6/mChoDid1.mito.fa.gz");
    next if ($filename eq "species/Choloepus_didactylus/mChoDid1/assembly_berlin_vgp_1.5/VGP1.5.mChoDid1.berlin.txt");
    next if ($filename eq "species/Corvus_hawaiiensis/bCorHaw1/assembly_vgp_standard_2.0/Export");
    next if ($filename eq "species/Elephas_maximus/mEleMax1/assembly_curated/mEleMax1_1.pretext.savestate_May5_beta");
    next if ($filename eq "species/Falco_rusticolus/bFalRus1/assembly/bFalRus1.pri.asm.20200215.pretext");
    next if ($filename eq "species/Falco_rusticolus/bFalRus1/assembly/bFalRus1.pri.asm.20200401.pretext");
    next if ($filename eq "species/Heterohyrax_brucei/mHetBru1/haps_rapid_prtxt_XL.tpf");
    next if ($filename eq "species/Heterohyrax_brucei/mHetBru1/rapid_prtxt_XL_mHetBru1.agp");
    next if ($filename eq "species/Heterohyrax_brucei/mHetBru1/rapid_prtxt_XL_mHetBru1.tpf");
    next if ($filename eq "species/Melopsittacus_undulatus/bMelUnd1/pat_contigs_less20.fasta.gz");
    next if ($filename eq "species/Melopsittacus_undulatus/bMelUnd1/pat_contigs_over20.fasta.gz");
    next if ($filename eq "species/Mesoplodon_densirostris/mMesDen1/genomic_data/temp/forward_mMesDen1_reads.fastqsanger.gz");
    next if ($filename eq "species/Mesoplodon_densirostris/mMesDen1/genomic_data/temp/reverse_mMesDen1_reads.fastqsanger.gz");
    next if ($filename eq "species/Nyctibius_grandis/bNycGra1/bNycGra1_c2p2.fasta.gz");
    next if ($filename eq "species/Pan_troglodytes/mPanTro1/genomic_data/10x/path.list");
    next if ($filename eq "species/Pan_troglodytes/mPanTro1/genomic_data/pacbio/srx.list");
    next if ($filename eq "species/Pristis_pectinata/sPriPec2/assembly_vgp_standard_1.5/sPriPec2.pri.untrimmed.asm.20190802.fasta.gz");
    next if ($filename eq "species/Taeniopygia_guttata/bTaeGut1/bTaeGut1_s3q2.qv");
    next if ($filename eq "species/Taeniopygia_guttata/bTaeGut2/assembly_vgp_standard_2.0/Map");
    next if ($filename eq "species/Taeniopygia_guttata/bTaeGut2/assembly_vgp_standard_2.0/Map");
    next if ($filename eq "species/Taeniopygia_guttata/bTaeGut2/assembly_vgp_standard_2.0/Map");
    next if ($filename eq "species/Taeniopygia_guttata/bTaeGut2/assembly_vgp_standard_2.0/Map");
    next if ($filename eq "species/Taeniopygia_guttata/bTaeGut3/genomic_data/pat.txt");
    next if ($filename eq "species/Taeniopygia_guttata/bTaeGut4/genomic_data/mat.txt");
    next if ($filename eq "species/Zalophus_californianus/mZalCal1/assembly_rockefeller_1.6/mZalCal1_u1.fasta.gz");
    next if ($filename eq "species/Zalophus_californianus/mZalCal1/assembly_rockefeller_1.6/mZalCal1_u2.fasta.gz");

    #  Duplicate data
    #    fToxJac2 appears to have all hifi reads in 'ccsAll.bam' and filtered hifi reads in 'ccs.bam'.

    next if ($filename eq "species/Toxotes_jaculatrix/fToxJac2/genomic_data/pacbio_hifi/m54345U_200927_005508.ccsAll.bam");
    next if ($filename eq "species/Toxotes_jaculatrix/fToxJac2/genomic_data/pacbio_hifi/m54345U_200927_005508.ccsAll.bam.pbi");
    next if ($filename eq "species/Toxotes_jaculatrix/fToxJac2/genomic_data/pacbio_hifi/m64046_201211_121116.ccsAll.bam");
    next if ($filename eq "species/Toxotes_jaculatrix/fToxJac2/genomic_data/pacbio_hifi/m64046_201211_121116.ccsAll.bam.pbi");

    next if ($filename eq "species/Tauraco_erythrolophus/bTauEry1/genomic_data/arima/concatenated-bTauEry1_ARI16_4_USPD16090947_HV3M3CCXY_L3_1.clean.bam");
    next if ($filename eq "species/Tauraco_erythrolophus/bTauEry1/genomic_data/arima/concatenated-bTauEry1_ARI16_4_USPD16090947_HV3M3CCXY_L3_1.clean.bam.bai");

    #  These look like obsolete drafts.

    next if ($filename eq "species/Balaenoptera_musculus/mBalMus1/assembly_RU_repolished/mBalMus1.alt.cur.20200528.fasta");
    next if ($filename eq "species/Balaenoptera_musculus/mBalMus1/assembly_RU_repolished/mBalMus1.pri.cur.20200528.fasta");

    next if ($filename eq "species/Bos_taurus/mBosTau1/assembly_curated/mBosTau1.mat.alt.cur.20200213.fasta.gz");
    next if ($filename eq "species/Bos_taurus/mBosTau1/assembly_curated/mBosTau1.pat.alt.cur.20200213.fasta.gz");

    next if ($filename eq "species/Bos_taurus/mBosTau1/assembly_curated/purged/mBosTau1.mat.cur.20200213.purged.fasta.gz");
    next if ($filename eq "species/Bos_taurus/mBosTau1/assembly_curated/purged/mBosTau1.pat.cur.20200213.purged.fasta.gz");

    next if ($filename eq "species/Callithrix_jacchus/mCalJac1/assembly_nhgri_rockefeller_trio_1.6/VGP_mCalJac1_pat_X.fa.gz");

    next if ($filename eq "species/Corvus_hawaiiensis/bCorHaw1/assembly_vgp_standard_2.0/bCorHaw1.alt.asm.20212707.fasta");
    next if ($filename eq "species/Corvus_hawaiiensis/bCorHaw1/assembly_vgp_standard_2.0/bCorHaw1.pri.asm.20212907.fasta");

    next if ($filename eq "species/Falco_rusticolus/bFalRus1/assembly/bFalRus1.alt.asm.20200215.fasta");
    next if ($filename eq "species/Falco_rusticolus/bFalRus1/assembly/bFalRus1.alt.asm.20200401.fasta");
    next if ($filename eq "species/Falco_rusticolus/bFalRus1/assembly/bFalRus1.pri.asm.20200215.fasta");
    next if ($filename eq "species/Falco_rusticolus/bFalRus1/assembly/bFalRus1.pri.asm.20200401.fasta");

    next if ($filename eq "species/Pogoniulus_pusillus/bPogPus1/bPogPus1.alt.asm.20200207.fasta.gz");
    next if ($filename eq "species/Pogoniulus_pusillus/bPogPus1/bPogPus1.pri.asm.20200207.fasta.gz");

    next if ($filename eq "species/Scatophagus_argus/fScaArg1/assembly_vgp_standard_1.7/fScaArg1.alt.20210531.fa.gz");
    next if ($filename eq "species/Scatophagus_argus/fScaArg1/assembly_vgp_standard_1.7/fScaArg1.pri.20210531.fa.gz");

    next if ($filename eq "species/Sebastes_umbrosus/fSebUmb1/assembly_vgp_standard_1.6/fSebUmb1-t3.alt.asm.20200701.fasta.gz");
    next if ($filename eq "species/Sebastes_umbrosus/fSebUmb1/assembly_vgp_standard_1.6/fSebUmb1-t3.pri.asm.20200701.fasta.gz");

    next if ($filename eq "species/Suncus_etruscus/mSunEtr1/assembly_curated/mSunEtr1.alt.cur.20211125.nopipe.fasta.gz");

    next if ($filename eq "species/Taeniopygia_guttata/bTaeGut2/assembly_vgp_trio_2.0/bTaeGut2_trio.rebinned.hap1.s2.fasta.gz");
    next if ($filename eq "species/Taeniopygia_guttata/bTaeGut2/assembly_vgp_trio_2.0/bTaeGut2_trio.rebinned.hap2.s2.fasta.gz");

    next if ($filename eq "species/Tamandua_tetradactyla/mTamTet1/assembly_vgp_standard_1.7/mTamTet1_t3.alt.asm.20210813.fasta.gz");
    next if ($filename eq "species/Tamandua_tetradactyla/mTamTet1/assembly_vgp_standard_1.7/mTamTet1_t3.pri.asm.20210813.fasta.gz");

    print FILT "$_\n";   #  Finally, a file we want to process!

    my @v = split '/', $_;
    $speciesList{$v[1]}++;
}

close(RAW);
close(FILT);

#  Emit the list of speices.

foreach my $s (sort keys %speciesMeta) {
    delete $speciesList{$s};
    print SPLI "$s\n";
}
close(SPLI);

#  Emit a list of species without metadata.

if (scalar(keys %speciesList) > 0) {
    print "No metadata for:\n";

    foreach my $s (sort keys %speciesList) {
        print "  $s\n";
    }
}

exit(0);
