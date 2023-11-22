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
#y %speciesMeta;   #  Defined if 'Genus_species/metadata.yaml' is found.
my %individuals;

my $lsRaw   = $ARGV[1];
my $lsFilt  = $ARGV[2];
my $spList  = undef;

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
if (($lsRaw  eq "") ||
    ($lsFilt eq "")) {
    print "usage: $0 --raw <input.ls> --filtered <output.ls> [--species-list <L>]\n";
    print "       $0 --raw downloads/genomeark.ls.raw --filtered downloads/genomeark.ls --species-list downloads/species-list\n";
    exit(1);
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


    if ($filename =~ m!species/(\w+_\w+)/([a-zA-Z]*)(\d)/!) {
        #print "$1/$2 -> $2$3 $filename\n";
        $individuals{"$1/$2"}{"$2$3"}++;
    }

    #  Subdirectories.

    next if ($filename =~ m!/$!);          #  Why are you giving me directories?

    next if ($filename =~ m!^working!);
    next if ($filename =~ m!^galaxy!);
    next if ($filename =~ m!^species/insects!);
    next if ($filename =~ m!^species/Mormyrids!);

    next if ($filename =~ m!/intermediate!i);
    next if ($filename =~ m!/troubleshooting!i);
    next if ($filename =~ m!/annotation!i);
    next if ($filename =~ m!/test!i);
    next if ($filename =~ m!/analyses!i);
    next if ($filename =~ m!/qc!i);
    next if ($filename =~ m!/comparative_analyses!i);
    next if ($filename =~ m!/transcriptomic_data!i);
    next if ($filename =~ m!/bam_to_fasta!i);
    next if ($filename =~ m!/other_data!i);              # t2t primates

    next if ($filename =~ m!/merfin!i);
    next if ($filename =~ m!/genomescope!i);
    next if ($filename =~ m!/Stats!i);

    #  Various log files.

    next if ($filename =~ m!\.DS_Store$!);

    next if ($filename =~ m!nohup\.(out|err)$!i);
    next if ($filename =~ m!log(|\.gz|\.xz)$!i);
    next if ($filename =~ m!out(|\.gz|\.xz)$!i);

    next if ($filename =~ m!transferdone$!i);

    next if ($filename =~ m!readme(|\.txt)$!i);
    next if ($filename =~ m!report(|\.txt)$!i);
    next if ($filename =~ m!summary(|\.txt)$!i);
    next if ($filename =~ m!manifest(|\.txt)$!i);

    next if ($filename =~ m!/genomic_data.*txt$!);

    next if ($filename =~ m!version(|\.txt)$!i);
    next if ($filename =~ m!md5$!i);
    next if ($filename =~ m!md5sum!i);
    next if ($filename =~ m!checksum!i);

    next if ($filename =~ m!upload_report!i);

    #  Scripts.

    next if ($filename =~ m!\.sh$!);
    next if ($filename =~ m!/scripts!);

    #  Descriptions.

    next if ($filename =~ m!\.html(|\.gz|\.xz)$!i);
    next if ($filename =~ m!\.json(|\.gz|\.xz)$!i);
    next if ($filename =~ m!\.md(|\.gz|\.xz)$!i);

    next if ($filename =~ m!\.pdf(|\.gz|\.xz)$!);
    next if ($filename =~ m!\.rtf(|\.gz|\.xz)$!);
    next if ($filename =~ m!\.xml(|\.gz|\.xz)$!);

    #  Images and analysee.

    next if ($filename =~ m!/evaluation!i);

    next if ($filename =~ m!\.png(|\.gz|\.xz)$!);
    next if ($filename =~ m!\.csv(|\.gz|\.xz)$!);
    next if ($filename =~ m!\.tsv(|\.gz|\.xz)$!);
    next if ($filename =~ m!\.gff3(|\.gz|\.xz)$!);

    next if ($filename =~ m!/assembly.*/.*/.*!);    #  ANYTHING in a subdirectory.

    next if ($filename =~ m!/FASTK!i);              #  species/Erethizon_dorsatum/mEreDor1/FASTK

    next if ($filename =~ m!gfastats$!);
    next if ($filename =~ m!bnstats(|\.yaml|\.yml)$!);
    next if ($filename =~ m!pb_?stats(|\.yaml|\.yml)$!);
    next if ($filename =~ m!hic_stats(|\.yaml|\.yml)$!);

    next if ($filename =~ m!/pretext!i);
    next if ($filename =~ m!pretext(|\.gz|\.xz)$!i);

    next if ($filename =~ m!/rapid_prtxt!i);

    #  Archives.

    next if ($filename =~ m!\.zip$!);
    next if ($filename =~ m!\.tar(|\.gz|\.xz)$!);

    #  Raw sequence/metadata and some Verkko generated files.

    next if ($filename =~ m!\.fai$!);
    next if ($filename =~ m!\.gzi$!);

    next if ($filename =~ m!fast5$!);
    next if ($filename =~ m!pod5$!);

    next if ($filename =~ m!/assembly.*/.*gb$!);   #  GenBank files in assembly directories.

    next if ($filename =~ m!(rdna|mito|ebv)-exemplar\.\d+\.fasta(|\.gz)$!);
    next if ($filename =~ m!unassigned\.\d+\.fasta(|\.gz)$!);
    next if ($filename =~ m!disconnected\.\d+\.fasta(|\.gz)$!);
    next if ($filename =~ m!analysis-\w+\.\d+\.fasta(|\.gz)$!);
    next if ($filename =~ m!\.gfa(|\.gz)$!);
    next if ($filename =~ m!\.gaf(|\.gz)$!);
    next if ($filename =~ m!\.bed(|\.gz)$!);

    #next if ($filename =~ m!mother$!);
    #next if ($filename =~ m!father$!);



    next if ($filename eq "species/Mesoplodon_densirostris/mMesDen1/assembly_curated/mMesDen1\.pri\.cur\.20220330\.png\.gz");

    #  Some specific crud that is either large or useless.

    next if ($filename =~ m!/aBomBom1/genomic_data/bionano/exp_refineFinal1/!);
    next if ($filename =~ m!/aBomBom1/genomic_data/pacbio/fasta!);
    next if ($filename =~ m!/bBucAby1/Test!);
    next if ($filename =~ m!/bCalAnn1/genomic_data/nanopore.*clean\.fastq\.gz$!);
    next if ($filename =~ m!/bGeoTri1/genomic_data/pacbio/old/!);
    next if ($filename =~ m!/fAloSap1/vgp_assembly_2\.0/evaluation!);
    next if ($filename =~ m!/fAngAng1/assembly_vgp_standard_1\.6/Scaffolding!);
    next if ($filename =~ m!/mCalJac1/SDA/!);
    next if ($filename =~ m!/mZalCal1/assembly_rockefeller_1\.6/longranger!);
    next if ($filename =~ m!/rCheMyd1/assembly_vgp_standard_1\.6/evaluation!);    #  LOTS of BUSCO intermediates
    next if ($filename =~ m!/sCarCar1/rawData/!);
    next if ($filename =~ m!/bHirRus1/Hirundo/!);
    next if ($filename =~ m!/mCalJac1/chrY/!);
    next if ($filename =~ m!/bSylAtr1/bSylAtr1_s4\.arrow!);
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
    next if ($filename eq "species/Cyclura_pinguis/rCycPin1/genomic_data/pacbio_hifi/slurm-42545496.out");
    next if ($filename eq "species/Apodemus_sylvaticus/mApoSyl1/genomic_data/pacbio_hifi/m64094e_220402_212648.ccs.rmdup.bam");
    next if ($filename eq "species/Apodemus_sylvaticus/mApoSyl1/genomic_data/pacbio_hifi/m64094e_220402_212648.ccs.rmdup.bam.pbi");
    next if ($filename eq "species/Pelobates_fuscus/aPelFus1/assembly_vgp_HiC_2.0/aPelFus1.HiC.pctg.20230308.fasta.gz");

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

    next if ($filename eq "species/Ara_ararauna/bAraAra1/assembly_vgp_HiC_2.0/bAraAra1.HiC.hap1.20220601.fasta.gz");   #  Misnamed HiC.hap1
    next if ($filename eq "species/Ara_ararauna/bAraAra1/assembly_vgp_HiC_2.0/bAraAra1.HiC.hap2.20220601.fasta.gz");   #  instead of hap1.asm

    next if ($filename eq "species/Eschrichtius_robustus/mEscRob2/assembly_vgp_standard_2.0/mEscRob2.standard.pri.20221005.fasta.gz");

    next if ($filename eq "species/Eubalaena_glacialis/mEubGla1/assembly_curated/mEubGla1.hap1.dup.20230118.fasta.gz");
    next if ($filename eq "species/Eubalaena_glacialis/mEubGla1/assembly_curated/mEubGla1.hap2.dup.20230118.fasta.gz");

    next if ($filename eq "species/Gypaetus_barbatus/bGypBar2/assembly_vgp_standard_2.0/bGypBar2.alt.asm.fasta.gz");           #  Missing date, but now
    next if ($filename eq "species/Gypaetus_barbatus/bGypBar2/assembly_vgp_standard_2.0/bGypBar2.pri.asm.20211119.fasta.gz");  #  obsolete draft

    next if ($filename eq "species/Gastrophryne_carolinensis/aGasCar1/assembly_vgp_HiC_2.0/aGasCar1.HiC.hap1.s2.fasta.gz");
    next if ($filename eq "species/Gastrophryne_carolinensis/aGasCar1/assembly_vgp_HiC_2.0/aGasCar1.HiC.hap2.s2.fasta.gz");

    next if ($filename eq "species/Malaclemys_terrapin/rMalTer1/assembly_curated/rMalTer1.alt.cur.20220708.fasta.gz");        #  Empty and no pri!
    next if ($filename eq "species/Malaclemys_terrapin/rMalTer1/assembly_curated/rMalTer1.hap2.alt.cur.20220708.fasta.gz");   #  Both empty and malformed

    next if ($filename eq "species/Gastrophryne_carolinensis/aGasCar1/assembly_curated/aGasCar1.hap1.alt.cur.20220713.fasta.gz");  #  Empty too.
    next if ($filename eq "species/Gastrophryne_carolinensis/aGasCar1/assembly_curated/aGasCar1.hap2.alt.cur.20220713.fasta.gz");  #  Empty.

    next if ($filename eq "species/Tamandua_tetradactyla/mTamTet1/assembly_curated/mTamTet1.pri.cur.20220203.noMT.fasta.gz'.");

    next if ($filename eq "species/Gorilla_gorilla/mGorGor1/genomic_data/pacbio_hifi/m54329U_210319_174352.hifi_reads.no-kinetics.bam");
    next if ($filename eq "species/Gorilla_gorilla/mGorGor1/genomic_data/pacbio_hifi/m64076_210326_192259.hifi_reads.no-kinetics.bam");
    next if ($filename eq "species/Pan_troglodytes/mPanTro3/genomic_data/pacbio_hifi/m64076_210810_005444.hifi_reads.no-kinetics.bam");
    next if ($filename eq "species/Pongo_abelii/mPonAbe1/genomic_data/pacbio_hifi/m54329U_210228_013345.hifi_reads.no-kinetics.bam");
    next if ($filename eq "species/Pongo_abelii/mPonAbe1/genomic_data/pacbio_hifi/m54329U_210404_020346.hifi_reads.no-kinetics.bam");
    next if ($filename eq "species/Pongo_abelii/mPonAbe1/genomic_data/pacbio_hifi/m64076_210225_020019.hifi_reads.no-kinetics.bam");
    next if ($filename eq "species/Pongo_abelii/mPonAbe1/genomic_data/pacbio_hifi/m64076_210330_204128.hifi_reads.no-kinetics.bam");



    next if ($filename eq "species/Gorilla_gorilla/mGorGor1/assembly_verkko_1.1-0.2-freeze/mGorGor1.analysis.20221111.fasta.gz");
    next if ($filename eq "species/Gorilla_gorilla/mGorGor1/assembly_verkko_1.1-0.2-freeze/mGorGor1.mat+MT.20221111.fasta.gz");

    next if ($filename eq "species/Pan_paniscus/mPanPan1/assembly_verkko_1.1-0.1-freeze/mPanPan1.analysis.20221111.fasta.gz");
    next if ($filename eq "species/Pan_paniscus/mPanPan1/assembly_verkko_1.1-0.1-freeze/mPanPan1.mat+MT.20221111.fasta.gz");

    next if ($filename eq "species/Pan_troglodytes/mPanTro3/assembly_verkko_1.1-hic-freeze/mPanTro3.analysis.20221111.fasta.gz");
    next if ($filename eq "species/Pan_troglodytes/mPanTro3/assembly_verkko_1.1-hic-freeze/mPanTro3.hap1+MT.20221111.fasta.gz");
    next if ($filename eq "species/Pan_troglodytes/mPanTro3/assembly_verkko_1.1-hic-freeze/mPanTro3.hap2+MT.20221111.fasta.gz");

    next if ($filename eq "species/Pongo_abelii/mPonAbe1/assembly_verkko_1.1-hic-freeze/mPonAbe1.analysis.20221111.fasta.gz");
    next if ($filename eq "species/Pongo_abelii/mPonAbe1/assembly_verkko_1.1-hic-freeze/mPonAbe1.hap1+MT.20221111.fasta.gz");
    next if ($filename eq "species/Pongo_abelii/mPonAbe1/assembly_verkko_1.1-hic-freeze/mPonAbe1.hap2+MT.20221111.fasta.gz");

    next if ($filename eq "species/Chanos_chanos/fChaCha1/genomic_data/pacbio_hifi/v2/ChanosChanos_HiFi_1cell_CCSreads_Q20.fastq.gz");
    next if ($filename eq "species/Chanos_chanos/fChaCha1/genomic_data/pacbio_hifi/v2/ChanosChanos_HiFi_PrimaryContigs_Polished.fasta.gz");
    next if ($filename eq "species/Chanos_chanos/fChaCha1/genomic_data/pacbio_hifi/v2/ChanosChanos_HiFi_PrimaryContigs_Polished.stats");

    next if ($filename eq "species/Pongo_pygmaeus/mPonPyg2/assembly_verkko_1.1-hic-freeze/mPonPyg2.analysis.20221111.fasta.gz");
    next if ($filename eq "species/Pongo_pygmaeus/mPonPyg2/assembly_verkko_1.1-hic-freeze/mPonPyg2.hap1+MT.20221111.fasta.gz");
    next if ($filename eq "species/Pongo_pygmaeus/mPonPyg2/assembly_verkko_1.1-hic-freeze/mPonPyg2.hap2+MT.20221111.fasta.gz");

    next if ($filename eq "species/Symphalangus_syndactylus/mSymSyn1/assembly_verkko_1.1-hic-freeze/mSymSyn1.analysis.20221111.fasta.gz");
    next if ($filename eq "species/Symphalangus_syndactylus/mSymSyn1/assembly_verkko_1.1-hic-freeze/mSymSyn1.hap1+MT.20221111.fasta.gz");
    next if ($filename eq "species/Symphalangus_syndactylus/mSymSyn1/assembly_verkko_1.1-hic-freeze/mSymSyn1.hap2+MT.20221111.fasta.gz");

    next if ($filename eq "species/Symphalangus_syndactylus/mSymSyn1/assembly_verkko_1.1-hic-freeze/mSymSyn1.analysis.20221111.fasta.gz");
    next if ($filename eq "species/Symphalangus_syndactylus/mSymSyn1/assembly_verkko_1.1-hic-freeze/mSymSyn1.hap1+MT.20221111.fasta.gz");
    next if ($filename eq "species/Symphalangus_syndactylus/mSymSyn1/assembly_verkko_1.1-hic-freeze/mSymSyn1.hap2+MT.20221111.fasta.gz");

    next if ($filename eq "species/Rissa_tridactyla/bRisTri1/assembly_vgp_trio_2.0/bRisTri1.trio.hap1.20220720.fasta.gz");
    next if ($filename eq "species/Rissa_tridactyla/bRisTri1/assembly_vgp_trio_2.0/bRisTri1.trio.hap2.20220720.fasta.gz");


    print FILT "$_\n";   #  Finally, a file we want to process!

    my @v = split '/', $_;
    $speciesList{$v[1]}++;
}

close(RAW);
close(FILT);

#  Emit the list of speices we found.  This is what drives scan-bucket and update-pages.

foreach my $s (sort keys %speciesList) {
    print SPLI "$s\n";
}
close(SPLI);

#  Create metadata templates if needed.

foreach my $ii (sort keys %individuals) {
    my ($species, $short_name, $name) = split '/', $ii;

    $name = $species;
    $name =~ s/_/ /g;

    #next   if (exists($speciesMeta{$species}));
    next   if (-e "genomeark-metadata/species/$species.yaml");
    next   if (-e "genomeark-metadata/species/$species.yaml.template");

    print " Create genomeark-metadata/species/$species.yaml.template ($short_name)\n";

    open(F, "> genomeark-metadata/species/$species.yaml.template") or die;

    print F "species:\n";
    print F "  name: $name\n";
    print F "  short_name: $short_name\n";
    print F "  common_name:\n";
    print F "  taxon_id:\n";
    print F "  order:\n";
    print F "    name:\n";
    print F "  family:\n";
    print F "    name:\n";
    print F "  individuals:\n";
    foreach my $indiv (sort keys %{$individuals{$ii}}) {
        print F "  - short_name: $indiv\n";
        print F "    biosample_id:\n";
        print F "    description:\n";
        print F "    provider:\n";
    }

    close(F);
}

exit(0);
