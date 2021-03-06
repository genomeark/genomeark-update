#!/bin/sh
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --time=24:00:00

updateFileList="no"
updateGenBank="yes"

module load aws
module load samtools

#  awscli v2 pipes all output through a pager, which we don't want.
export AWS_PAGER=""

mkdir -p downloads

#
#  Download the genomeark file list.  This is obnoxiously slow.
#  Once downloaded, filter out the crud if the filtered version is older.
#

if [ $updateFileList = "yes" -a ! -e "downloads/genomeark.ls.raw" ] ; then
  echo "Fetching list of files in s3://genomeark/species/."
  aws --no-sign-request s3 ls --recursive s3://genomeark/species/ > downloads/genomeark.ls.raw.WORKING \
  && \
  mv downloads/genomeark.ls.raw.WORKING downloads/genomeark.ls.raw

  if [ ! -e "downloads/genomeark.ls.raw" ] ; then
    echo "Failed to fetch list of files in s3://genomeark/."
    exit 1
  fi
fi

if [ ! -e "downloads/genomeark.ls" -o \
          "downloads/genomeark.ls" -ot "downloads/genomeark.ls.raw" ] ; then
  echo "Filtering filenames."
  perl scripts/filter-aws-bucket-list.pl < downloads/genomeark.ls.raw > downloads/genomeark.ls
fi

#
#  Update genbank accessions, if requested.
#    Versions:
#      edirect/8.60
#      edirect/10.0
#      edirect/14.5 - doesn't know the -q option.
#
#  Extract elements 'Genbank', 'AssemblyName' and 'AssemblyType' from all
#  'DocumentSummary' elements IF it contains a `Synonym` element (which is
#  where the 'Genbank' element is).
#
#  PRJNA533106 = ebp
#  PRJNA489243 = vgp
#
if [ $updateGenBank = "yes" ] ; then
  echo "Updating genbank mappings."

  if [ -z $NCBI_API_KEY ] ; then
    echo NCBI_API_KEY not set.
    exit
  fi

  if [ -e "/work/software/ncbi-edirect/xtract" ] ; then
    export PATH=/work/software/ncbi-edirect:$PATH
  fi
  if [ -e "/usr/local/apps/edirect/10.0/xtract" ] ; then
    module load edirect/10.0
  fi

  for pp in PRJNA489243 PRJNA533106 ; do
    if [ ! -e "downloads/genbank.$pp.map.raw" ] ; then
      echo "Fetching downloads/genbank.$pp.xml."

      esearch -db bioproject -q $pp \
      | \
      elink -db bioproject -target assembly -name bioproject_assembly_all \
      | \
      esummary \
      > downloads/genbank.$pp.xml

      echo "Parsing downloads/genbank.$pp.xml."

      xtract \
        -input downloads/genbank.$pp.xml \
        -pattern DocumentSummary -if Synonym -element Genbank AssemblyName AssemblyType \
      | \
      sort -k2,2 \
      > downloads/genbank.$pp.map.raw
    fi

    if [ ! -e "downloads/genbank.$pp.map" -o \
              "downloads/genbank.$pp.map" -ot "downloads/genbank.$pp.map.raw" ] ; then
      echo "Fixing downloads/genbank.$pp.xml."
      perl scripts/genbank-fixup.pl < downloads/genbank.$pp.map.raw > projects/genbank.$pp.map
    fi
  done
fi

#
#  Generate summaries of raw data and assemblies.
#  This runs in parallel across species.
#

#for x in `cat genomeark-metadata/species-list` ; do
#  echo $x
#done

#perl scripts/scan-bucket.pl $@
