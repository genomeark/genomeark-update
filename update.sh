#!/bin/sh
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --time=2-0
#SBATCH --output=update.out

set -e

export AWS_PAGER=""  #  awscli v2 pipes all output through a pager, which we don't want.

updateFileList="yes"
updateGenBank="no"
updateMetadata="yes"
updateStats="yes"
updatePages="yes"
species=""

if [ -e "/work/software/ncbi-edirect/xtract" ] ; then
  export PATH=/work/software/ncbi-edirect:$PATH
else
  module load aws
  module load samtools
  module load edirect
fi

#  Pull in the submodules if they appear to not exist.
if [ ! -e "seqrequester/src/Makefile" ] ;  then git submodule update --init seqrequester;        fi
if [ ! -e "genomeark-metadata/species" ] ; then git submodule update --init genomeark-metadata;  fi
if [ ! -e "genomeark.github.io/404.md" ] ; then git submodule update --init genomeark.github.io; fi

#  Build seqrequester, if needed.
if [ ! -e "seqrequester/build/bin/seqrequester" ] ; then
  echo Building seqrequester.
  gmake -C seqrequester/src -j 4 > seqrequester.err 2>&1
fi


########################################
#
#  Parse the command line; turn on/off processing steps and create the list
#  of species to process.
#
while [ $# -gt 0 ] ; do
  opt=$1

  if   [ "$opt" = "--no-ls" ] ;         then updateFileList="no";                                                                                    shift
  elif [ "$opt" = "--no-genbank" ] ;    then                       updateGenBank="no";                                                               shift
  elif [ "$opt" = "--no-metadata" ] ;   then                                            updateMetadata="no";                                         shift
  elif [ "$opt" = "--no-stats" ] ;      then                                                                  updateStats="no";                      shift
  elif [ "$opt" = "--no-pages" ] ;      then                                                                                     updatePages="no";   shift

  elif [ "$opt" = "--only-ls" ] ;       then updateFileList="yes"; updateGenBank="no";  updateMetadata="no";  updateStats="no";  updatePages="no";   shift
  elif [ "$opt" = "--only-genbank" ] ;  then updateFileList="no";  updateGenBank="yes"; updateMetadata="no";  updateStats="no";  updatePages="no";   shift
  elif [ "$opt" = "--only-metadata" ] ; then updateFileList="no";  updateGenBank="no";  updateMetadata="yes"; updateStats="no";  updatePages="no";   shift
  elif [ "$opt" = "--only-stats" ] ;    then updateFileList="no";  updateGenBank="no";  updateMetadata="no";  updateStats="yes"; updatePages="no";   shift
  elif [ "$opt" = "--only-pages" ] ;    then updateFileList="no";  updateGenBank="no";  updateMetadata="no";  updateStats="no";  updatePages="yes";  shift

  elif [ "$opt" = "--ls" ] ;            then updateFileList="yes";                                                                                   shift
  elif [ "$opt" = "--genbank" ] ;       then                       updateGenBank="yes";                                                              shift
  elif [ "$opt" = "--metadata" ] ;      then                                            updateMetadata="yes";                                        shift
  elif [ "$opt" = "--stats" ] ;         then                                                                  updateStats="yes";                     shift
  elif [ "$opt" = "--pages" ] ;         then                                                                                     updatePages="yes";  shift

  else
    species="$species $opt"
    shift
  fi
done


mkdir -p downloads

########################################
#
#  Update genomeark metadata.
#
if [ $updateMetadata = "yes" ] ; then
  cd genomeark-metadata
  git fetch > ../metadata-fetch.out 2>&1
  git merge > ../metadata-merge.out 2>&1
fi


########################################
#
#  Download the genomeark file list.  This is obnoxiously slow.
#
if [ $updateFileList = "yes" ] ; then
  echo "Fetching list of files in s3://genomeark/species/."

  if [ -e downloads/genomeark.ls.raw ] ; then
    mv downloads/genomeark.ls.raw downloads/genomeark.ls.raw.$$
  fi

  aws --no-sign-request s3 ls --recursive s3://genomeark/species/ > downloads/genomeark.ls.raw.WORKING \
  && \
  mv downloads/genomeark.ls.raw.WORKING downloads/genomeark.ls.raw

  if [ ! -e "downloads/genomeark.ls.raw" ] ; then
    echo "Failed to fetch list of files in s3://genomeark/."
    exit 1
  fi
fi

#  ALWAYS filter the downloaded list if it is newer than the filtered
#  version.
#
if [ ! -e "downloads/genomeark.ls" -o \
          "downloads/genomeark.ls" -ot "downloads/genomeark.ls.raw" ] ; then
  echo "Filtering AWS filenames."
  perl scripts/filter-aws-bucket-list.pl --raw downloads/genomeark.ls.raw --filtered downloads/genomeark.ls --species-list downloads/species-list
fi



########################################
#
#  Update genbank accessions using Enterez Direct ('ncbi-edirect').
#  Versions 12.0 and above are known to work.
#
#    The 'xtract' command extract elements 'Genbank', 'AssemblyName' and
#    'AssemblyType' from all 'DocumentSummary' elements IF it contains a
#    `Synonym` element (which is where the 'Genbank' element is).
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

  for pp in PRJNA489243 PRJNA533106 ; do
    if [ -e "downloads/genbank.$pp.map.raw" ] ; then
      mv downloads/genbank.$pp.map.raw downloads/genbank.$pp.map.raw.$$
    fi

    echo "Fetching downloads/genbank.$pp.xml."

    esearch -db bioproject -query $pp \
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

    if [ ! -e "downloads/genbank.$pp.map" -o \
              "downloads/genbank.$pp.map" -ot "downloads/genbank.$pp.map.raw" ] ; then
      echo "Fixing downloads/genbank.$pp.xml."
      perl scripts/genbank-fixup.pl < downloads/genbank.$pp.map.raw > genbank/genbank.$pp.map
    fi
  done
fi



########################################
#
#  Fetch and summarize data, then
#  update markdown files in ../genomeark.github.io/
#

if [ "$updateStats" = "yes" ] ; then
  perl scripts/scan-bucket.pl $@ > scan-bucket.out 2> scan-bucket.err
fi

if [ "$updatePages" = "yes" ] ; then
  perl scripts/update-pages.pl $@ > update-pages.out 2> update-pages.err
fi

exit 0
