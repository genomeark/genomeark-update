#!/bin/sh

#  Use genbank-fetch.sh to pull the genbank XML and parse out a few fields.
#  Use genbank-fixup.pl to clean up the result and write genbank.map.

if [ -z $NCBI_API_KEY ] ; then
  echo NCBI_API_KEY not set.
  exit
fi

if [ $# -eq 0 ] ; then
  echo "usage: $0 [-ebp | -vgp]"
  exit 1
fi
for pp in $@ ; do
  if   [ $pp != "-ebp" -a $pp != "-vgp" ] ; then
    echo "usage: $0 [-ebp | -vgp]"
    exit 1
  fi
done



#  Versions:
#    edirect/8.60
#    edirect/10.0
#    edirect/14.5 - doesn't know the -q option.

module load edirect/10.0

for pp in $@ ; do
  if   [ $pp = "-ebp" ] ; then
    project="PRJNA533106"
  elif [ $pp = "-vgp" ] ; then
    project="PRJNA489243"
  else
    echo "usage: $0 [-ebp | -vgp]"
    exit 1
  fi

  if [ ! -e genbank.$project.xml ] ; then
    echo Fetching genbank.$project.xml.

    esearch -db bioproject -q $project \
    | \
    elink -db bioproject -target assembly -name bioproject_assembly_all \
    | \
    esummary \
    > genbank.$project.xml
  fi

  if [ ! -e genbank.$project.xml.map ] ; then
    echo Parsing genbank.$project.xml.

    #  Extract elements 'Genbank', 'AssemblyName' and 'AssemblyType' from all
    #  'DocumentSummary' elements IF it contains a `Synonym` element (which is
    #  where the 'Genbank' element is).

    xtract \
      -input genbank.$project.xml \
      -pattern DocumentSummary -if Synonym -element Genbank AssemblyName AssemblyType \
    | \
    sort -k2,2 \
    > genbank.$project.map.raw
  fi
done

exit 0
