This repository contains scripts and data for generating Jekyll pages in
https://github.com/genomeark/genomeark.github.io which are then displayed on
https://genomeark.org/.

## Brief Instructions

```
mkdir genomeark-update
cd genomeark-update
git clone git@github.com:genomeark/genomeark-update .
git submodule update --init --recursive
cd seqrequester/src
gmake
cd ../../
sbatch update.sh
```

This will make a copy of the primary https://github.com/genomeark/genomeark-update repository, which contains
the summary data for the genomeark bucket and scripts to create said data.  It
also includes, as submodules:
 - https://github.com/genomeark/genomeark.github.io - the Jekyll markdown for genomeark.org.
 - https://github.com/genomeark/genomeark-metadata - metadata for each and every species in genomeark.
 - https://github.com/genomeark/genomeark-analyses - a few extra analyses on some species.
 - https://github.com/marbl/seqrequester - a tool used for summarizing data.

To run the update, either `sh update.sh` or (much preferred!) submit to Slurm
with `sbatch update.sh`.  When lots of new genomic data is added to the
bucket, an update can take several hours.

`update.sh` does a few book-keeping tasks, notably, fetching the current list
of files in the bucket and filtering out lots and lots of crud, then calls
two perl scripts (scan-bucket.pl and update-pages.pl) to fetch and summarize
data, and the to convert those summaries into Jekyll markdown.

Once `update.sh` finishes, you need to manually commit any changes to the five
repositories listed above.

## Details

Most of the real work is done in `scan-bucket.pl`.  For a given species, it
attempts to summarize each file in the filtered file list.  A large
collection of rules is used to decide what type of data is in the file based
only on the file name and path.

For assemblies, it will download the whole fasta file and compute N50
statistics.

For read data, it will download a subset of the data (2 GB by default) and
compute the number of sequences and bases in the subset, then extrapolate to
the whole (compressed) size.  There is a second extrapolation -- that
probably should be removed -- that will download up to three files for each
data type then extrapolate to the whole set of files.

All of these summaries are saved in the genomeark-update repository, in directory
`species/`.  Downloaded data is saved in `downloads/species` and can be safely deleted.

'update-pages.pl' takes the pre-computed summaries and massages them into
mardown files for the web to display.  The complication in this process is
picking which assembly is the most recent one.

Once `update.sh` finishes, you're left with a bunch of changes to commit.

In genomeark-metadata, there might be new 'template' files.  These are created
when there is no existing metadata.  I've been committing them to the repo.  They will
be removed automagically if someone does add a real metadata file.

In genomeark.github.io, there will be lots of changes to the markdown files,
probably a few files renamed, and certainly new files created.  The lazy way
is to `git add _*/*.md` them.  I usually poke through the changes to make sure
nothing is messed up before committing.

In genomeark-update, there will be lots of temporary data in `downloads/`
that can be deleted.  The summary data is in `species/`; same idea, scan for
weird then add everything.

When genomeark.github.io is pushed, github will regenerate the pages for
genomeark.org.  You can track progress at
https://github.com/genomeark/genomeark.github.io/actions.

For testing the process, you can delete directory
`genomeark-update/species/Species_name/` then `update.sh Species_name` should
regenerate the summaries (unless I'm forgetting a cache somewhere).  The
scripts run no git commands, so they're pretty safe to fiddle with - worst
case, you can just delete your local clone and start again.






