This repository contains scripts and data for generating Jekyll pages in
https://github.com/genomeark/genomeark.github.io which are then displayed on
https://genomeark.github.io/.

Brief Instructions:
```
git clone git@github.com:genomeark/genomeark-update .
git submodule update --init --recursive
cd seqrequester/src
gmake
cd ../../
sh update.sh
```

`update.sh` needs `git`, `aws`, `samtools` and maybe NCBI edirect.

`update.sh` does:
  fetch any genomeark-metadata changes
  fetch the complete list of files in the bucket
  filter lots and lots of crud from the list
  fetch genbank data (disabled by default)
  call scan-bucket.pl
  call update-pages.pl

Most of the real work is done in scan-bucket.pl.  For a given species, it
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

All of these summaries are saved in github.

'update-pages.pl' takes those pre-computed summaries and massages them into
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
regenerate the summaries (unless I'm forgetting a cache somewhere)..  The
scripts run no git commands, so they're pretty safe to fiddle with - worst
case, you can just delete your local clone and start again.






