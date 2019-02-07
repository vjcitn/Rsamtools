# Migration of Rsamtools to Rhtslib


## What is the migration about?

Before this migration (i.e. for versions < 1.99.0), Rsamtools contained a
copy of an old (9-10 year old?) version of the samtools and tabix C code.

The latest version of samtools (as of Feb 6, 2019) is 1.9:
http://www.htslib.org/

In the recent years the samtools code has been split into samtools + htslib.
Most of the code that used to be in samtools is now in htslib. The latest
version of htslib is also 1.9.

Bioconductor package Rhtslib contains htslib 1.7, which is a decently recent
version of htslib.

After this migration (i.e. for versions >= 1.99.0), Rsamtools no longer
contains a copy of the samtools and tabix C code but compiles and links
against Rhtslib instead (i.e. against htslib 1.7).


## Current migration status

In Rsamtools 1.99.0 all the code from old Rsamtools was migrated to Rhtslib,
except `asBcf()` and `razip()`:

- `asBcf()` is currently disabled but will need to be migrated later if
  it turns out that some other package needs it (it doesn't seem to be
  the case though).

- `razip()` was simply removed. This is because recent samtools and
  htslib no longer support RAZF. From the samtools NEWS file:

    RAZF and razip are superseded by BGZF/bgzip and have been removed from
    samtools.

  This happened in samtools 1.0 (15 August, 2014).

Also, even though `applyPileups()` was migrated (`pileupbam.c`), 3 tests in
`inst/unitTests/test_applyPileups.R` are currently failing
(`test_applyPileups_byPosition`, `test_applyPileups_byPosition_yieldAll`,
and `test_applyPileups_byRange`).

`NEWS` file needs to be updated.


## Status of Bioconductor packages that depend directly on Rsamtools

As of Feb 6, 2019, 159 Bioconductor packages (154 software, 3 data-experiment,
and 2 workflow packages) depend **directly** on Rsamtools (i.e. via their
Depends, Imports, or LinkingTo field). So we've only checked manually a few
of them. This was on a 64-bit Ubuntu 16.04.5 LTS laptop running R 3.6
(2018-12-04 r75758) and Bioconductor 3.9 (current devel).

### All software packages with "LinkingTo: Rsamtools" were tested

  - VariantAnnotation: migrated (changes not committed yet),
    passes `R CMD check`

  - ArrayExpressHTS: migrated (changes not committed yet),
    passes `R CMD check`

  - BitSeq: migrated (changes not committed yet),
    passes `R CMD check`

  - DiffBind: migrated (changes not committed yet),
    passes `R CMD check`

  - h5vc: migrated (changes not committed yet, they also include
    fixing pre-existing build ERROR currently visible on the build report),
    passes `R CMD check`

  - podkat: migrated (changes not committed yet, include fixing
    `inst/examples/example1.vcf.gz` to make it compatible with bcftools 1.7),
    passes `R CMD check`

  - qrqc: migrated (changes not committed yet), passes `R CMD check`

  - QuasR: migrated (changes not committed yet), passes `R CMD check`

  - seqbias: migrated (changes not committed yet), passes `R CMD check`

  - TransView: migrated (changes not committed yet), passes `R CMD check`

### All software packages that use applyPileups() were tested

  - AllelicImbalance: contains old BCF file (`ERP000101.bcf`, in
    `extdata/ERP000101_subset/`) that needs to be fixed (it no longer
    works with recent bcftools).

  - biovizBase: no change needed, passes `R CMD check`

  - compEpiTools: no change needed, passes `R CMD check`

  - SICtools: `test_snpDiff()` (from `test_snpDiff.R`) FAILS!

  - VariantTools: no change needed, passes `R CMD check`

### A random sample of a few other software packages were tested

  - GenomicAlignments: no change needed, passes `R CMD check`

  - gmapR: does NOT compile (needs to be migrated)

### All workflow packages were tested

  - rnaseqGene: no change needed, vignette builds

  - sequencing: vignette does NOT build (because of AnnotationHub razip'ed
    FASTA file `AH18522`)


## Status of CRAN packages that depend directly on Rsamtools

As of Feb 6, 2019, 8 CRAN packages depend **directly** on Rsamtools (i.e.
via their Depends, Imports, or LinkingTo field).

  - BinQuasi: not tested yet

  - Brundle: not tested yet

  - ExomeDepth: not tested yet

  - hoardeR: not tested yet

  - NIPTeR: not tested yet

  - PlasmaMutationDetector: not tested yet

  - RAPIDR: not tested yet

  - spp: not tested yet


## What to do about razip-compressed FASTA files and their index files?

### The problem

RAZF is no longer supported in recent samtools, and old razip-compressed
FASTA files and their index files are now causing problems. Here is an example
from the sequencing workflow:
```
library(AnnotationHub)
ah <- AnnotationHub()
fa <- ah[["AH18522"]]
library(Rsamtools)
idx <- scanFaIndex(fa)  # still works
long <- idx[width(idx) > 82000]
getSeq(fa, param=long)  # ERROR! ('open' index failed)
```

### What to do about it?

The new way to go is to compress with bgzip and to recreate the index.

#### From the Unix command line

```
/path/to/samtools-1.7/htslib/bgzip myfile.fa       # creates myfile.fa,gz

/path/to/samtools-1.7/samtools faidx myfile.fa,gz  # generates index files
                                                   # myfile.fa.gz.fai
                                                   # and myfile.fa.gz.gzi
```

Note that the compression step is actually optional i.e. one can use
`samtools faidx` directly on the uncompressed FASTA file:

```
gunzip myfile.fa,gz  # uncompress to get myfile.fa back
/path/to/samtools-1.7/samtools faidx myfile.fa     # generates index file
                                                   # myfile.fa.fai
```

Note that in this case, only the `.fai` index file is generated
(no `.gzi` index file).

#### From R

```
library(Rsamtools)
bgzip("myfile.fa")
indexFa("myfile.fa.bgz")
```

Then create a FaFile object that can be used with `getSeq()` as usual:
```
fa <- FaFile("myfile.fa.bgz")
```

