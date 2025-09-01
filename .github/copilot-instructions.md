# Rsamtools - R/Bioconductor Package

Rsamtools is an R/Bioconductor package that provides an interface to the `samtools`, `bcftools`, and `tabix` utilities for manipulating SAM (Sequence Alignment / Map), FASTA, binary variant call (BCF) and compressed indexed tab-delimited (tabix) files.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

### System Requirements and Dependencies
- Install system dependencies first:
  - `sudo apt-get update && sudo apt-get install -y r-base r-base-dev build-essential`
- Install required Bioconductor packages via apt:
  - `sudo apt-get install -y r-bioc-biocgenerics r-bioc-s4vectors r-bioc-iranges r-bioc-biostrings r-bioc-genomicranges r-bioc-rhtslib r-bioc-biocparallel`
  - `sudo apt-get install -y r-cran-knitr r-cran-rmarkdown r-cran-bitops`
- **CRITICAL DEPENDENCY**: This package requires Rhtslib >= 3.3.1, but Ubuntu packages may provide older versions (e.g., 2.4.1). Full installation requires proper Bioconductor repository setup.

### Building and Checking
- **Basic build** (FAST - 0.6 seconds): `R CMD build --no-build-vignettes .`
  - This creates `Rsamtools_X.Y.Z.tar.gz`
  - NEVER CANCEL: Even though this is fast, always set timeout to 60+ seconds for safety
- **Full build with vignettes** (requires complete dependency installation): `R CMD build .`
  - NEVER CANCEL: Set timeout to 60+ minutes as it requires package installation
- **Package check** (12-15 seconds): `R CMD check --no-vignettes --no-examples --no-tests --no-install Rsamtools_X.Y.Z.tar.gz`
  - NEVER CANCEL: Set timeout to 30+ minutes for complete checks

### Testing
- **Unit tests location**: `inst/unitTests/` contains extensive test suite
- **Test runner**: Execute via `tests/Rsamtools_unit_tests.R` which calls `BiocGenerics:::testPackage('Rsamtools')`
- **Manual test execution** requires full package installation with all dependencies
- **NEVER CANCEL**: Test suite can take 15+ minutes. Set timeout to 30+ minutes.

### Code Structure and Navigation
- **C/C++ source**: `src/` contains 68 C/C++ source files that compile against htslib via Rhtslib
- **R source**: `R/` contains 35+ R files with core functionality
- **Build configuration**: `src/Makevars` uses GNU make and links against Rhtslib
- **Key R files**:
  - `methods-BamFile.R` - BAM file handling
  - `methods-FaFile.R` - FASTA file handling  
  - `methods-TabixFile.R` - Tabix file handling
  - `scanBam.R` - Core BAM scanning functionality
  - `pileup.R` - Pileup operations

### Sample Data and Examples
- **Test data**: `inst/extdata/` contains sample BAM, FASTA, BCF, and VCF files:
  - `ex1.bam`, `ex1.bam.bai` - Example BAM file and index
  - `ce2dict1.fa`, `ce2dict1.fa.fai` - Example FASTA and index
  - `ex1.bcf.gz`, `ex1.vcf.gz` - Example variant files
- **Unit test cases**: `inst/unitTests/cases/` contains additional test files

### Validation Scenarios
After making changes, ALWAYS validate by:
1. **Build test**: `R CMD build --no-build-vignettes .` (should complete in ~1 second)
2. **Basic check**: `R CMD check --no-vignettes --no-examples --no-tests --no-install package.tar.gz` (should complete in ~12 seconds)
3. **Code compilation test**: Check that C/C++ code compiles without errors during build
4. **Documentation check**: Verify no new Rd warnings are introduced

### Common Commands and Timings
```bash
# Build package (0.6 seconds) - NEVER CANCEL, timeout 60+ seconds
R CMD build --no-build-vignettes .

# Check package (12 seconds) - NEVER CANCEL, timeout 30+ minutes  
R CMD check --no-vignettes --no-examples --no-tests --no-install Rsamtools_*.tar.gz

# Clean build artifacts
rm -f Rsamtools_*.tar.gz && rm -rf Rsamtools.Rcheck/
```

### Known Build Issues and Workarounds
- **Vignette building fails** without complete package installation due to dependency requirements
- **Package installation requires** all Bioconductor dependencies, particularly:
  - Seqinfo (may be part of GenomeInfoDb)  
  - Rhtslib >= 3.3.1
  - BiocParallel, GenomicRanges, Biostrings
- **Network connectivity issues** may prevent Bioconductor package installation via `install.packages()` or `BiocManager::install()`
- **Use apt packages** when available: `sudo apt-get install r-bioc-*` for Ubuntu/Debian systems

### Migration Notes
- Review `migration_notes.md` for historical context on migration from bundled samtools to Rhtslib
- Package previously contained embedded samtools/htslib code but now links against external Rhtslib
- Some functionality like `razip()` was removed as RAZF is deprecated in favor of BGZF/bgzip

### Development Workflow
1. Make minimal changes to source files
2. Build and check package: `R CMD build --no-build-vignettes . && R CMD check --no-vignettes --no-examples --no-tests --no-install Rsamtools_*.tar.gz`
3. Review build warnings and errors carefully
4. Test compilation of C/C++ code by examining build output
5. Clean up build artifacts before committing

### Common File Locations (for time-saving reference)

#### Repository Structure
```
.
├── DESCRIPTION          # Package metadata and dependencies
├── NAMESPACE           # Package exports and imports  
├── R/                  # R source code (35+ files)
├── src/                # C/C++ source (68 files) 
├── inst/
│   ├── extdata/        # Sample data files (BAM, FASTA, etc.)
│   ├── unitTests/      # Unit test suite
│   └── scripts/        # Utility scripts
├── tests/              # Test runner
├── vignettes/          # Package documentation
├── man/                # Generated documentation
└── migration_notes.md  # Historical migration information
```

#### Key Files to Check After Changes
- Always verify `src/Makevars` if modifying build configuration
- Check `NAMESPACE` if adding/removing exported functions
- Review `inst/unitTests/test_*.R` files if adding new functionality
- Update `NEWS` if making user-visible changes
