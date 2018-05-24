#ifndef _TABIXFILE_H
#define _TABIXFILE_H

#include <Rdefines.h>
#include "htslib/tbx.h"

#ifdef MIGRATE_ME

// MIGRATION NOTE: The tabix_t struct is no longer defined in SAMtools/HTSlib.
// See src/tabix/tabix.h in the RELEASE_3_7 branch of Rsamtools for its
// definition. The problem is that it's not clear what we should replace it
// with. For example the tbx_t struct defined in htslib 1.7 (see htslib/tbx.h)
// looks VERY different from this tabix_t struct.
typedef struct {
    tabix_t *tabix;
    hts_itr_t iter;
} _TABIX_FILE;

#define TABIXFILE(b) ((_TABIX_FILE *) R_ExternalPtrAddr(b))

typedef SEXP SCAN_FUN(tabix_t *tabix, hts_itr_t iter, const int size,
                      SEXP state, SEXP rownames);

SCAN_FUN tabix_as_character;

SCAN_FUN tabix_count;

#endif  /* MIGRATE_ME */

SEXP tabixfile_init();
SEXP tabixfile_open(SEXP filename, SEXP indexname);
SEXP tabixfile_close(SEXP ext);
SEXP tabixfile_isopen(SEXP ext);

SEXP index_tabix(SEXP filename, SEXP format,
                 SEXP seq, SEXP begin, SEXP end,
                 SEXP skip, SEXP comment, SEXP zeroBased);
SEXP header_tabix(SEXP ext);
SEXP scan_tabix(SEXP ext, SEXP space, SEXP yield, SEXP fun, 
                SEXP state, SEXP rownames);

#endif
