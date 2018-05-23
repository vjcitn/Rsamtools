#ifndef SCAN_BAM_DATA_H
#define SCAN_BAM_DATA_H

#include "khash.h"
#include "Rdefines.h"

#ifdef MIGRATE_ME

#include "bam_data.h"

#endif  /* MIGRATE_ME */

KHASH_SET_INIT_STR(str)

typedef struct {
    int *flag, *rname, *strand, *pos, *qwidth, *mapq, *mrnm, *mpos, *isize,
        *partition, *mates;
    char **qname;
    const char **cigar, **seq, **qual;
    khash_t(str) *cigarhash;
    int icnt, ncnt,
        mates_flag, partition_id; /* set prior to parsing 1 bam record */
    SEXP result;
} _SCAN_BAM_DATA, *SCAN_BAM_DATA;


SCAN_BAM_DATA _init_SCAN_BAM_DATA(SEXP result);
void _Free_SCAN_BAM_DATA(SCAN_BAM_DATA sbd);

#ifdef MIGRATE_ME

int _grow_SCAN_BAM_DATA(BAM_DATA bd, int len);
void _finish1range_SCAN_BAM_DATA(SCAN_BAM_DATA sbd, bam_header_t *header,
				 int irange);
SEXP _scan_bam_result_init(SEXP template_list, SEXP names, SEXP space,
                           BAM_FILE bfile);
SEXP _get_or_grow_SCAN_BAM_DATA(BAM_DATA bd, int len);

#endif  /* MIGRATE_ME */

#endif
