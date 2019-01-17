#include <stdlib.h>
#include <htslib/bgzf.h>
#include "tabixfile.h"
#include "utilities.h"

static SEXP TABIXFILE_TAG = NULL;

static const int TBX_INIT_SIZE = 32767;

static void _tabixfile_close(SEXP ext)
{
    _TABIX_FILE *tfile = TABIXFILE(ext);
    if (NULL != tfile->file) {
        hts_close(tfile->file);
        tfile->file = NULL;
    }
    if (NULL != tfile->index) {
        tbx_destroy(tfile->index);
        tfile->index = NULL;
    }
    if (NULL != tfile->iter) {
        tbx_itr_destroy(tfile->iter);
        tfile->iter = NULL;
    }
}

static void _tabixfile_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
        return;
    _tabixfile_close(ext);
    _TABIX_FILE *tfile = TABIXFILE(ext);
    Free(tfile);
    R_SetExternalPtrAddr(ext, NULL);
}

SEXP tabixfile_init()
{
    TABIXFILE_TAG = install("TabixFile");
    return R_NilValue;
}

SEXP tabixfile_open(SEXP filename, SEXP indexname)
{
    if (!IS_CHARACTER(filename) || 1L != Rf_length(filename))
        Rf_error("'filename' must be character(1)");
    if (!IS_CHARACTER(indexname) || 1L != Rf_length(indexname))
        Rf_error("'indexname' must be character(1)");

    _TABIX_FILE *tfile = Calloc(1, _TABIX_FILE);

    const char *fn = translateChar(STRING_ELT(filename, 0));
    tfile->file = hts_open(fn, "r");
    if (tfile->file == NULL) {
        Free(tfile);
        Rf_error("failed to open file: %s", fn);
    }

    const char *fnidx = translateChar(STRING_ELT(indexname, 0));
    tfile->index = tbx_index_load2(fn, fnidx);
    if (tfile->index == NULL) {
        hts_close(tfile->file);
        Free(tfile);
        Rf_error("failed to open index file: %s", fnidx);
    }

    tfile->iter = NULL;

    SEXP ext = PROTECT(R_MakeExternalPtr(tfile, TABIXFILE_TAG, filename));
    R_RegisterCFinalizerEx(ext, _tabixfile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}

SEXP tabixfile_close(SEXP ext)
{
    _checkext(ext, TABIXFILE_TAG, "close");
    _tabixfile_close(ext);
    return ext;
}

SEXP tabixfile_isopen(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != TABIXFILE(ext)) {
        _checkext(ext, TABIXFILE_TAG, "isOpen");
        if (TABIXFILE(ext)->file)
            ans = ScalarLogical(TRUE);
    }
    return ans;
}

SEXP index_tabix(SEXP filename, SEXP format, SEXP seq, SEXP begin, SEXP end,
                 SEXP skip, SEXP comment, SEXP zeroBased)
{
    tbx_conf_t conf = tbx_conf_gff;

    if (!IS_CHARACTER(filename) || 1L != Rf_length(filename))
        Rf_error("'filename' must be character(1)");

    const char *fn = translateChar(STRING_ELT(filename, 0));

    if (1L == Rf_length(format)) {
        const char *txt = CHAR(STRING_ELT(format, 0));
        if (strcmp(txt, "gff") == 0)
            conf = tbx_conf_gff;
        else if (strcmp(txt, "bed") == 0)
            conf = tbx_conf_bed;
        else if (strcmp(txt, "sam") == 0)
            conf = tbx_conf_sam;
        else if (strcmp(txt, "vcf") == 0 || strcmp(txt, "vcf4") == 0)
            conf = tbx_conf_vcf;
        else if (strcmp(txt, "psltbl") == 0)
            conf = tbx_conf_psltbl;
        else
            Rf_error("format '%s' unrecognized", txt);
    } else {
        if (!IS_INTEGER(seq) || 1L != Rf_length(seq))
            Rf_error("'seq' must be integer(1)");
        conf.sc = INTEGER(seq)[0];
        if (!IS_INTEGER(begin) || 1L != Rf_length(begin))
            Rf_error("'begin' must be integer(1)");
        conf.bc = INTEGER(begin)[0];
        if (!IS_INTEGER(end) || 1L != Rf_length(end))
            Rf_error("'end' must be integer(1)");
        conf.ec = INTEGER(end)[0];
    }

    if (IS_INTEGER(skip) && 1L == Rf_length(skip))
        conf.line_skip = INTEGER(skip)[0];
    if (IS_CHARACTER(comment) && 1L == Rf_length(comment))
        conf.meta_char = CHAR(STRING_ELT(comment, 0))[0];
    if (IS_LOGICAL(zeroBased) && 1L == Rf_length(zeroBased) &&
        TRUE == LOGICAL(zeroBased)[0])
        conf.preset |= TBX_UCSC;

    if (1 != bgzf_is_bgzf(fn))
        Rf_error("file does not appear to be bgzip'd");
    if (-1L == tbx_index_build(fn, 0, &conf))
        Rf_error("index build failed");

    return filename;
}

static const char *_tabix_read(htsFile *file, tbx_t *index, hts_itr_t *iter,
                               int *len)
{
    static kstring_t ksbuf = {0, 0, NULL};

    if (tbx_itr_next(file, index, iter, &ksbuf) < 0)
        return NULL;
    *len = ksbuf.l;
    return ksbuf.s;
}

SEXP _header_lines(htsFile *file, tbx_t *index, const tbx_conf_t * conf)
{
    const int GROW_BY = 100;
    SEXP lns;
    int i_lns = 0, pidx;

    /* Create iterator that will iterate over the entire file. */
    hts_itr_t *iter = tbx_itr_queryi(index, HTS_IDX_REST, 0, 0);
    if (NULL == iter)
        Rf_error("[internal] failed to obtain tabix iterator");

    uint64_t curr_off = 0;
    const char *s;
    int len;

    PROTECT_WITH_INDEX(lns = NEW_CHARACTER(0), &pidx);
    BGZF *fp = hts_get_bgzfp(file);
    if (NULL == fp)
        Rf_error("[internal] hts_get_bgzfp() failed");
    curr_off = bgzf_tell(fp);
    while (NULL != (s = _tabix_read(file, index, iter, &len))) {
        if ((int) (*s) != conf->meta_char)
            break;
	curr_off = bgzf_tell(fp);
        if (0 == (i_lns % GROW_BY)) {
            lns = Rf_lengthgets(lns, Rf_length(lns) + GROW_BY);
            REPROTECT(lns, pidx);
        }
        SET_STRING_ELT(lns, i_lns++, mkChar(s));
    }
    tbx_itr_destroy(iter);
    if (bgzf_seek(fp, curr_off, SEEK_SET) < 0)
        Rf_error("[internal] bgzf_seek() failed");

    lns = Rf_lengthgets(lns, i_lns);
    UNPROTECT(1);

    return lns;
}

SEXP header_tabix(SEXP ext)
{
    _checkext(ext, TABIXFILE_TAG, "headerTabix");
    _TABIX_FILE *tfile = TABIXFILE(ext);
    htsFile *file = tfile->file;
    tbx_t *index = tfile->index;

    SEXP result = PROTECT(NEW_LIST(5)), tmp, nms;
    nms = NEW_CHARACTER(Rf_length(result));
    Rf_namesgets(result, nms);
    SET_STRING_ELT(nms, 0, mkChar("seqnames"));
    SET_STRING_ELT(nms, 1, mkChar("indexColumns"));
    SET_STRING_ELT(nms, 2, mkChar("skip"));
    SET_STRING_ELT(nms, 3, mkChar("comment"));
    SET_STRING_ELT(nms, 4, mkChar("header"));

    /* seqnames */
    int n;
    const char **seqnames = tbx_seqnames(index, &n);
    if (n < 0)
        Rf_error("'seqnamesTabix' found <0 (!) seqnames");
    tmp = NEW_CHARACTER(n);
    SET_VECTOR_ELT(result, 0, tmp);
    for (int i = 0; i < n; ++i)
        SET_STRING_ELT(tmp, i, mkChar(seqnames[i]));
    free(seqnames);

    const tbx_conf_t conf = index->conf;

    /* indexColumns */
    tmp = NEW_INTEGER(3);
    SET_VECTOR_ELT(result, 1, tmp);
    INTEGER(tmp)[0] = conf.sc;
    INTEGER(tmp)[1] = conf.bc;
    INTEGER(tmp)[2] = conf.ec;
    nms = NEW_CHARACTER(3);
    Rf_namesgets(tmp, nms);
    SET_STRING_ELT(nms, 0, mkChar("seq"));
    SET_STRING_ELT(nms, 1, mkChar("start"));
    SET_STRING_ELT(nms, 2, mkChar("end"));

    /* skip */
    SET_VECTOR_ELT(result, 2, ScalarInteger(conf.line_skip));

    /* comment */
    char comment[2];
    comment[0] = (char) conf.meta_char;
    comment[1] = '\0';
    SET_VECTOR_ELT(result, 3, ScalarString(mkChar(comment)));

    /* header lines */
    SET_VECTOR_ELT(result, 4, _header_lines(file, index, &conf));

    UNPROTECT(1);
    return result;
}

SEXP tabix_as_character(htsFile *file, tbx_t *index, hts_itr_t *iter,
                        const int size, SEXP state, SEXP rownames)
{
    const double SCALE = 1.6;
    const tbx_conf_t conf = index->conf;

    int irec = 0, nrec, pidx;
    nrec = NA_INTEGER == size ? TBX_INIT_SIZE : size;
    SEXP record = NEW_CHARACTER(nrec);
    PROTECT_WITH_INDEX(record, &pidx);

    if (R_NilValue != rownames)
        Rf_error("[internal] expected 'NULL' rownames in tabix_as_character");
    if (R_NilValue != state)
        Rf_error("[internal] expected 'NULL' state in tabix_as_character");

    kstring_t ksbuf = {0, 0, NULL};

    while (tbx_itr_next(file, index, iter, &ksbuf) >= 0) {
        if (conf.meta_char == ksbuf.s[0])
            continue;

        if (irec == nrec) {
            nrec = nrec * SCALE;
            record = Rf_lengthgets(record, nrec);
            REPROTECT(record, pidx);
        }

        SET_STRING_ELT(record, irec, mkCharLen(ksbuf.s, ksbuf.l));

        irec += 1;

        if (NA_INTEGER != size && irec == nrec)
            break;
    }

    free(ksbuf.s);
    record = Rf_lengthgets(record, irec);
    UNPROTECT(1);
    return record;
}

SEXP tabix_count(htsFile *file, tbx_t *index, hts_itr_t *iter,
                 const int size, SEXP state, SEXP rownames)
{
    const tbx_conf_t conf = index->conf;
    const char *line;
    int linelen, irec = 0;

    if (R_NilValue != rownames)
        Rf_error("[internal] expected 'NULL' rownames in tabix_count");
    if (R_NilValue != state)
        Rf_error("[internal] expected 'NULL' state in tabix_count");

    while (NULL != (line = _tabix_read(file, index, iter, &linelen))) {
        if (conf.meta_char == *line)
            continue;
        irec += 1;
    }

    return ScalarInteger(irec);
}

SEXP scan_tabix(SEXP ext, SEXP regions, SEXP yield, SEXP fun,
                SEXP state, SEXP rownames)
{
    _checkparams(regions, R_NilValue, R_NilValue);
    if (!IS_INTEGER(yield) || 1L != Rf_length(yield))
        Rf_error("'yieldSize' must be integer(1)");

    _checkext(ext, TABIXFILE_TAG, "scanTabix");
    _TABIX_FILE *tfile = TABIXFILE(ext);
    htsFile *file = tfile->file;
    tbx_t *index = tfile->index;
    SCAN_FUN *scan = (SCAN_FUN *) R_ExternalPtrAddr(fun);

    SEXP space = VECTOR_ELT(regions, 0);
    const int nregions = Rf_length(space);
    SEXP result, elt;

    PROTECT(result = NEW_LIST(nregions == 0 ? 1 : nregions));
    if (0 == nregions) {
        hts_itr_t *iter = tfile->iter;
        if (NULL == iter) {
            /* Create iterator that will iterate over the entire file. */
            iter = tbx_itr_queryi(index, HTS_IDX_REST, 0, 0);
            if (NULL == iter)
                Rf_error("[internal] failed to obtain tabix iterator");
            tfile->iter = iter;
        }
        elt = scan(file, index, iter, INTEGER(yield)[0], state, rownames);
        SET_VECTOR_ELT(result, 0, elt);
    } else {
        const int
            *start = INTEGER(VECTOR_ELT(regions, 1)),
            *end = INTEGER(VECTOR_ELT(regions, 2));

        for (int i = 0; i < nregions; ++i) {
            int ibeg, iend, tid;
            hts_itr_t *iter;
            const char *tid_name;

            ibeg = start[i] == 0 ? 0 : start[i] - 1;
            iend = end[i];
            tid_name = CHAR(STRING_ELT(space, i));
            tid = tbx_name2id(index, tid_name);
            if (0 > tid)
                Rf_error("'%s' not present in tabix index", tid_name);
            iter = tbx_itr_queryi(index, tid, ibeg, iend);

            elt = scan(file, index, iter, NA_INTEGER, state, rownames);
            SET_VECTOR_ELT(result, i, elt);

            tbx_itr_destroy(iter);
        }
    }

    UNPROTECT(1);
    return result;
}

