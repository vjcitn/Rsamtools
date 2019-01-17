#include <string.h>
#include <errno.h>
#include <htslib/kstring.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <hts_internal.h>
#include "bcffile.h"
#include "utilities.h"

enum {
    BCF_HDR_REF = 0, BCF_HDR_SAMPLE, BCF_HDR_HEADER, BCF_HDR_LAST
};

static const char *BCF_HDR_NM[] = { "Reference", "Sample", "Header" };

enum {
    BCF_TID = 0, BCF_POS, BCF_ID, BCF_REF, BCF_ALT, BCF_QUAL,
    BCF_FLT, BCF_INFO, BCF_FMT, BCF_GENO, BCF_RECS_PER_RANGE,
    BCF_LAST
};

enum {
    BCF_TYPE_Integer = 1, BCF_TYPE_Float, BCF_TYPE_Character,
    BCF_TYPE_String, BCF_TYPE_Flag, BCF_TYPE_Last
};

static SEXP BCFFILE_TAG = NULL;

static const int BCF_BUFSIZE_GROW = 100000;	/* initial # records */

static htsFile *_bcf_tryopen(const char *fn, const char *mode)
{
    return vcf_open(fn, mode);
}

static const char *_find_index(const char *fn)
{
    static char fnidx2[999];

    const char *fnidx = hts_idx_getfn(fn, ".csi");
    if (fnidx == NULL) {
        fnidx = hts_idx_getfn(fn, ".tbi");
        if (fnidx == NULL)
            return NULL;
    }
    int size = snprintf(fnidx2, sizeof(fnidx2), "%s", fnidx);
    if (size >= sizeof(fnidx2))
        Rf_error("Rsamtools internal error in _find_index(): "
                 "fnidx2 string buffer too small");
    return fnidx2;
}

static void _bcf_close(htsFile * bcf, int errmsg)
{
    int err = vcf_close(bcf);
    if ((0 != err) && errmsg)
        Rf_error("_bcf_close error (%d)", err);
}

static void _bcffile_close(SEXP ext)
{
    _BCF_FILE *bfile = BCFFILE(ext);
    if (NULL != bfile->file)
        vcf_close(bfile->file);
    if (NULL != bfile->index)
        tbx_destroy(bfile->index);
    bfile->file = NULL;
    bfile->index = NULL;
}

static void _bcffile_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
        return;
    _bcffile_close(ext);
    _BCF_FILE *bfile = BCFFILE(ext);
    Free(bfile);
    R_SetExternalPtrAddr(ext, NULL);
}

SEXP bcffile_init()
{
    BCFFILE_TAG = install("BcfFile");
    return R_NilValue;
}

/*
SEXP bcffile_open(SEXP filename, SEXP indexname, SEXP filemode)
{
    _checknames(filename, indexname, filemode);

    _BCF_FILE *bfile = Calloc(1, _BCF_FILE);

    bfile->file = NULL;
    if (0 != Rf_length(filename)) {
        const char *fn = translateChar(STRING_ELT(filename, 0));
        bfile->file = _bcf_tryopen(fn, CHAR(STRING_ELT(filemode, 0)));
        if (NULL == bfile->file) {
            Free(bfile);
            Rf_error("'open' BCF failed\n  filename: %s", fn);
        }
    }

    bfile->index = NULL;
    if (0 != Rf_length(indexname) &&
        hts_get_format(bfile->file)->format != vcf)
    {
        const char *cindex = translateChar(STRING_ELT(indexname, 0));
        bfile->index = _bcf_idx_load(cindex);
        if (NULL == bfile->index) {
            _bcf_close(bfile->file, 0);
            Free(bfile);
            Rf_error("'open' BCF index failed\n  indexname: %s\n", cindex);
        }
    }

    SEXP ext = PROTECT(R_MakeExternalPtr(bfile, BCFFILE_TAG, filename));
    R_RegisterCFinalizerEx(ext, _bcffile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}
*/

SEXP bcffile_open(SEXP filename, SEXP indexname, SEXP filemode)
{
    _checknames(filename, indexname, filemode);
    if (Rf_length(filename) != 1)
        Rf_error("'filename' must have length 1");

    _BCF_FILE *bfile = Calloc(1, _BCF_FILE);

    const char *fn = translateChar(STRING_ELT(filename, 0));
    bfile->file = _bcf_tryopen(fn, CHAR(STRING_ELT(filemode, 0)));
    if (bfile->file == NULL) {
        Free(bfile);
        Rf_error("'open' VCF/BCF failed\n  filename: %s", fn);
    }
    bfile->index = NULL;
    // Rf_length(indexname) will be 0 when scanBcfHeader() is called on a
    // file path e.g. scanBcfHeader("chr22.vcf.gz")
    if (Rf_length(indexname) == 1) {
        const char *cindex = translateChar(STRING_ELT(indexname, 0));
        const char *fnidx = _find_index(cindex);
        if (fnidx == NULL) {
            _bcf_close(bfile->file, 0);
            Free(bfile);
            Rf_error("no VCF/BCF index found\n  filename: %s", fn);
        }
        bfile->index = tbx_index_load2(fn, fnidx);
        if (bfile->index == NULL) {
            _bcf_close(bfile->file, 0);
            Free(bfile);
            Rf_error("'open' VCF/BCF index failed\n  index file: %s\n", fnidx);
        }
    }

    SEXP ext = PROTECT(R_MakeExternalPtr(bfile, BCFFILE_TAG, filename));
    R_RegisterCFinalizerEx(ext, _bcffile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}

SEXP bcffile_close(SEXP ext)
{
    _checkext(ext, BCFFILE_TAG, "close");
    _bcffile_close(ext);
    return ext;
}

SEXP bcffile_isopen(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != BCFFILE(ext)) {
        _checkext(ext, BCFFILE_TAG, "isOpen");
        if (BCFFILE(ext)->file)
            ans = ScalarLogical(TRUE);
    }
    return ans;
}

/* BROKEN!! Always returns TRUE on an open BCF file! */
SEXP bcffile_isvcf(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != BCFFILE(ext)) {
        _checkext(ext, BCFFILE_TAG, "isVcf");
        if (BCFFILE(ext)->file &&
            hts_get_format(BCFFILE(ext)->file)->format == vcf)
        {
            ans = ScalarLogical(TRUE);
        }
    }
    return ans;
}

/* implementation */
static int _bcf_ans_grow(SEXP ans, R_len_t sz, int n_smpl)
{
    R_len_t n = sz;
    if (0 <= sz)
        n += Rf_length(VECTOR_ELT(ans, BCF_TID));
    else
        n *= -1;
    for (int i = 0; i < BCF_LAST; ++i) {
        SEXP elt = VECTOR_ELT(ans, i);
        switch (i) {
        case BCF_GENO:
            for (int i = 0; i < Rf_length(elt); ++i) {
                SEXP g = VECTOR_ELT(elt, i);
                SEXP dim = GET_DIM(g);
                if (R_NilValue == dim) {
                    g = Rf_lengthgets(g, n);
                    SET_VECTOR_ELT(elt, i, g);	/* protect */
                } else {
                    PROTECT(dim);
                    g = Rf_lengthgets(g, n_smpl * n);
                    SET_VECTOR_ELT(elt, i, g);	/* protect */
                    INTEGER(dim)[0] = n_smpl;
                    INTEGER(dim)[1] = n;
                    SET_DIM(g, dim);
                    UNPROTECT(1);
                }
            }
            break;
        case BCF_RECS_PER_RANGE:
            break;
        default:
            elt = Rf_lengthgets(elt, n);
            SET_VECTOR_ELT(ans, i, elt);	/* protect */
            break;
        }
    }
    return n;
}

SEXP scan_bcf_header(SEXP ext)
{
    _checkext(ext, BCFFILE_TAG, "scanBcfHeader");
    htsFile *bcf = BCFFILE(ext)->file;
    if (0 != bgzf_seek(bcf->fp.bgzf, 0, SEEK_SET)) {
        Rf_error("internal: failed to 'seek' on bcf file");
    }
    bcf_hdr_t *hdr = vcf_hdr_read(bcf);
    if (NULL == hdr)
        Rf_error("no 'header' line \"#CHROM POS ID...\"?");

    SEXP ans = PROTECT(NEW_LIST(BCF_HDR_LAST));
    SEXP ans_elt;
    int i;

    /* Get reference seqnames from header */
    int nseqs, seqname_len;
    const char **seqnames, *seqname;
    seqnames = bcf_hdr_seqnames(hdr, &nseqs);
    SET_VECTOR_ELT(ans, BCF_HDR_REF, NEW_STRING(nseqs));
    ans_elt = VECTOR_ELT(ans, BCF_HDR_REF);
    for (i = 0; i < nseqs; i++) {
        seqname = seqnames[i];
        seqname_len = _delete_trailing_LF_or_CRLF(seqname, -1);
        SET_STRING_ELT(ans_elt, i, mkCharLen(seqname, seqname_len));
    }
    free(seqnames);

    /* Get sample names from header */
    int nsamples, samplename_len;
    const char *samplename;
    nsamples = bcf_hdr_nsamples(hdr);
    SET_VECTOR_ELT(ans, BCF_HDR_SAMPLE, NEW_STRING(nsamples));
    ans_elt = VECTOR_ELT(ans, BCF_HDR_SAMPLE);
    for (i = 0; i < nsamples; i++) {
        samplename = hdr->samples[i];
        samplename_len = _delete_trailing_LF_or_CRLF(samplename, -1);
        SET_STRING_ELT(ans_elt, i, mkCharLen(samplename, samplename_len));
    }

    /* Get header lines from header */
    SET_VECTOR_ELT(ans, BCF_HDR_HEADER, NEW_STRING(hdr->nhrec));
    ans_elt = VECTOR_ELT(ans, BCF_HDR_HEADER);
    kstring_t str = {0, 0, NULL};
    for (i = 0; i < hdr->nhrec; i++) {
        str.l = 0;
        bcf_hrec_format(hdr->hrec[i], &str);
        str.l = _delete_trailing_LF_or_CRLF(str.s, str.l);
        SET_STRING_ELT(ans_elt, i, mkCharLen(str.s, str.l));
    }
    free(str.s);

    /* Set names on 'ans' */
    SEXP ans_names = NEW_CHARACTER(3);
    SET_NAMES(ans, ans_names);         /* protect */
    for (i = 0; i < BCF_HDR_LAST; ++i)
        SET_STRING_ELT(ans_names, i, mkChar(BCF_HDR_NM[i]));

    bcf_hdr_destroy(hdr);
    UNPROTECT(1);
    return ans;
}

/* See vcf_format() in Rhtslib/src/htslib-1.7/vcf.c for how to extract
   information from a VCF/BCF line represented by a bcf1_t structure.
   The _get_* functions below are following what vcf_format() does. */

static SEXP _get_CHROM(bcf1_t *bcf1, bcf_hdr_t *hdr)
{
    const char *key = hdr->id[BCF_DT_CTG][bcf1->rid].key;
    return mkChar(key);
}

static SEXP _get_ID(bcf1_t *bcf1)
{
    if (bcf1->d.id == NULL)
        return R_NaString;
    return mkChar(bcf1->d.id);
}

static SEXP _get_REF(bcf1_t *bcf1)
{
    if (bcf1->n_allele == 0)
        return R_NaString;
    return mkChar(bcf1->d.allele[0]);
}

static SEXP _get_ALT(bcf1_t *bcf1, kstring_t *ksbuf)
{
    ksbuf->l = 0;
    if (bcf1->n_allele > 1) {
        for (int i = 1; i < bcf1->n_allele; i++) {
            if (i > 1)
                kputc(',', ksbuf);
            kputs(bcf1->d.allele[i], ksbuf);
        }
    } else {
        kputc('.', ksbuf);
    }
    return mkCharLen(ksbuf->s, ksbuf->l);
}

static SEXP _get_FILTER(bcf1_t *bcf1, bcf_hdr_t *hdr, kstring_t *ksbuf)
{
    ksbuf->l = 0;
    if (bcf1->d.n_flt > 0) {
        for (int i = 0; i < bcf1->d.n_flt; i++) {
            if (i > 0)
                kputc(';', ksbuf);
            const char *key = hdr->id[BCF_DT_ID][bcf1->d.flt[i]].key;
            kputs(key, ksbuf);
        }
    } else {
        kputc('.', ksbuf);
    }
    return mkCharLen(ksbuf->s, ksbuf->l);
}

static SEXP _get_INFO(bcf1_t *bcf1, bcf_hdr_t *hdr, kstring_t *ksbuf)
{
    ksbuf->l = 0;
    int first = 1;
    if (bcf1->n_info != 0) {
        for (int i = 0; i < bcf1->n_info; i++) {
            bcf_info_t *z = &bcf1->d.info[i];
            if (!z->vptr)
                continue;
            if (first) {
                first = 0;
            } else {
                kputc(';', ksbuf);
            }
            if (z->key >= hdr->n[BCF_DT_ID]) {
                free(ksbuf->s);
                bcf_destroy(bcf1);
                error("Invalid VCF/BCF, the INFO index is too large");
            }
            const char *key = hdr->id[BCF_DT_ID][z->key].key;
            kputs(key, ksbuf);
            if (z->len <= 0)
                continue;
            kputc('=', ksbuf);
            if (z->len == 1) {
                int type = z->type;
                switch (type)
                {
                    case BCF_BT_INT8:
                        if (z->v1.i == bcf_int8_missing)
                            kputc('.', ksbuf);
                        else
                            kputw(z->v1.i, ksbuf);
                        break;
                    case BCF_BT_INT16:
                        if (z->v1.i == bcf_int16_missing)
                            kputc('.', ksbuf);
                        else
                            kputw(z->v1.i, ksbuf);
                        break;
                    case BCF_BT_INT32:
                        if (z->v1.i == bcf_int32_missing)
                            kputc('.', ksbuf);
                        else
                            kputw(z->v1.i, ksbuf);
                        break;
                    case BCF_BT_FLOAT:
                        if (bcf_float_is_missing(z->v1.f))
                            kputc('.', ksbuf);
                        else
                            kputd(z->v1.f, ksbuf);
                        break;
                    case BCF_BT_CHAR:
                        kputc(z->v1.i, ksbuf);
                        break;
                    default:
                        free(ksbuf->s);
                        bcf_destroy(bcf1);
                        error("Unexpected type %d", type);
                }
            } else {
                bcf_fmt_array(ksbuf, z->len, z->type, z->vptr);
            }
        }
    }
    if (first)
        kputc('.', ksbuf);
    return mkCharLen(ksbuf->s, ksbuf->l);
}

static SEXP _get_FORMAT(bcf1_t *bcf1, bcf_hdr_t *hdr, kstring_t *ksbuf)
{
    ksbuf->l = 0;
    int first = 1;
    if (bcf1->n_sample != 0) {
        int n_fmt = (int) bcf1->n_fmt;
        if (n_fmt != 0) {
            for (int i = 0; i < n_fmt; i++) {
                bcf_fmt_t *fmt = bcf1->d.fmt + i;
                if (fmt->p == NULL)
                    continue;
                if (first) {
                    first = 0;
                } else {
                    kputc(':', ksbuf);
                }
                int id = fmt->id;
                if (id < 0) {
                    free(ksbuf->s);
                    bcf_destroy(bcf1);
                    error("Invalid VCF/BCF: internal FORMAT id %d cannot "
                          "be matched to a FORMAT tag in the header", id);
                }
                const char *tag = hdr->id[BCF_DT_ID][id].key;
                kputs(tag, ksbuf);
            }
        }
    }
    if (first)
        kputc('.', ksbuf);
    return mkCharLen(ksbuf->s, ksbuf->l);
}

static void _load_int_vector(int *in, int n_sample,
                             const bcf_fmt_t *fmt, const char *tag)
{
    const void *data;
    int val;

    if (fmt->n != 1)
        error("FORMAT tag '%s' has an unexpected "
              "number of values per sample", tag);
    for (int j = 0; j < n_sample; j++) {
        data = fmt->p + j * fmt->size;
        switch (fmt->type) {
            case BCF_BT_INT8:
                val = ((uint8_t *) data)[0];
            break;
            case BCF_BT_INT16:
                val = ((uint16_t *) data)[0];
            break;
            default:
                error("FORMAT tag '%s' has an unexpected C type (%d)",
                      tag, fmt->type);
        }
        *(in++) = val;
    }
    return;
}

static SEXP _make_int_matrix(int nrow, int n_sample,
                             const bcf_fmt_t *fmt, const char *tag)
{
    const void *data;
    int val;
    SEXP ans;

    if (fmt->n != nrow)
        error("FORMAT tag '%s' has an unexpected "
              "number of values per sample", tag);
    ans = PROTECT(allocMatrix(INTSXP, nrow, n_sample));
    for (int j = 0; j < n_sample; j++) {
        data = fmt->p + j * fmt->size;
        for (int i = 0; i < nrow; i++) {
            switch (fmt->type) {
                case BCF_BT_INT8:
                    val = ((uint8_t *) data)[i];
                break;
                case BCF_BT_INT16:
                    val = ((uint16_t *) data)[i];
                break;
                default:
                    UNPROTECT(1);
                    error("FORMAT tag '%s' has an unexpected C type (%d)",
                          tag, fmt->type);
            }
            INTEGER(ans)[j * nrow + i] = val;
        }
    }
    UNPROTECT(1);
    return ans;
}

static SEXP _make_double_matrix(int nrow, int n_sample,
                                const bcf_fmt_t *fmt, const char *tag)
{
    const void *data;
    double val;
    SEXP ans;

    if (fmt->n != nrow)
        error("FORMAT tag '%s' has an unexpected "
              "number of values per sample", tag);
    if (fmt->type != BCF_BT_FLOAT)
        error("FORMAT tag '%s' has an unexpected C type (%d)", tag, fmt->type);
    ans = PROTECT(allocMatrix(REALSXP, nrow, n_sample));
    for (int j = 0; j < n_sample; j++) {
        data = fmt->p + j * fmt->size;
        for (int i = 0; i < nrow; i++) {
            val = ((float *) data)[i];
            REAL(ans)[j * nrow + i] = val;
        }
    }
    UNPROTECT(1);
    return ans;
}

static void _load_GENO(SEXP geno, const int n, bcf1_t *bcf1, bcf_hdr_t *hdr)
{
    /* FIXME: more flexible geno not supported by bcftools */
    int n_fmt = (int) bcf1->n_fmt;
    if (n_fmt == 0)
        return;
    int n_sample = bcf1->n_sample;
    SEXP geno_names = GET_NAMES(geno);
    for (int i = 0; i < n_fmt; i++) {
        const bcf_fmt_t *fmt = bcf1->d.fmt + i;
        if (fmt->p == NULL)
            continue;
        int id = fmt->id;
        if (id < 0) {
            bcf_destroy(bcf1);
            error("Invalid VCF/BCF: FORMAT tag id=%d not present "
                  "in header", id);
        }

        const char *tag = hdr->id[BCF_DT_ID][id].key;
        int t;
        for (t = 0; t < Rf_length(geno_names); t++) {
            const char *nm = CHAR(STRING_ELT(geno_names, t));
            if (strcmp(nm, tag) == 0)
                break;
        }
        if (Rf_length(geno_names) <= t)
            Rf_error("FORMAT tag '%s' not supported yet", tag);
        SEXP geno_elt = VECTOR_ELT(geno, t);

        int off = n * n_sample;  // could easily overflow! -> FIXME
        if (strcmp(tag, "PL") == 0) {
            /* 'geno_elt' is a list of matrices of integers */
            const int nrow = bcf1->n_allele * (bcf1->n_allele + 1) / 2;
            SEXP pl = _make_int_matrix(nrow, n_sample, fmt, tag);
            SET_VECTOR_ELT(geno_elt, n, pl);  // protect
        } else if (strcmp(tag, "DP") == 0 ||
                   strcmp(tag, "GQ") == 0 ||
                   strcmp(tag, "SP") == 0)
        {
            /* 'geno_elt' is an integer matrix with 1 row per sample */
            _load_int_vector(INTEGER(geno_elt) + off, n_sample, fmt, tag);
        } else if (strcmp(tag, "GT") == 0) {
            /* 'geno_elt' is a character matrix with 1 row per sample */
            if (fmt->type != BCF_BT_INT8)
                error("FORMAT tag '%s' has an unexpected C type (%d)",
                      tag, fmt->type);
            if (fmt->n != 1)
                error("FORMAT tag '%s' has an unexpected "
                      "number of values per sample", tag);
            char s[4];
            s[3] = '\0';
            for (int j = 0; j < n_sample; j++) {
                const void *data = fmt->p + j * fmt->size;
                int y = ((uint8_t *) data)[0];
                if (y >> 7 & 1)
                    SET_STRING_ELT(geno_elt, off++, mkChar("./."));
                else {
                    s[0] = '0' + (y >> 3 & 7);
                    s[1] = "/|"[y >> 6 & 1];
                    s[2] = '0' + (y & 7);
                    SET_STRING_ELT(geno_elt, off++, mkChar(s));
                }
            }
        } else if (strcmp(tag, "GL") == 0) {
            /* 'geno_elt' is a list of matrices of doubles */
            const int nrow = bcf1->n_allele * (bcf1->n_allele + 1) / 2;
            SEXP gl = _make_double_matrix(nrow, n_sample, fmt, tag);
            SET_VECTOR_ELT(geno_elt, n, gl);  // protect
        }
    }
}

static void _scan_bcf_line(bcf1_t *bcf1, bcf_hdr_t *hdr,
                           SEXP ans, int n, kstring_t *ksbuf)
{
    bcf_unpack(bcf1, BCF_UN_ALL);

    SEXP chrom = PROTECT(_get_CHROM(bcf1, hdr));
    SET_STRING_ELT(VECTOR_ELT(ans, BCF_TID), n, chrom);
    UNPROTECT(1);

    INTEGER(VECTOR_ELT(ans, BCF_POS))[n] = bcf1->pos + 1;

    SEXP id = PROTECT(_get_ID(bcf1));
    SET_STRING_ELT(VECTOR_ELT(ans, BCF_ID), n, id);
    UNPROTECT(1);

    SEXP ref = PROTECT(_get_REF(bcf1));
    SET_STRING_ELT(VECTOR_ELT(ans, BCF_REF), n, ref);
    UNPROTECT(1);

    SEXP alt = PROTECT(_get_ALT(bcf1, ksbuf));
    SET_STRING_ELT(VECTOR_ELT(ans, BCF_ALT), n, alt);
    UNPROTECT(1);

    REAL(VECTOR_ELT(ans, BCF_QUAL))[n] = bcf1->qual;

    SEXP filter = PROTECT(_get_FILTER(bcf1, hdr, ksbuf));
    SET_STRING_ELT(VECTOR_ELT(ans, BCF_FLT), n, filter);
    UNPROTECT(1);

    SEXP info = PROTECT(_get_INFO(bcf1, hdr, ksbuf));
    SET_STRING_ELT(VECTOR_ELT(ans, BCF_INFO), n, info);
    UNPROTECT(1);

    SEXP format = PROTECT(_get_FORMAT(bcf1, hdr, ksbuf));
    SET_STRING_ELT(VECTOR_ELT(ans, BCF_FMT), n, format);
    UNPROTECT(1);

    _load_GENO(VECTOR_ELT(ans, BCF_GENO), n, bcf1, hdr);
    return;
}

static int _scan_bcf_lines(htsFile *bcf, bcf_hdr_t *hdr, SEXP ans, int n)
{
    bcf1_t *bcf1 = bcf_init();  /* free'd in bcf_destroy */
    if (NULL == bcf1)
        Rf_error("_scan_bcf_lines: failed to allocate memory");
    int sz = Rf_length(VECTOR_ELT(ans, BCF_TID));
    kstring_t ksbuf = {0, 0, NULL};
    while (vcf_read(bcf, hdr, bcf1) >= 0) {
        if (n >= sz)
            sz = _bcf_ans_grow(ans, BCF_BUFSIZE_GROW, bcf_hdr_nsamples(hdr));
        if (n >= sz) {
            free(ksbuf.s);
            bcf_destroy(bcf1);
            Rf_error("_scan_bcf_lines: failed to increase size; out of memory?");
        }
        _scan_bcf_line(bcf1, hdr, ans, n, &ksbuf);
        ++n;
    }
    free(ksbuf.s);
    bcf_destroy(bcf1);
    return n;
}

static int _scan_bcf_region(htsFile *bcf, bcf_hdr_t *hdr, tbx_t *idx,
                            const char *spc, int start, int end,
                            SEXP ans, int n)
{
    int tid = tbx_name2id(idx, spc);
    if (tid == -1)
        Rf_error("'space' not in file: %s", spc);

    hts_itr_t *iter = tbx_itr_queryi(idx, tid, start - 1, end);
    if (iter == NULL)  /* invalid 'tid', should never happen */
        Rf_error("'space' not in file: %s", spc);

    bcf1_t *bcf1 = bcf_init();  /* free'd in bcf_destroy */
    if (NULL == bcf1) {
        tbx_itr_destroy(iter);
        Rf_error("_scan_bcf_region: failed to allocate memory");
    }
    int sz = Rf_length(VECTOR_ELT(ans, BCF_TID));
    kstring_t ksbuf = {0, 0, NULL};
    while (tbx_itr_next(bcf, idx, iter, &ksbuf) >= 0) {
        if (vcf_parse1(&ksbuf, hdr, bcf1) < 0) {
            free(ksbuf.s);
            bcf_destroy(bcf1);
            tbx_itr_destroy(iter);
            Rf_error("_scan_bcf_region: parse error");
        }
        if (n >= sz)
            sz = _bcf_ans_grow(ans, BCF_BUFSIZE_GROW, bcf_hdr_nsamples(hdr));
        if (n >= sz) {
            free(ksbuf.s);
            bcf_destroy(bcf1);
            tbx_itr_destroy(iter);
            Rf_error("_scan_bcf_region: failed to increase size; out of memory?");
        }
        _scan_bcf_line(bcf1, hdr, ans, n, &ksbuf);
        ++n;
    }
    free(ksbuf.s);
    bcf_destroy(bcf1);
    tbx_itr_destroy(iter);
    return n;
}

SEXP scan_bcf(SEXP ext, SEXP regions, SEXP tmpl)
{
    _checkparams(regions, R_NilValue, R_NilValue);
    _checkext(ext, BCFFILE_TAG, "scanBcf");
    htsFile *bcf = BCFFILE(ext)->file;
    if (0 != bgzf_seek(bcf->fp.bgzf, 0, SEEK_SET)) {
        Rf_error("internal: failed to 'seek' on bcf file");
    }
    bcf_hdr_t *hdr = vcf_hdr_read(bcf);
    if (NULL == hdr)
        Rf_error("no 'header' line \"#CHROM POS ID...\"?");

    int n = 0;
    SEXP ans = PROTECT(Rf_duplicate(tmpl));

    if (regions == R_NilValue) {
        SET_VECTOR_ELT(ans, BCF_RECS_PER_RANGE, NEW_INTEGER(1));
        n = _scan_bcf_lines(bcf, hdr, ans, n);
        INTEGER(VECTOR_ELT(ans, BCF_RECS_PER_RANGE))[0] = n;
    } else {
        tbx_t *idx = BCFFILE(ext)->index;
        SEXP space = VECTOR_ELT(regions, 0);
        const int *start = INTEGER(VECTOR_ELT(regions, 1)),
                  *end = INTEGER(VECTOR_ELT(regions, 2)),
                  nregions = Rf_length(space);
        SEXP nrec = NEW_INTEGER(nregions);
        SET_VECTOR_ELT(ans, BCF_RECS_PER_RANGE, nrec);
        for (int i = 0; i < nregions; ++i) {
            const char *spc = CHAR(STRING_ELT(space, i));
            n = _scan_bcf_region(bcf, hdr, idx, spc, start[i], end[i], ans, n);
            if (i == 0)
                INTEGER(nrec)[i] = n;
            else
                INTEGER(nrec)[i] = n - INTEGER(nrec)[i - 1];
        }
    }
    _bcf_ans_grow(ans, -1 * n, bcf_hdr_nsamples(hdr));

    UNPROTECT(1);
    return ans;
}

#ifdef MIGRATE_ME

static int _as_bcf(htsFile * fin, const char *dict, htsFile * fout)
{
    bcf1_t *bcf1 = bcf_init();  /* free'd in bcf_destroy */
    if (NULL == bcf1)
        Rf_error("_as_bcf: failed to allocate memory");
    bcf_hdr_t *hin, *hout;
    int r, count = 0;

    hin = hout = vcf_hdr_read(fin);
    vcf_dictread(fin, hin, dict);
    vcf_hdr_write(fout, hout);
    while (0 <= (r = vcf_read(fin, hin, bcf1))) {
        if (NULL == bcf1->ref)
            Rf_error("cannot (yet) coerce VCF files without FORMAT");
        vcf_write(fout, hout, bcf1);
        count++;
    }

    if (hin != hout)
        bcf_hdr_destroy(hout);
    bcf_hdr_destroy(hin);
    bcf_destroy(bcf1);

    return r >= -1 ? count : -1 * count;
}

SEXP as_bcf(SEXP file, SEXP dictionary, SEXP destination)
{
    if (!IS_CHARACTER(file) || 1 != LENGTH(file))
        Rf_error("'file' must be character(1)");
    if (!IS_CHARACTER(dictionary) || 1 != LENGTH(dictionary))
        Rf_error("'dictionary' must be character(1)");
    if (!IS_CHARACTER(destination) || 1 != LENGTH(destination))
        Rf_error("'destination' must be character(1)");

    htsFile *fin = _bcf_tryopen(translateChar(STRING_ELT(file, 0)), "r");
    if (NULL == fin)
        Rf_error("failed to open VCF 'file'");

    htsFile *fout = _bcf_tryopen(translateChar(STRING_ELT(destination, 0)), "wb");
    if (NULL == fout)
        Rf_error("failed to open BCF 'destination'");

    int count = _as_bcf(fin, translateChar(STRING_ELT(dictionary, 0)), fout);

    _bcf_close(fin, 0);
    _bcf_close(fout, 0);
    if (count < 0)
        Rf_error("truncated input file at record %d", -1 * count);

    return destination;
}

#endif  /* MIGRATE_ME */

SEXP index_bcf(SEXP file)
{
    if (!IS_CHARACTER(file) || 1 != LENGTH(file))
        Rf_error("'file' must be character(1)");
    const char *fbcf = translateChar(STRING_ELT(file, 0));
    int status = bcf_index_build(fbcf, 0);
    if (0 != status)
        Rf_error("failed to build index");
    char *fidx = (char *) R_alloc(strlen(fbcf) + 5, sizeof(char));
    sprintf(fidx, "%s.bci", fbcf);
    return mkString(fidx);
}
