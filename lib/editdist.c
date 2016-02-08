/* 
Contains 3 primary functions

### PRIMER ALIGNMENT
bounded_edit_distance - edit distance bounded to the 5' end with restrictions to the 3' end

edit_distance - standard edit distance

### BARCODE ALIGNMENT
hammingdist_distance - hamming distance 
*/

#include "Python.h"

#if defined(_MSC_VER)
typedef unsigned __int8 u_int8_t;
#endif

#define EDITDIST_VERSION    "0.6"

/* $Id: editdist.c,v 1.5 2007/05/03 23:36:36 djm Exp $ */
/* $Id: editdist.c,v 0.4 2013/12/31 mls $ */
/* $Id: editdist.c,v 0.4 2014/4/17 mls added hamming distnace */
/* $Id: editdist.c,v 0.6 2015/6/13 mls change to allow for faster processing */

#ifndef MIN3
# define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#endif

typedef struct Tuple {
    int dist; // The edit distance
    int spos; // Matching Start Position
    int epos; // Matching End Position
} Tuple;

/* 
compute the Levenstein distance between a (primer) and b (sequence read)
pegged to the 5' end and bounded by edit distance k
*/
struct Tuple
bounded_edit_distance(const char *primer, int primerlen, const char *seq, int seqlen, int f, int k, int m)
{
    // f is the pre-primer float value, k is max error and m is end matches
    unsigned int x, i, j, l, lastdiag, olddiag, cmin, endmatchcount;
    unsigned int column[primerlen - m + 1];
    Tuple val = { k+1, 0, primerlen}; // dist and positions
    /* a (primer) should always be < b (read) */
    if (primerlen > seqlen) { // primer should never be greater than the seq
        val.dist = -2;
        return (val);
    }

    if (primerlen == 0){
        val.dist = -2;
        return (val);
    }
    for (x = 0; x <= f; x++) {
        for (i = 1; i <= primerlen - m; i++)
            column[i] = i;
        for (i = 1; i <= primerlen - m + k ; i++) { // outer loop is the read
            column[0] = i;
            cmin = i;
            for (j = 1, lastdiag = i-1; j <= primerlen - m ; j++) { // inner loop is the primer
                olddiag = column[j];
                column[j] = MIN3(column[j] + 1, column[j-1] + 1, lastdiag + (primer[j-1] == seq[x+i-1] ? 0 : 1));
                lastdiag = olddiag;
                if (column[j] < cmin) cmin = column[j];
            }
            if (cmin > k) break; // if the smallest value in the column is > max error break
            if (column[primerlen - m] <= val.dist ){
                endmatchcount=0;
                for (l = 1; l <= m; l++){
                    if (primer[primerlen - l] != seq[x + i + m - l]){
                        break;
                    }
                    endmatchcount++;
                }
                if (endmatchcount == m){
                    val.dist = column[primerlen - m ]; // bottom right node, global alignment
                    val.spos = x;
                    val.epos = x + i + m;
                } 
            }
        }
    }
    return (val);
}


PyDoc_STRVAR(bounded_editdist_distance_doc,
"bounded_edit_distance(a, b, f, k, m) -> int, int\n\
    Calculates the bounded Levenshtein's edit distance between strings \"a\" and \"b\" with pre-string fload f, bound \"k\" and \"m\" matching bases at end\n");

static PyObject *
bounded_editdist_distance(PyObject *self, PyObject *args)
{
    Tuple r;
    char *primer, *seq;
    int primerlen, seqlen, f, k, m;

    if (!PyArg_ParseTuple(args, "s#s#ii", &primer, &primerlen, &seq, &seqlen, &f, &k, &m))
                return NULL;
    r = bounded_edit_distance(primer, primerlen, seq, seqlen, f, k, m);
    if (r.dist== -1) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        return NULL;
    }
    if (r.dist== -2) {
        PyErr_SetString(PyExc_SystemError, "Bad Arguments");
        return NULL;
    }
    return Py_BuildValue("iii", r.dist, r.spos, r.epos);
}

PyDoc_STRVAR(bounded_editdist_distance_list_doc,
"bounded_edit_distance(a, b, k, m) -> int, int\n\
    Calculates the bounded Levenshtein's edit distance between  a list of strings and a string with bound \"k\" and \"m\" matching bases at end\n");

static PyObject *
bounded_editdist_distance_list(PyObject *self, PyObject *args)
{
    Tuple r, s;
    Py_ssize_t i;
    char *seq;
    int seqlen, c = 0, f = 0, k = 0, m = 0; // c = current index of best, k = maxdist, m = finalmatch

    PyListObject* primer_list_o = NULL;

    if (!PyArg_ParseTuple(args, "O!s#iii", &PyList_Type, &primer_list_o, &seq, &seqlen, &f, &k, &m)){
        return NULL;
    }
    Py_ssize_t primer_list_o_length = PyList_GET_SIZE(primer_list_o);

    for (i = 0; i < primer_list_o_length; i++) {
        PyObject * primer_o = PyList_GET_ITEM(primer_list_o, i);
        if (PyString_Check(primer_o)) {
            r = bounded_edit_distance(PyString_AS_STRING(primer_o), (int)PyString_GET_SIZE(primer_o), seq, seqlen, f, k, m);
            if (r.dist== -1) {
                PyErr_SetString(PyExc_MemoryError, "Out of memory");
                return NULL;
            } else if (r.dist== -2) {
                PyErr_SetString(PyExc_SystemError, "Bad Arguments");
                return NULL;
            } else if ( i == 0){
                c = 0;
                s = r;
            } else if ( r.dist < s.dist ){
                c = (int)i;
                s = r;
            } else if ( s.dist == 0){ // can't get better than a perfect match!
                break;
            }
        } else {
            PyErr_SetString(PyExc_TypeError,
                "first argument must be a list of strings");
            return NULL;
        }
    }

    return Py_BuildValue("iiii", c, s.dist, s.spos, s.epos);
}

/*
Compute the Levenstein distance between a and b 
*/
static int
edit_distance(const char *a, int alen, const char *b, int blen)
{
    int tmplen, i, j;
    const char *tmp;
    int *current, *previous, *tmpl, add, del, chg, r;

    /* Swap to reduce worst-case memory requirement */
    if (alen > blen) {
        tmp = a;
        a = b;
        b = tmp;
        tmplen = alen;
        alen = blen;
        blen = tmplen;
    }

    if (alen == 0)
        return (blen);

    if ((previous = calloc(alen + 1, sizeof(*previous))) == NULL) {
        free(previous);
        return (-1);
    }
    if ((current = calloc(alen + 1, sizeof(*current))) == NULL) {
        free(current);
        return (-1);
    }

    for (i = 0; i < alen + 1; i++)
        previous[i] = i;

    for (i = 1; i < blen + 1; i++) {
        if (i > 1) {
            memset(previous, 0, (alen + 1) * sizeof(*previous));
            tmpl = previous;
            previous = current;
            current = tmpl;
        }
        current[0] = i;
        for (j = 1; j < alen + 1; j++) {
            add = previous[j] + 1;
            del = current[j - 1] + 1;
            chg = previous[j - 1];
            if (a[j - 1] != b[i - 1])
                chg++;
            current[j] = MIN3(add, del, chg);
        }
    }
    r = current[alen];
    free(previous);
    free(current);
    return (r);
}

PyDoc_STRVAR(editdist_distance_doc,
"distance(a, b) -> int\n\
    Calculates Levenshtein's edit distance between strings \"a\" and \"b\"\n");

static PyObject *
editdist_distance(PyObject *self, PyObject *args)
{
    char *a, *b;
    int alen, blen, r;

    if (!PyArg_ParseTuple(args, "s#s#", &a, &alen, &b, &blen))
                return NULL;
    r = edit_distance(a, alen, b, blen);
    if (r == -1) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        return NULL;
    }
    return PyInt_FromLong(r);
}

/*
Compute the hamming distance between a and b 
*/
static int
hammingdist_distance(const char *a, int alen, const char *b, int blen, int score)
{
    int i, diff;

    /* a and b should be same length, otherwise return error */
    if (alen != blen){
        return (-2);
    }
    /* length should not be 0, otherwise return error */
    if (alen == 0)
        return (-2);

    diff = 0;
    for (i = 0; i < alen; i++){
        if (a[i] != b[i])
            diff++;
        if (diff > score){
            return (diff);
        }
    }
    return (diff);
}

PyDoc_STRVAR(hamming_distance_doc,
"distance(a, b) -> int\n\
    Calculates hamming distance between two equal length strings \"a\" and \"b\"\n");

static PyObject *
hamming_distance(PyObject *self, PyObject *args)
{
    char *a, *b;
    int alen, blen, r;

    if (!PyArg_ParseTuple(args, "s#s#", &a, &alen, &b, &blen))
                return NULL;

    r = hammingdist_distance(a, alen, b, blen, blen);
    if (r == -2) {
        PyErr_SetString(PyExc_SystemError, "Bad Arguments");
        return NULL;
    }
    return PyInt_FromLong(r);
}

PyDoc_STRVAR(hamming_distance_list_doc,
"distance(a, b) -> int\n\
    Calculates hamming distance between a string and a list of strings\n");

static PyObject *
hamming_distance_list(PyObject *self, PyObject *args)
{
    char *b;
    Py_ssize_t i;
    int blen = 0, c = 0, s = 0, r = 0;

    PyListObject* barcode_list_o = NULL;

    if (!PyArg_ParseTuple(args, "O!s#i", &PyList_Type, &barcode_list_o, &b, &blen, &s)){
        return NULL;
    }
    Py_ssize_t barcode_list_o_length = PyList_GET_SIZE(barcode_list_o);
    for (i = 0; i < barcode_list_o_length; i++) {
        PyObject * barcode_o = PyList_GET_ITEM(barcode_list_o, i);
        if (PyString_Check(barcode_o)) {
            r = hammingdist_distance(PyString_AS_STRING(barcode_o), (int)PyString_GET_SIZE(barcode_o), b, blen, s);
            if (r == -2) {
                PyErr_SetString(PyExc_SystemError, "Bad Arguments");
                return NULL;
            } else if ( r < s ){
                c = (int)i;
                s = r;
            } else if ( s == 0 ){ // can't get better than a perfect match!
                break;
            }
        } else {
            PyErr_SetString(PyExc_TypeError,
                "first argument must be a list of strings");
            return NULL;
        }
    }
//    Py_DECREF(barcode_list_o);
    return Py_BuildValue("ii", c, s);
}

static PyMethodDef editdist_methods[] = {
    {   "distance", (PyCFunction)editdist_distance,
        METH_VARARGS,   editdist_distance_doc       },
    {   "bounded_distance", (PyCFunction)bounded_editdist_distance,
        METH_VARARGS,   bounded_editdist_distance_doc       },
    {   "bounded_distance_list", (PyCFunction)bounded_editdist_distance_list,
        METH_VARARGS,   bounded_editdist_distance_list_doc       },
    {   "hamming_distance", (PyCFunction)hamming_distance,
        METH_VARARGS,    hamming_distance_doc       },
    {   "hamming_distance_list", (PyCFunction)hamming_distance_list,
        METH_VARARGS,    hamming_distance_list_doc  },
    {   NULL, NULL, 0, NULL }  /* sentinel */
};

PyDoc_STRVAR(module_doc, "Calculate Hamming distance, Levenshtein's edit distance, and a edge bounded Levenshtein's edit distance.\n");

PyMODINIT_FUNC
initeditdist(void)
{
    PyObject *m;

    m = Py_InitModule3("editdist", editdist_methods, module_doc);
    PyModule_AddStringConstant(m, "__version__", EDITDIST_VERSION);
}
