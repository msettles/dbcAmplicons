/* Modifications made by Matt Settles and fall under the Copyright below
 * Copyright 2013, Institute for Bioninformatics and Evolutionary Studies
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *    http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Python.h"

#if defined(_MSC_VER)
typedef unsigned __int8 u_int8_t;
#endif

#define TRIM_VERSION    "0.1"

/* $Id: trim.c,v 0.1 2013/1/29 mls $ */

#ifndef MIN
# define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

typedef struct Tuple {
    int left_trim;
    int right_trim;
} Tuple;

/* 
Compute the trim points for a pair of reads, given min Quality
*/
struct Tuple
trim_sequence(const char *qual1, int q1len, const char *qual2, int q2len, int minQ)
{
    // qual1 and qual2 may be different lengths
    int i, j;

    Tuple val = { q1len , q2len }; // initialize to the read lengths

    if (q1len == 0 | q2len == 0)
        return (val);

    for (i = q1len-1; i > 0 ; i--) {
        if (((int)qual1[i] - 33) >= minQ ){
            val.left_trim = i+1; // add 1 as i is the first good base
            break;
        }
    }
    for (j = q2len-1; j > 0 ; j--) {
        if (((int)qual2[j] - 33) >= minQ ){
            val.right_trim = j+1; // add 1 as i is the first good base
            break;
        }
    }
    return (val);
}


PyDoc_STRVAR(trimseq_doc,
"trim(qual1, qual2, minQ) -> int, int\n\
    Determines the first postion from both qualities (from the end) at which the cooresponding quality value is >= minQ and returns those positions\n");

static PyObject *
trimseq(PyObject *self, PyObject *args)
{
    Tuple r;
    char *a, *b;
    int alen, blen, m;

    if (!PyArg_ParseTuple(args, "s#s#i", &a, &alen, &b, &blen, &m))
                return NULL;
    r = trim_sequence(a, alen, b, blen, m);
    if (r.left_trim== -1) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        return NULL;
    }
    return Py_BuildValue("{sisi}", "left_trim",r.left_trim, "right_trim",r.right_trim);
}

static PyMethodDef trim_methods[] = {
    {   "trim", (PyCFunction)trimseq,
        METH_VARARGS,   trimseq_doc       },
    {NULL, NULL, 0, NULL }  /* sentinel */
};

PyDoc_STRVAR(module_doc, "Trim sequence read pairs to minimum a read quality\n");

PyMODINIT_FUNC
inittrim(void)
{
    PyObject *m;

    m = Py_InitModule3("trim", trim_methods, module_doc);
    PyModule_AddStringConstant(m, "__version__", TRIM_VERSION);
}
