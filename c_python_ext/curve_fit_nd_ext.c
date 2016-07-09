
#include <Python.h>

#include "curve_fit_nd.h"

PyDoc_STRVAR(M_Curve_fit_nd_doc,
"TODO\n"
);


PyDoc_STRVAR(M_Curve_fit_nd_curve_from_points_doc,
".. function:: curve_from_points(points, error)\n"
"\n"
"   Returns the newly calculated curve.\n"
"\n"
"   :arg line: Points representing a line\n"
"   :type line: list\n"
"   :arg error: Error threshold.\n"
"   :type error: float\n"
"   :return: The point of intersection or None if no intersection is found\n"
"   :rtype: list of float tuples\n"
);
static PyObject *M_Curve_fit_nd_curve_from_points(PyObject *self, PyObject *args)
{
	(void)self;

	const char *error_prefix = "curve_from_points";
	PyObject *points;
	PyObject *points_fast;
	double error_threshold;


	if (!PyArg_ParseTuple(
	        args, "Od:curve_from_points",
	        &points,
	        &error_threshold) ||
	    !(points_fast = PySequence_Fast(points, error_prefix)))
	{
		return NULL;
	}

	const unsigned int points_len = PySequence_Fast_GET_SIZE(points_fast);
	if (points_len == 0) {
		Py_DECREF(points_fast);
		return PyList_New(0);
	}

	PyObject **points_array = PySequence_Fast_ITEMS(points_fast);
	double *points_data = NULL;
	unsigned int dims = 0;

	for (unsigned int i = 0; i < points_len; i++) {
		PyObject *item = points_array[i];
		PyObject *item_fast = PySequence_Fast(item, "curve_from_points item");
		if (item_fast == NULL) {
			if (points_data != NULL) {
				PyMem_Free(points_data);
			}
			Py_DECREF(points_fast);
			return NULL;
		}

		{
			unsigned int item_dims = PySequence_Fast_GET_SIZE(item_fast);
			if (i == 0) {
				if (item_dims == 0) {
					PyErr_SetString(PyExc_ValueError, "empty item");
					Py_DECREF(points_fast);
					Py_DECREF(item_fast);
					return NULL;
				}
				else {
					dims = item_dims;
					points_data = PyMem_Malloc((size_t)points_len * dims * sizeof(double));
				}
			}

			if (item_dims != dims) {
				PyErr_SetString(PyExc_ValueError, "item size mismatch");
				Py_DECREF(points_fast);
				Py_DECREF(item_fast);
				PyMem_Free(points_data);
				return NULL;
			}
		}

		PyObject **item_array = PySequence_Fast_ITEMS(item_fast);
		for (unsigned int j = 0; j < dims; j++) {
			const double number = PyFloat_AsDouble(item_array[j]);
			if ((number == -1.0) && PyErr_Occurred()) {
				Py_DECREF(points_fast);
				Py_DECREF(item_fast);
				PyMem_Free(points_data);
				return NULL;
			}
			points_data[(i * dims) + j] = number;
		}
		Py_DECREF(item_fast);
	}

	Py_DECREF(points_fast);

	double *cubic_array = NULL;
	unsigned int cubic_array_len = 0;

	if (curve_fit_cubic_to_points_db(
	        points_data, points_len, dims, error_threshold, 0,
	        NULL, 0,

	        &cubic_array, &cubic_array_len,
	        NULL, NULL, NULL) != 0)
	{

		PyErr_SetString(PyExc_ValueError, "error fitting the curve");
		PyMem_Free(points_data);
		return NULL;
	}

	PyMem_Free(points_data);

	PyObject *ret = PyList_New(cubic_array_len);
	double *c = cubic_array;

	for (unsigned int i = 0; i < cubic_array_len; i++) {
		PyObject *item = PyTuple_New(3);
		for (unsigned int h = 0; h < 3; h++) {
			PyObject *item_point = PyTuple_New(dims);
			for (unsigned int j = 0; j < dims; j++) {
				PyTuple_SET_ITEM(item_point, j, PyFloat_FromDouble(*c));
				c++;
			}
			PyTuple_SET_ITEM(item, h, item_point);
		}
		PyList_SET_ITEM(ret, i, item);
	}
	assert(c == (cubic_array + (cubic_array_len * dims * 3)));

	return ret;
}

static struct PyMethodDef M_Curve_fit_nd_methods[] = {
	{"curve_from_points", (PyCFunction) M_Curve_fit_nd_curve_from_points, METH_VARARGS, M_Curve_fit_nd_curve_from_points_doc},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef M_Curve_fit_nd_module_def = {
	PyModuleDef_HEAD_INIT,
	"curve_fit_nd",  /* m_name */
	M_Curve_fit_nd_doc,  /* m_doc */
	0,  /* m_size */
	M_Curve_fit_nd_methods,  /* m_methods */
	NULL,  /* m_reload */
	NULL,  /* m_traverse */
	NULL,  /* m_clear */
	NULL,  /* m_free */
};

PyMODINIT_FUNC PyInit_curve_fit_nd(void)
{
	PyObject *mod;

	mod = PyModule_Create(&M_Curve_fit_nd_module_def);

	return mod;
}

