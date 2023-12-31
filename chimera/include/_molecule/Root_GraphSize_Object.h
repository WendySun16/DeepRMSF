#ifndef molecule_Root_GraphSize_object_h
# define molecule_Root_GraphSize_object_h
# if defined(_MSC_VER) && (_MSC_VER >= 1020)
#  pragma once
# endif

# define PY_SSIZE_T_CLEAN 1
#include <Python.h>
# include "Mol.h"
#include <otf/WrapPy2.h>
namespace molecule {

MOLECULE_IMEX extern PyTypeObject Root_GraphSize_objectType;

struct Root_GraphSize_object: public PyObject {
	double _inst_data[(sizeof (Root::GraphSize) + sizeof (double) - 1) / sizeof (double)];
	bool _initialized;
	Root::GraphSize* _inst() { return reinterpret_cast<Root::GraphSize*>(_inst_data); }
};

MOLECULE_IMEX extern Root::GraphSize* getInst(Root_GraphSize_object* self);

} // namespace molecule

namespace otf {

template <> inline bool
WrapPyType<molecule::Root::GraphSize>::check(PyObject* _o, bool noneOk)
{
	if (noneOk && _o == Py_None)
		return true;
	return PyObject_TypeCheck(_o, &molecule::Root_GraphSize_objectType);
}

template <> MOLECULE_IMEX PyObject* pyObject(molecule::Root::GraphSize* _o);
template <> inline PyObject* pyObject(molecule::Root::GraphSize _o) { return pyObject(&_o); }
template <> inline PyObject* pyObject(molecule::Root::GraphSize const* _o) { return pyObject(const_cast<molecule::Root::GraphSize*>(_o)); }

} // namespace otf

#endif
