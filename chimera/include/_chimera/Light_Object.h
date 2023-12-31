#ifndef chimera_Light_object_h
# define chimera_Light_object_h
# if defined(_MSC_VER) && (_MSC_VER >= 1020)
#  pragma once
# endif

# define PY_SSIZE_T_CLEAN 1
#include <Python.h>
# include "Light.h"
#include <otf/WrapPy2.h>
namespace chimera {

CHIMERA_IMEX extern PyTypeObject Light_objectType;

// help subclasses initialize cached attributes
CHIMERA_IMEX extern PyObject* Light_attrInitcolor(PyObject*, void*);

struct Light_object: public PyObject {
	PyObject* _inst_dict;
	otf::WrapPyObj* _inst_data;
	Light* _inst() { return static_cast<Light*>(_inst_data); }
	PyObject* _weaklist;
};

CHIMERA_IMEX extern Light* getInst(Light_object* self);

} // namespace chimera

namespace otf {

template <> inline bool
WrapPyType<chimera::Light>::check(PyObject* _o, bool noneOk)
{
	if (noneOk && _o == Py_None)
		return true;
	return PyObject_TypeCheck(_o, &chimera::Light_objectType);
}

template <> CHIMERA_IMEX PyObject* pyObject(chimera::Light* _o);
template <> inline PyObject* pyObject(chimera::Light const* _o) { return pyObject(const_cast<chimera::Light*>(_o)); }

} // namespace otf

#endif
