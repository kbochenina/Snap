#include "stdafx.h"
#include "GenPy.h"
#include <stdlib.h>

void AddPath(const char * path)
{
	PyObject * pName = PyString_FromString(path), *syspath;
	// reference to Python search path
	syspath = PySys_GetObject("path");
	// add path to syspath
	if (PyList_Insert(syspath, 0, pName))
		printf("Error inserting extra path into sys.path list\n");
	// reset sys.path object
	if (PySys_SetObject("path", syspath))
		printf("Error setting sys.path object\n");
}

void PyInit()
{
	char *path = "D:\\ITMO\\Graphs\\Software\\NetworkX\\networkx-1.9.1\\networkx\\generators";
	Py_Initialize(); // инициализация интерпретатора 
	AddPath(path);
	path = "D:\\ITMO\\Graphs\\Software\\NetworkX\\networkx-1.9.1";
	AddPath(path);
	path = "D:\\ITMO\\Graphs\\Software\\NetworkX\\networkx-1.9.1\\networkx\\readwrite";
	AddPath(path);
	path = "D:\\ITMO\\Graphs\\Software\\NetworkX\\networkx-1.9.1\\networkx\\classes";
	AddPath(path);
}

int ParseArgument(const TStr& arg, const TStr& argType, PyObject* argPy)
{
	if (argType == "int")
	{
		argPy = PyInt_FromLong(arg.GetInt());
	}
	else if (argType == "double")
	{
		argPy = PyFloat_FromDouble(arg.GetFlt());
	}
	else if (argType == "string")
	{
		argPy = PyString_FromString(arg.CStr());
	}
	if (!argPy) {
        fprintf(stderr, "Cannot convert argument\n");
        return 1;
    }
	return 0;
}

int CallPyFunction(const char *moduleName, const char *funcName, const TStrV& args, const TStrV& argTypes, PyObject** pyObjects, PyObject* res)
{
	PyObject *pName, *pModule, *pFunc, *pArgs;
	bool err = false;
	// get PyObject representation of moduleName
	pName = PyString_FromString(moduleName);
	// import module
	pModule = PyImport_Import(pName);
	// we don't need pName anymore
	Py_DECREF(pName);
	// if module was loaded 
	if (pModule != nullptr) {
		// get pointer to function
		pFunc = PyObject_GetAttrString(pModule, funcName);
		// check function for existence
		if (pFunc && PyCallable_Check(pFunc)) {
			// a number of arguments
			int argc = args.Len();
			// tuple of arguments
			pArgs = PyTuple_New(argc);
			// index of current PyObject in array of PyObjects
			int argPyObjIndex = 0;
			// parsing arguments
			for (size_t i = 0; i < argc; i++)
			{
				PyObject* arg;
				if (argTypes[i] != "pyobject"){
					arg = new PyObject;
					if (!ParseArgument(args[i], argTypes[i], arg))
						err = true;
				}
				else
				{
					// non-safe
					arg = pyObjects[argPyObjIndex++];
				}
				PyTuple_SetItem(pArgs, i, arg);
			}
			res = PyObject_CallObject(pFunc, pArgs);
			if (res == nullptr)
			{
				PyErr_Print();
                fprintf(stderr,"Call failed\n");
                err = true;
			}
		}
		else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\" in module \"%s\"\n", funcName, moduleName);
			err = true;
        }
	}
	else {
		PyErr_Print();
		fprintf(stderr, "Failed to load \"%s\"\n", moduleName);
		err = true;
    }
	Py_DECREF(pArgs);
	Py_DECREF(pFunc);
	Py_DECREF(pModule);
	if (err == true) return 1;
	return 0;
}

int GenPy()
{
	TStr genName("random_graphs");
	char *funcName = "fast_gnp_random_graph";
	int argc = 2;
	PyObject *pName, *pModule, *pFunc;
	PyObject *pArgs, *pValue, *pGraph;
	// add path to NetworkX folders and initialize interpretator
	PyInit();

	
	pName = PyString_FromString(genName.CStr());
	// import module with name GenName
	pModule = PyImport_Import(pName);
	Py_DECREF(pName);
	if (pModule != nullptr) {
        pFunc = PyObject_GetAttrString(pModule, funcName);
		if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(argc);
			pValue = PyInt_FromLong(1024);
			
			if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }	
			PyTuple_SetItem(pArgs, 0, pValue);
			pValue = PyFloat_FromDouble(0.25);
			if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }	
			PyTuple_SetItem(pArgs, 1, pValue);
           
            pGraph = PyObject_CallObject(pFunc, pArgs);

			
			pName = PyString_FromString("pajek");
			pModule = PyImport_Import(pName);
			if (pModule != nullptr)
			{
				funcName = "write_pajek";
				pFunc = PyObject_GetAttrString(pModule, funcName);
				if (pFunc && PyCallable_Check(pFunc))
				{
					pArgs = PyTuple_New(2);
					PyTuple_SetItem(pArgs, 0, pGraph);
					pValue = PyString_FromString("random_graph.dat");
					PyTuple_SetItem(pArgs, 1, pValue);
					PyObject_CallObject(pFunc, pArgs);
				}
			}
			
	        

			//PyObject_CallFunction(pValue, "is_directed");
			
			

            Py_DECREF(pArgs);
            if (pValue != nullptr) {
                printf("Result of call: %ld\n", PyInt_AsLong(pValue));
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", genName.CStr());
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
	}

	else {
		PyErr_Print();
		fprintf(stderr, "Failed to load \"%s\"\n", genName.CStr());
    }

	Py_Finalize(); // очищение памяти, отданной интерпретатору
	system("pause");
	return 0;
}