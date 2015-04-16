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

int ParseArgument(const TStr& arg, const TStr& argType, PyObject** argPy)
{
	if (argType == "int")
	{
		*argPy = PyInt_FromLong(arg.GetInt());
	}
	else if (argType == "double")
	{
		*argPy = PyFloat_FromDouble(arg.GetFlt());
	}
	else if (argType == "string")
	{
		*argPy = PyString_FromString(arg.CStr());
	}
	if (!*argPy) {
        fprintf(stderr, "Cannot convert argument\n");
        return 0;
    }
	return 1;
}

int CallPyFunction(const char *moduleName, const char *funcName, const TStrV& args, const TStrV& argTypes,  PyObject** res, PyObject*** pyObjects = nullptr)
{
	PyObject *pName, *pModule, *pFunc, *pArgs;
	bool err = false;
	// get PyObject representation of moduleName
	pName = PyString_FromString(moduleName);
	// import module
	pModule = PyImport_Import(pName);
	// we don't need pName anymore
	//Py_DECREF(pName);
	// if module was loaded 
	if (pModule != nullptr) {
		// get pointer to function
		pFunc = PyObject_GetAttrString(pModule, funcName);
		// check function for existence
		if (pFunc && PyCallable_Check(pFunc)) {
			// a number of arguments
			int argc = argTypes.Len();
			// tuple of arguments
			pArgs = PyTuple_New(argc);
			// index of current PyObject in array of PyObjects
			int argPyObjIndex = 0, argStrIndex = 0;
			// parsing arguments
			for (size_t i = 0; i < argc; i++)
			{
				PyObject **arg = new PyObject*[1];
				//printf("argtypes[%d] = %s\n", i, argTypes[i].CStr());
				if (argTypes[i] != "pyobject"){
					if (!ParseArgument(args[argStrIndex++], argTypes[i], arg))
						err = true;
				}
				else
				{
					// non-safe
					*arg = pyObjects[argPyObjIndex++][0];
				}
				PyTuple_SetItem(pArgs, i, *arg);
				Py_DECREF(arg);
			}
			*res = PyObject_CallObject(pFunc, pArgs);
			if (res == nullptr)
			{
				PyErr_Print();
                fprintf(stderr,"Call failed\n");
                err = true;
			}
			//Py_DECREF(pArgs);
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
	//Py_DECREF(pArgs);
	Py_DECREF(pFunc);
	Py_DECREF(pModule);
	if (err == true) return 1;
	return 0;
}


void Get(PyObject **G, const char* funcName, PyObject***list)
{
	PyObject *** pyObjects = new PyObject**[1];
	pyObjects[0] = G;
	*list = new PyObject*[1];
	char * moduleName = "function";
	int argc = 1;
	TStrV args(0); TStrV argTypes(argc);
	argTypes[0] = "pyobject";
	CallPyFunction(moduleName, funcName, args, argTypes, *list, pyObjects);
	Py_DECREF(pyObjects);	
}

void GetNodes(PyObject **G, PyObject***list)
{
	return Get(G, "nodes", list);
}

void GetEdges(PyObject **G, PyObject***list)
{
	return Get(G, "edges", list);
}



int GenPy(PUNGraph &res)
{
	int argc = 2;
	PyObject **G = new PyObject*[1];
	// add path to NetworkX folders and initialize interpretator
	PyInit();
		
	char *moduleName = "random_graphs";
	char *funcName = "fast_gnp_random_graph";
	TStrV args(argc);
	args[0] = "1024"; args[1] = "0.003";
	TStrV argTypes(argc);
	argTypes[0] = "int"; argTypes[1] = "double";
	CallPyFunction(moduleName, funcName, args, argTypes, G);
	
	PyObject*** nodes = new PyObject**[1];
	GetNodes(G, nodes);
	int nodesCount = PyList_Size(*(nodes[0]));
	printf("nodesCount = %d, ", nodesCount);
	res = PUNGraph::TObj::New();
    res->Reserve(nodesCount, nodesCount*nodesCount);
	for (size_t i = 0; i < nodesCount; i++)
		res->AddNode(i);
	Py_DECREF(nodes);

	PyObject*** edges = new PyObject**[1];
	GetEdges(G, edges);
	int edgesCount = PyList_Size(*(edges[0]));
	printf("edgesCount = %d\n", edgesCount);
	for (size_t i = 0; i < edgesCount; i++)
	{
		PyObject* item = PySequence_Fast_GET_ITEM(*(edges[0]), i);
		int v1, v2;
		PyObject* node = PySequence_Fast_GET_ITEM(item,0);
		v1 = PyInt_AsLong(node);
		node = PySequence_Fast_GET_ITEM(item,1);
		v2 = PyInt_AsLong(node);
		res->AddEdge(v1,v2);
	}

	Py_DECREF(G);
	Py_DECREF(edges);
	Py_Finalize(); // очищение памяти, отданной интерпретатору
	return 0;
}