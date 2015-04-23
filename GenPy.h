#include "Python.h"
#include "stdlib.h"
#include <map>

void PyInit();

int GenPy(PUNGraph & G, ofstream& TFile, const TStr& args);