#include "stdafx.h"

void Error(const TStr& FuncName, const TStr& ErrorMsg){
	printf("Error in %s: %s\n", FuncName.CStr(), ErrorMsg.CStr());
	system("pause");
	exit(1);
}