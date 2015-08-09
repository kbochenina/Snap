#include "stdafx.h"



int main(int argc, char* argv[]) {
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Two stage generation. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	Try
	const TStr SettingsFNm = Env.GetIfArgPrefixStr("-s:", "parameters.txt", "Settings filename");
	vector<TStr> commandLineArgs;
	// read command line arguments for all generators
	ReadParameters(SettingsFNm, commandLineArgs);
	KroneckerByConf(commandLineArgs);
	
	Catch
    system("pause");
	return 0;
}
