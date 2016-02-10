#include "tbconfig.h"
#include "stdlib.h"
#include "makeduts.h"

void setupDir(char* tosetup);
//void resetConfig(TbConfig &config);
void setupRunProc(TbConfig &config, void (*choiceSet)(TbConfig &, bool, bool, bool, bool, char*));

void parseArgs(int argc, char* argv[], TbConfig &config);
void printUsage(char* argv[]);
void rangeSupplied(bool gotRange, char* argc[]);
