#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
using namespace std;

extern int var_num_int;
extern int var_num_mon;

double var_dbl(string, string);
int var_int(string, string);
string var_str(string, string);
const char* var_cstr(string, string);

void var_init(string);

#endif
