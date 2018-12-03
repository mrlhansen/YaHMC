#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
using namespace std;

double var_dbl(string);
void var_dbl(string, double);

int var_int(string);
void var_int(string, int);

string var_str(string);
void var_str(string, string);
const char* var_cstr(string);

void var_init(string);

#endif
