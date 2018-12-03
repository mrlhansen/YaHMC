#ifndef LOGGER_H
#define LOGGER_H

#include <string>
using namespace std;

#define CRITICAL 0
#define WARNING  10
#define INFO     20
#define NOTICE   30
#define DEBUG    40

void lprintf(const char*, int, const char*, ...);
void logger_init(string, int);
void logger_exit();

#endif
