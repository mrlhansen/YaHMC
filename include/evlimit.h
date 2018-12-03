#ifndef EVLIMIT_H
#define EVLIMIT_H

#include <inverter.h>

double smallest_eigenvalue(double, mvm_operator);
double largest_eigenvalue(double, mvm_operator);
void evlimits(mvm_operator, double, double*, double*);

#endif
