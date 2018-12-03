#ifndef RATIONAL_H
#define RATIONAL_H

#include <spinorfield.h>
#include <inverter.h>

typedef struct {
	int n;
	int d;
	int order;
	double prec;
	double min;
	double max;
	double a0;
	double pole[48];
	double root[48];
	double residue[48];
} rational_app;

void rational_find(rational_app&, int, int, double, double, double);
void rational_calc(rational_app&, mvm_operator, double, double, SpinorField&, SpinorField&);

#endif
