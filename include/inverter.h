#ifndef INVERTER_H
#define INVERTER_H

#include <spinorfield.h>

typedef void (*mvm_operator)(double, SpinorField&, SpinorField&);

typedef struct _inv_par {
	int num;
	double prec;
	double mass;
	double *shift;
	int log;
	mvm_operator mvm;
	_inv_par(): num(0), prec(0), mass(0), shift(0), log(0), mvm(0) {}
} inv_par;

int bicgstab_inverter(inv_par&, SpinorField&, SpinorField&);
int cg_inverter(inv_par&, SpinorField&, SpinorField&);
int cg_mshift(inv_par&, SpinorField*, SpinorField&);
int minres_inverter(inv_par&, SpinorField&, SpinorField&);
int auto_inverter(inv_par&, SpinorField&, SpinorField&);

#endif
