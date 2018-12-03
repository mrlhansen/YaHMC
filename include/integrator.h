#ifndef INTEGRATOR_H
#define INTEGRATOR_H

typedef struct _ip {
	int level;
	int steps;
	void (*integrator)(double, struct _ip*);
	struct _ip *next;
} int_par;

void o2lf_multilevel(double, int_par*);
void o2mn_multilevel(double, int_par*);
void o4mn_multilevel(double, int_par*);

#endif
