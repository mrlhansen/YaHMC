#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <complex.h>
#include <spinorfield.h>

#define SEMWALL 0x01
#define POINT   0x02

typedef struct {
	int col[4];
	complex val[4];
} dirac_t;

// propagator.cpp
int random_timeslice();
void create_semwall_source(int, int, SpinorField&);
void create_point_source(int, int, SpinorField&);
void calc_propagator(double, double, SpinorField&, SpinorField&);

// mesons.cpp
void meson_init(double, double, int, int);
void meson_measure();

// sf_mpcac.cpp
void sf_mpcac_measure(double);

#endif
