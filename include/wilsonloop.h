#ifndef WILSONLOOP_H
#define WILSONLOOP_H

#include <suN.h>

suNg clover_fund(int id, int mu, int nu);
suNf clover_repr(int id, int mu, int nu);
double wilson_loop(int, int, int, int, int);
double avr_wilson(int, int);

#endif
