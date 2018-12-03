#ifndef WILSONFLOW_H
#define WILSONFLOW_H

#include <spinorfield.h>

void wf_zeta(double);
void wf_update(double);
void wf_iterate(double);
double wf_e_plaq();
double wf_e_sym();
double wf_charge();
void wf_init();

void wf_fermion_iterate(double);
void wf_fermion_init(SpinorField*, int);
void wf_fermion_free();

#endif
