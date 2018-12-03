#ifndef FORCES_H
#define FORCES_H

#include <spinorfield.h>

void clover_logdet_force(double, double);
void clover_fermion_force(double, SpinorField&, SpinorField&);
void wilson_fermion_force(double, double, double, SpinorField&, SpinorField&);
void wilson_fermion_force_q(double, double, double, SpinorField&);

void begin_fermion_force();
void end_fermion_force(double);

void wilson_gauge_force(double, double);
void improved_gauge_force(double, double, double, double);

#endif
