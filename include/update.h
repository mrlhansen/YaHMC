#ifndef UPDATE_H
#define UPDATE_H

void hmc_init();
bool detect_improved_gauge();

void update_momenta(double, int);
void update_links(double);
void update(double);

#endif
