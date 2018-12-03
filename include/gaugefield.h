#ifndef GAUGEFIELD_H
#define GAUGEFIELD_H

#include <suN.h>

#define link(id,mu) \
	gfield_fund[4*(id)+mu]

#define fermion_link(id,mu) \
	gfield_repr[4*(id)+mu]

#define momentum(id,mu) \
	afield_momenta[4*(id)+mu]

extern suNg *gfield_fund;
extern suNg *gfield_copy;
extern suNf *gfield_repr;
extern suNa *afield_momenta;

#define CFG_VERSION 0x20
#define CFG_IDENT 0x47464359

typedef struct {
	int ident;
	int version;
	int dim_t;
	int dim_x;
	int dim_y;
	int dim_z;
	int nc;
	int repr;
	int num_cnfg;
	int num_accept;
	double plaquette;
	char reserved[208];
} __attribute__((packed)) header_t;

void gauge_field_init();
void gauge_field_backup();
void gauge_field_restore();
void gauge_field_save(const char*);
void gauge_field_load(const char*);
void gauge_field_set_random();
void gauge_field_set_unit();

#endif
