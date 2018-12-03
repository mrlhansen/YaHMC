#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
using namespace std;

#define spinor_for(id,sf) \
	for(int id = sf.offset; id < sf.offset + sf.sites; id++)

#define sites_for(id) \
	for(int id = 0; id < local.vol4; id++)

#define extended_sites_for(id) \
	for(int id = 0; id < outer.vol4; id++)

#define odd_sites_for(id) \
	for(int id = (local.vol4 / 2); id < local.vol4; id++)

#define even_sites_for(id) \
	for(int id = 0; id < (local.vol4 / 2); id++)

#define fw_index(id,mu) \
	site[id].fw[mu]

#define bk_index(id,mu) \
	site[id].bk[mu]

typedef struct {
	int fw[4];
	int bk[4];
	int coord[4];
	int parity;
} site_t;

typedef struct _dim {
	int &dim_t;
	int &dim_x;
	int &dim_y;
	int &dim_z;
	int dim[4];
	int vol3;
	int vol4;
	_dim(): dim_t(dim[0]), dim_x(dim[1]), dim_y(dim[2]), dim_z(dim[3]) {}
} dim_t;

typedef struct {
	int pid;        // Neighbor PID
	int tag_send;   // Tag for sending
	int tag_recv;   // Tag for receiving
	int size;       // Buffer size
	int *idx_send;  // Index for elements in send buffer
	int *idx_recv;  // Index for elements in recv buffer
	void *buf_send; // Send buffer
	void *buf_recv; // Recv buffer
} gd_block_t;

typedef vector<gd_block_t> gd_t;
extern gd_t gd_list[5];

extern site_t *site;
extern int volume;

extern dim_t global;
extern dim_t outer;
extern dim_t local;
extern dim_t inner;

extern int cp_t;
extern int cp_x;
extern int cp_y;
extern int cp_z;

extern int np_t;
extern int np_x;
extern int np_y;
extern int np_z;

int pid_from_coords(int, int, int, int);
int pid_from_offset(int, int, int, int, int);
int extended_index(int, int, int, int);
int extended_index_offset(int, int, int, int, int);
int extended_index_global(int, int, int, int);
int local_index(int, int, int, int);
int global_time(int);
void global_index(int, int, int, int, int*, int*);
void geometry_init(int);

#endif
