#include <global.h>

// Memory layout:
// Standard: | inner | border | halo |
// Even/odd: | inner even | border even | inner odd | border odd | halo even | halo odd |

// List of communication types:
// Type 0 (gauge field): sides+corners, all sites, full border
// Type 1 (spinor field): sides, all sites, single border
// Type 2 (spinor field): sides, even sites, single border
// Type 3 (spinor field): sides, odd sites, single border
// Type 4 (clover force): sides+corners, all sites, single border

// Global variables
site_t *site;
gd_t gd_list[5];
int volume;

// Coordinates for local process
int cp_t;
int cp_x;
int cp_y;
int cp_z;

// Parallel dimensions
int np_t;
int np_x;
int np_y;
int np_z;

// Lattice dimensions
dim_t global;
dim_t outer;
dim_t local;
dim_t inner;

// Block size for borders
static int sz_local_border_t;
static int sz_local_border_x;
static int sz_local_border_y;
static int sz_local_border_z;

static int sz_outer_border_t;
static int sz_outer_border_x;
static int sz_outer_border_y;
static int sz_outer_border_z;

// Border size in all directions
static int border;
static int border_t;
static int border_x;
static int border_y;
static int border_z;

// Auxiliary functions
static int minimum(int a, int b)
{
	return (a < b) ? a : b;
}

static int site_parity(int t, int x, int y, int z)
{
	return (t + x + y + z) % 2;
}

// Communication functions
int pid_from_coords(int t, int x, int y, int z)
{
	return (t * np_x * np_y * np_z) + (x * np_y * np_z) + (y * np_z) + z;
}

int pid_from_offset(int id, int to, int xo, int yo, int zo)
{
	int t, x, y, z;

	z = id;
	t = z / (np_x * np_y * np_z);
	z = z % (np_x * np_y * np_z);
	x = z / (np_y * np_z);
	z = z % (np_y * np_z);
	y = z / np_z;
	z = z % np_z;

	while(to < 0) to += np_t;
	while(xo < 0) xo += np_x;
	while(yo < 0) yo += np_y;
	while(zo < 0) zo += np_z;

	t = (t + to) % np_t;
	x = (x + xo) % np_x;
	y = (y + yo) % np_y;
	z = (z + zo) % np_z;

	return pid_from_coords(t, x, y, z);
}

static void setup_comm_one(int ndx, int parity, int size)
{
	int rp, sp, factor;
	int bdr_t, sz_t, ts, to;
	int bdr_x, sz_x, xs, xo;
	int bdr_y, sz_y, ys, yo;
	int bdr_z, sz_z, zs, zo;

	gd_block_t gft, gbt;
	gd_block_t gfx, gbx;
	gd_block_t gfy, gby;
	gd_block_t gfz, gbz;

	bdr_t = minimum(border_t, size);
	bdr_x = minimum(border_x, size);
	bdr_y = minimum(border_y, size);
	bdr_z = minimum(border_z, size);

	sz_t = bdr_t * local.dim_x * local.dim_y * local.dim_z;
	sz_x = bdr_x * local.dim_y * local.dim_z * local.dim_t;
	sz_y = bdr_y * local.dim_z * local.dim_t * local.dim_x;
	sz_z = bdr_z * local.dim_t * local.dim_x * local.dim_y;

	if(parity == EVEN)
	{
		sp = EVEN;
		rp = ((size % 2) ? ODD : EVEN);
		factor = 2;
	}
	else if(parity == ODD)
	{
		sp = ODD;
		rp = ((size % 2) ? EVEN : ODD);
		factor = 2;
	}
	else
	{
		sp = BOTH;
		rp = BOTH;
		factor = 1;
	}

	sz_t /= factor;
	sz_x /= factor;
	sz_y /= factor;
	sz_z /= factor;

	if(sz_t)
	{
		gft.pid = pid_from_offset(mpi_rank, +1, 0, 0, 0);
		gft.size = sz_t;
		gft.tag_send = 0;
		gft.tag_recv = 1;
		gft.idx_send = new int[sz_t];
		gft.idx_recv = new int[sz_t];

		gbt.pid = pid_from_offset(mpi_rank, -1, 0, 0, 0);
		gbt.size = sz_t;
		gbt.tag_send = 1;
		gbt.tag_recv = 0;
		gbt.idx_send = new int[sz_t];
		gbt.idx_recv = new int[sz_t];

		gd_list[ndx].push_back(gft);
		gd_list[ndx].push_back(gbt);
	}

	if(sz_x)
	{
		gfx.pid = pid_from_offset(mpi_rank, 0, +1, 0, 0);
		gfx.size = sz_x;
		gfx.tag_send = 0;
		gfx.tag_recv = 1;
		gfx.idx_send = new int[sz_x];
		gfx.idx_recv = new int[sz_x];

		gbx.pid = pid_from_offset(mpi_rank, 0, -1, 0, 0);
		gbx.size = sz_x;
		gbx.tag_send = 1;
		gbx.tag_recv = 0;
		gbx.idx_send = new int[sz_x];
		gbx.idx_recv = new int[sz_x];

		gd_list[ndx].push_back(gfx);
		gd_list[ndx].push_back(gbx);
	}

	if(sz_y)
	{
		gfy.pid = pid_from_offset(mpi_rank, 0, 0, +1, 0);
		gfy.size = sz_y;
		gfy.tag_send = 0;
		gfy.tag_recv = 1;
		gfy.idx_send = new int[sz_y];
		gfy.idx_recv = new int[sz_y];

		gby.pid = pid_from_offset(mpi_rank, 0, 0, -1, 0);
		gby.size = sz_y;
		gby.tag_send = 1;
		gby.tag_recv = 0;
		gby.idx_send = new int[sz_y];
		gby.idx_recv = new int[sz_y];

		gd_list[ndx].push_back(gfy);
		gd_list[ndx].push_back(gby);
	}

	if(sz_z)
	{
		gfz.pid = pid_from_offset(mpi_rank, 0, 0, 0, +1);
		gfz.size = sz_z;
		gfz.tag_send = 0;
		gfz.tag_recv = 1;
		gfz.idx_send = new int[sz_z];
		gfz.idx_recv = new int[sz_z];

		gbz.pid = pid_from_offset(mpi_rank, 0, 0, 0, -1);
		gbz.size = sz_z;
		gbz.tag_send = 1;
		gbz.tag_recv = 0;
		gbz.idx_send = new int[sz_z];
		gbz.idx_recv = new int[sz_z];

		gd_list[ndx].push_back(gfz);
		gd_list[ndx].push_back(gbz);
	}

	for(int t = 0; t < local.dim_t; t++)
	for(int x = 0; x < local.dim_x; x++)
	for(int y = 0; y < local.dim_y; y++)
	for(int z = 0; z < local.dim_z; z++)
	{
		ts = t + border_t;
		xs = x + border_x;
		ys = y + border_y;
		zs = z + border_z;

		to = local.dim_t - t - 1;
		xo = local.dim_x - x - 1;
		yo = local.dim_y - y - 1;
		zo = local.dim_z - z - 1;

		parity = site_parity(ts, xs, ys, zs);
		parity = parity ? ODD : EVEN;

		if(t < bdr_t)
		{
			if(sp & parity)
			{
				*gbt.idx_send++ = extended_index(ts, xs, ys, zs);
			}
			if(rp & parity)
			{
				*gbt.idx_recv++ = extended_index(ts-bdr_t, xs, ys, zs);
			}
		}
		else if(to < bdr_t)
		{
			if(sp & parity)
			{
				*gft.idx_send++ = extended_index(ts, xs, ys, zs);
			}
			if(rp & parity)
			{
				*gft.idx_recv++ = extended_index(ts+bdr_t, xs, ys, zs);
			}
		}

		if(x < bdr_x)
		{
			if(sp & parity)
			{
				*gbx.idx_send++ = extended_index(ts, xs, ys, zs);
			}
			if(rp & parity)
			{
				*gbx.idx_recv++ = extended_index(ts, xs-bdr_x, ys, zs);
			}
		}
		else if(xo < bdr_x)
		{
			if(sp & parity)
			{
				*gfx.idx_send++ = extended_index(ts, xs, ys, zs);
			}
			if(rp & parity)
			{
				*gfx.idx_recv++ = extended_index(ts, xs+bdr_x, ys, zs);
			}
		}

		if(y < bdr_y)
		{
			if(sp & parity)
			{
				*gby.idx_send++ = extended_index(ts, xs, ys, zs);
			}
			if(rp & parity)
			{
				*gby.idx_recv++ = extended_index(ts, xs, ys-bdr_y, zs);
			}
		}
		else if(yo < bdr_y)
		{
			if(sp & parity)
			{
				*gfy.idx_send++ = extended_index(ts, xs, ys, zs);
			}
			if(rp & parity)
			{
				*gfy.idx_recv++ = extended_index(ts, xs, ys+bdr_y, zs);
			}
		}

		if(z < bdr_z)
		{
			if(sp & parity)
			{
				*gbz.idx_send++ = extended_index(ts, xs, ys, zs);
			}
			if(rp & parity)
			{
				*gbz.idx_recv++ = extended_index(ts, xs, ys, zs-bdr_z);
			}
		}
		else if(zo < bdr_z)
		{
			if(sp & parity)
			{
				*gfz.idx_send++ = extended_index(ts, xs, ys, zs);
			}
			if(rp & parity)
			{
				*gfz.idx_recv++ = extended_index(ts, xs, ys, zs+bdr_z);
			}
		}
	}
}

static void setup_comm_two(int ndx, int size)
{
	int bdr_t, ts, to;
	int bdr_x, xs, xo;
	int bdr_y, ys, yo;
	int bdr_z, zs, zo;
	int sz_tx, sz_ty, sz_tz;
	int sz_xy, sz_xz, sz_yz;

	gd_block_t nwtx, netx, setx, swtx;
	gd_block_t nwty, nety, sety, swty;
	gd_block_t nwtz, netz, setz, swtz;
	gd_block_t nwxy, nexy, sexy, swxy;
	gd_block_t nwxz, nexz, sexz, swxz;
	gd_block_t nwyz, neyz, seyz, swyz;

	bdr_t = minimum(border_t, size);
	bdr_x = minimum(border_x, size);
	bdr_y = minimum(border_y, size);
	bdr_z = minimum(border_z, size);

	sz_tx = bdr_t * bdr_x * local.dim_y * local.dim_z;
	sz_ty = bdr_t * bdr_y * local.dim_x * local.dim_z;
	sz_tz = bdr_t * bdr_z * local.dim_x * local.dim_y;
	sz_xy = bdr_x * bdr_y * local.dim_t * local.dim_z;
	sz_xz = bdr_x * bdr_z * local.dim_t * local.dim_y;
	sz_yz = bdr_y * bdr_z * local.dim_t * local.dim_x;

	if(sz_tx)
	{
		nwtx.pid = pid_from_offset(mpi_rank, -1, +1, 0, 0);
		nwtx.size = sz_tx;
		nwtx.tag_send = 0;
		nwtx.tag_recv = 1;
		nwtx.idx_send = new int[sz_tx];
		nwtx.idx_recv = new int[sz_tx];

		netx.pid = pid_from_offset(mpi_rank, +1, +1, 0, 0);
		netx.size = sz_tx;
		netx.tag_send = 0;
		netx.tag_recv = 1;
		netx.idx_send = new int[sz_tx];
		netx.idx_recv = new int[sz_tx];

		setx.pid = pid_from_offset(mpi_rank, +1, -1, 0, 0);
		setx.size = sz_tx;
		setx.tag_send = 1;
		setx.tag_recv = 0;
		setx.idx_send = new int[sz_tx];
		setx.idx_recv = new int[sz_tx];

		swtx.pid = pid_from_offset(mpi_rank, -1, -1, 0, 0);
		swtx.size = sz_tx;
		swtx.tag_send = 1;
		swtx.tag_recv = 0;
		swtx.idx_send = new int[sz_tx];
		swtx.idx_recv = new int[sz_tx];

		gd_list[ndx].push_back(nwtx);
		gd_list[ndx].push_back(netx);
		gd_list[ndx].push_back(setx);
		gd_list[ndx].push_back(swtx);
	}

	if(sz_ty)
	{
		nwty.pid = pid_from_offset(mpi_rank, -1, 0, +1, 0);
		nwty.size = sz_ty;
		nwty.tag_send = 0;
		nwty.tag_recv = 1;
		nwty.idx_send = new int[sz_ty];
		nwty.idx_recv = new int[sz_ty];

		nety.pid = pid_from_offset(mpi_rank, +1, 0, +1, 0);
		nety.size = sz_ty;
		nety.tag_send = 0;
		nety.tag_recv = 1;
		nety.idx_send = new int[sz_ty];
		nety.idx_recv = new int[sz_ty];

		sety.pid = pid_from_offset(mpi_rank, +1, 0, -1, 0);
		sety.size = sz_ty;
		sety.tag_send = 1;
		sety.tag_recv = 0;
		sety.idx_send = new int[sz_ty];
		sety.idx_recv = new int[sz_ty];

		swty.pid = pid_from_offset(mpi_rank, -1, 0, -1, 0);
		swty.size = sz_ty;
		swty.tag_send = 1;
		swty.tag_recv = 0;
		swty.idx_send = new int[sz_ty];
		swty.idx_recv = new int[sz_ty];

		gd_list[ndx].push_back(nwty);
		gd_list[ndx].push_back(nety);
		gd_list[ndx].push_back(sety);
		gd_list[ndx].push_back(swty);
	}

	if(sz_tz)
	{
		nwtz.pid = pid_from_offset(mpi_rank, -1, 0, 0, +1);
		nwtz.size = sz_tz;
		nwtz.tag_send = 0;
		nwtz.tag_recv = 1;
		nwtz.idx_send = new int[sz_tz];
		nwtz.idx_recv = new int[sz_tz];

		netz.pid = pid_from_offset(mpi_rank, +1, 0, 0, +1);
		netz.size = sz_tz;
		netz.tag_send = 0;
		netz.tag_recv = 1;
		netz.idx_send = new int[sz_tz];
		netz.idx_recv = new int[sz_tz];

		setz.pid = pid_from_offset(mpi_rank, +1, 0, 0, -1);
		setz.size = sz_tz;
		setz.tag_send = 1;
		setz.tag_recv = 0;
		setz.idx_send = new int[sz_tz];
		setz.idx_recv = new int[sz_tz];

		swtz.pid = pid_from_offset(mpi_rank, -1, 0, 0, -1);
		swtz.size = sz_tz;
		swtz.tag_send = 1;
		swtz.tag_recv = 0;
		swtz.idx_send = new int[sz_tz];
		swtz.idx_recv = new int[sz_tz];

		gd_list[ndx].push_back(nwtz);
		gd_list[ndx].push_back(netz);
		gd_list[ndx].push_back(setz);
		gd_list[ndx].push_back(swtz);
	}

	if(sz_xy)
	{
		nwxy.pid = pid_from_offset(mpi_rank, 0, -1, +1, 0);
		nwxy.size = sz_xy;
		nwxy.tag_send = 0;
		nwxy.tag_recv = 1;
		nwxy.idx_send = new int[sz_xy];
		nwxy.idx_recv = new int[sz_xy];

		nexy.pid = pid_from_offset(mpi_rank, 0, +1, +1, 0);
		nexy.size = sz_xy;
		nexy.tag_send = 0;
		nexy.tag_recv = 1;
		nexy.idx_send = new int[sz_xy];
		nexy.idx_recv = new int[sz_xy];

		sexy.pid = pid_from_offset(mpi_rank, 0, +1, -1, 0);
		sexy.size = sz_xy;
		sexy.tag_send = 1;
		sexy.tag_recv = 0;
		sexy.idx_send = new int[sz_xy];
		sexy.idx_recv = new int[sz_xy];

		swxy.pid = pid_from_offset(mpi_rank, 0, -1, -1, 0);
		swxy.size = sz_xy;
		swxy.tag_send = 1;
		swxy.tag_recv = 0;
		swxy.idx_send = new int[sz_xy];
		swxy.idx_recv = new int[sz_xy];

		gd_list[ndx].push_back(nwxy);
		gd_list[ndx].push_back(nexy);
		gd_list[ndx].push_back(sexy);
		gd_list[ndx].push_back(swxy);
	}

	if(sz_xz)
	{
		nwxz.pid = pid_from_offset(mpi_rank, 0, -1, 0, +1);
		nwxz.size = sz_xz;
		nwxz.tag_send = 0;
		nwxz.tag_recv = 1;
		nwxz.idx_send = new int[sz_xz];
		nwxz.idx_recv = new int[sz_xz];

		nexz.pid = pid_from_offset(mpi_rank, 0, +1, 0, +1);
		nexz.size = sz_xz;
		nexz.tag_send = 0;
		nexz.tag_recv = 1;
		nexz.idx_send = new int[sz_xz];
		nexz.idx_recv = new int[sz_xz];

		sexz.pid = pid_from_offset(mpi_rank, 0, +1, 0, -1);
		sexz.size = sz_xz;
		sexz.tag_send = 1;
		sexz.tag_recv = 0;
		sexz.idx_send = new int[sz_xz];
		sexz.idx_recv = new int[sz_xz];

		swxz.pid = pid_from_offset(mpi_rank, 0, -1, 0, -1);
		swxz.size = sz_xz;
		swxz.tag_send = 1;
		swxz.tag_recv = 0;
		swxz.idx_send = new int[sz_xz];
		swxz.idx_recv = new int[sz_xz];

		gd_list[ndx].push_back(nwxz);
		gd_list[ndx].push_back(nexz);
		gd_list[ndx].push_back(sexz);
		gd_list[ndx].push_back(swxz);
	}

	if(sz_yz)
	{
		nwyz.pid = pid_from_offset(mpi_rank, 0, 0, -1, +1);
		nwyz.size = sz_yz;
		nwyz.tag_send = 0;
		nwyz.tag_recv = 1;
		nwyz.idx_send = new int[sz_yz];
		nwyz.idx_recv = new int[sz_yz];

		neyz.pid = pid_from_offset(mpi_rank, 0, 0, +1, +1);
		neyz.size = sz_yz;
		neyz.tag_send = 0;
		neyz.tag_recv = 1;
		neyz.idx_send = new int[sz_yz];
		neyz.idx_recv = new int[sz_yz];

		seyz.pid = pid_from_offset(mpi_rank, 0, 0, +1, -1);
		seyz.size = sz_yz;
		seyz.tag_send = 1;
		seyz.tag_recv = 0;
		seyz.idx_send = new int[sz_yz];
		seyz.idx_recv = new int[sz_yz];

		swyz.pid = pid_from_offset(mpi_rank, 0, 0, -1, -1);
		swyz.size = sz_yz;
		swyz.tag_send = 1;
		swyz.tag_recv = 0;
		swyz.idx_send = new int[sz_yz];
		swyz.idx_recv = new int[sz_yz];

		gd_list[ndx].push_back(nwyz);
		gd_list[ndx].push_back(neyz);
		gd_list[ndx].push_back(seyz);
		gd_list[ndx].push_back(swyz);
	}

	for(int t = 0; t < local.dim_t; t++)
	for(int x = 0; x < local.dim_x; x++)
	for(int y = 0; y < local.dim_y; y++)
	for(int z = 0; z < local.dim_z; z++)
	{
		ts = t + border_t;
		xs = x + border_x;
		ys = y + border_y;
		zs = z + border_z;

		to = local.dim_t - t - 1;
		xo = local.dim_x - x - 1;
		yo = local.dim_y - y - 1;
		zo = local.dim_z - z - 1;

		if(t < bdr_t)
		{
			if(x < bdr_x)
			{
				*swtx.idx_send++ = extended_index(ts, xs, ys, zs);
				*swtx.idx_recv++ = extended_index(ts-bdr_t, xs-bdr_x, ys, zs);
			}
			else if(xo < bdr_x)
			{
				*nwtx.idx_send++ = extended_index(ts, xs, ys, zs);
				*nwtx.idx_recv++ = extended_index(ts-bdr_t, xs+bdr_x, ys, zs);
			}

			if(y < bdr_y)
			{
				*swty.idx_send++ = extended_index(ts, xs, ys, zs);
				*swty.idx_recv++ = extended_index(ts-bdr_t, xs, ys-bdr_y, zs);
			}
			else if(yo < bdr_y)
			{
				*nwty.idx_send++ = extended_index(ts, xs, ys, zs);
				*nwty.idx_recv++ = extended_index(ts-bdr_t, xs, ys+bdr_y, zs);
			}

			if(z < bdr_z)
			{
				*swtz.idx_send++ = extended_index(ts, xs, ys, zs);
				*swtz.idx_recv++ = extended_index(ts-bdr_t, xs, ys, zs-bdr_z);
			}
			else if(zo < bdr_z)
			{
				*nwtz.idx_send++ = extended_index(ts, xs, ys, zs);
				*nwtz.idx_recv++ = extended_index(ts-bdr_t, xs, ys, zs+bdr_z);
			}
		}
		else if(to < bdr_t)
		{
			if(x < bdr_x)
			{
				*setx.idx_send++ = extended_index(ts, xs, ys, zs);
				*setx.idx_recv++ = extended_index(ts+bdr_t, xs-bdr_x, ys, zs);
			}
			else if(xo < bdr_x)
			{
				*netx.idx_send++ = extended_index(ts, xs, ys, zs);
				*netx.idx_recv++ = extended_index(ts+bdr_t, xs+bdr_x, ys, zs);
			}

			if(y < bdr_y)
			{
				*sety.idx_send++ = extended_index(ts, xs, ys, zs);
				*sety.idx_recv++ = extended_index(ts+bdr_t, xs, ys-bdr_y, zs);
			}
			else if(yo < bdr_y)
			{
				*nety.idx_send++ = extended_index(ts, xs, ys, zs);
				*nety.idx_recv++ = extended_index(ts+bdr_t, xs, ys+bdr_y, zs);
			}

			if(z < bdr_z)
			{
				*setz.idx_send++ = extended_index(ts, xs, ys, zs);
				*setz.idx_recv++ = extended_index(ts+bdr_t, xs, ys, zs-bdr_z);
			}
			else if(zo < bdr_z)
			{
				*netz.idx_send++ = extended_index(ts, xs, ys, zs);
				*netz.idx_recv++ = extended_index(ts+bdr_t, xs, ys, zs+bdr_z);
			}
		}

		if(x < bdr_x)
		{
			if(y < bdr_y)
			{
				*swxy.idx_send++ = extended_index(ts, xs, ys, zs);
				*swxy.idx_recv++ = extended_index(ts, xs-bdr_x, ys-bdr_y, zs);
			}
			else if(yo < bdr_y)
			{
				*nwxy.idx_send++ = extended_index(ts, xs, ys, zs);
				*nwxy.idx_recv++ = extended_index(ts, xs-bdr_x, ys+bdr_y, zs);
			}

			if(z < bdr_z)
			{
				*swxz.idx_send++ = extended_index(ts, xs, ys, zs);
				*swxz.idx_recv++ = extended_index(ts, xs-bdr_x, ys, zs-bdr_z);
			}
			else if(zo < bdr_z)
			{
				*nwxz.idx_send++ = extended_index(ts, xs, ys, zs);
				*nwxz.idx_recv++ = extended_index(ts, xs-bdr_x, ys, zs+bdr_z);
			}
		}
		else if(xo < bdr_x)
		{
			if(y < bdr_y)
			{
				*sexy.idx_send++ = extended_index(ts, xs, ys, zs);
				*sexy.idx_recv++ = extended_index(ts, xs+bdr_x, ys-bdr_y, zs);
			}
			else if(yo < bdr_y)
			{
				*nexy.idx_send++ = extended_index(ts, xs, ys, zs);
				*nexy.idx_recv++ = extended_index(ts, xs+bdr_x, ys+bdr_y, zs);
			}

			if(z < bdr_z)
			{
				*sexz.idx_send++ = extended_index(ts, xs, ys, zs);
				*sexz.idx_recv++ = extended_index(ts, xs+bdr_x, ys, zs-bdr_z);
			}
			else if(zo < bdr_z)
			{
				*nexz.idx_send++ = extended_index(ts, xs, ys, zs);
				*nexz.idx_recv++ = extended_index(ts, xs+bdr_x, ys, zs+bdr_z);
			}
		}

		if(y < bdr_y)
		{
			if(z < bdr_z)
			{
				*swyz.idx_send++ = extended_index(ts, xs, ys, zs);
				*swyz.idx_recv++ = extended_index(ts, xs, ys-bdr_y, zs-bdr_z);
			}
			else if(zo < bdr_z)
			{
				*nwyz.idx_send++ = extended_index(ts, xs, ys, zs);
				*nwyz.idx_recv++ = extended_index(ts, xs, ys-bdr_y, zs+bdr_z);
			}
		}
		else if(yo < bdr_y)
		{
			if(z < bdr_z)
			{
				*seyz.idx_send++ = extended_index(ts, xs, ys, zs);
				*seyz.idx_recv++ = extended_index(ts, xs, ys+bdr_y, zs-bdr_z);
			}
			else if(zo < bdr_z)
			{
				*neyz.idx_send++ = extended_index(ts, xs, ys, zs);
				*neyz.idx_recv++ = extended_index(ts, xs, ys+bdr_y, zs+bdr_z);
			}
		}
	}
}

static void setup_communication()
{
	setup_comm_one(0, BOTH, border);
	setup_comm_two(0, border);
	setup_comm_one(1, BOTH, 1);
	setup_comm_one(2, EVEN, 1);
	setup_comm_one(3, ODD, 1);
	setup_comm_one(4, BOTH, 1);
	setup_comm_two(4, 1);
}

// Geometry functions
int extended_index(int t, int x, int y, int z)
{
	int offset = 0;
	int factor = 1;
	int hshift = 0;
	int lshift = 0;

#ifdef UPDATE_EO
	factor = 2;
	lshift = site_parity(t, x, y, z);
	hshift = local.vol4 + lshift * (outer.vol4 - local.vol4);
	lshift = lshift * local.vol4;
#endif

	// Halo in t-direction
	if(t < border_t || (outer.dim_t - t) <= border_t)
	{
		t = t % local.dim_t;
		offset += local.vol4 + hshift;
		offset += (t * outer.dim_x * outer.dim_y * outer.dim_z) + (x * outer.dim_y * outer.dim_z) + (y * outer.dim_z) + z;
		return (offset / factor);
	}
	else
	{
		t = t - border_t;
	}

	// Halo in x-direction
	if(x < border_x || (outer.dim_x - x) <= border_x)
	{
		x = x % local.dim_x;
		offset += local.vol4 + hshift + sz_outer_border_t;
		offset += (x * outer.dim_y * outer.dim_z * local.dim_t) + (y * outer.dim_z * local.dim_t) + (z * local.dim_t) + t;
		return (offset / factor);
	}
	else
	{
		x = x - border_x;
	}

	// Halo in y-direction
	if(y < border_y || (outer.dim_y - y) <= border_y)
	{
		y = y % local.dim_y;
		offset += local.vol4 + hshift + sz_outer_border_t + sz_outer_border_x;
		offset += (y * outer.dim_z * local.dim_t * local.dim_x) + (z * local.dim_t * local.dim_x) + (t * local.dim_x) + x;
		return (offset / factor);
	}
	else
	{
		y = y - border_y;
	}

	// Halo in z-direction
	if(z < border_z || (outer.dim_z - z) <= border_z)
	{
		z = z % local.dim_z;
		offset += local.vol4 + hshift + sz_outer_border_t + sz_outer_border_x + sz_outer_border_y;
		offset += (z * local.dim_t * local.dim_x * local.dim_y) + (t * local.dim_x * local.dim_y) + (x * local.dim_y) + y;
		return (offset / factor);
	}
	else
	{
		z = z - border_z;
	}

	// Border in t-direction
	if(t < border_t || (local.dim_t - t) <= border_t)
	{
		t = t % inner.dim_t;
		offset += inner.vol4 + lshift;
		offset += (t * local.dim_x * local.dim_y * local.dim_z) + (x * local.dim_y * local.dim_z) + (y * local.dim_z) + z;
		return (offset / factor);
	}
	else
	{
		t = t - border_t;
	}

	// Border in x-direction
	if(x < border_x || (local.dim_x - x) <= border_x)
	{
		x = x % inner.dim_x;
		offset += inner.vol4 + lshift + sz_local_border_t;
		offset += (x * local.dim_y * local.dim_z * inner.dim_t) + (y * local.dim_z * inner.dim_t) + (z * inner.dim_t) + t;
		return (offset / factor);
	}
	else
	{
		x = x - border_x;
	}

	// Border in y-direction
	if(y < border_y || (local.dim_y - y) <= border_y)
	{
		y = y % inner.dim_y;
		offset += inner.vol4 + lshift + sz_local_border_t + sz_local_border_x;
		offset += (y * local.dim_z * inner.dim_t * inner.dim_x) + (z * inner.dim_t * inner.dim_x) + (t * inner.dim_x) + x;
		return (offset / factor);
	}
	else
	{
		y = y - border_y;
	}

	// Border in z-direction
	if(z < border_z || (local.dim_z - z) <= border_z)
	{
		z = z % inner.dim_z;
		offset += inner.vol4 + lshift + sz_local_border_t + sz_local_border_x + sz_local_border_y;
		offset += (z * inner.dim_t * inner.dim_x * inner.dim_y) + (t * inner.dim_x * inner.dim_y) + (x * inner.dim_y) + y;
		return (offset / factor);
	}
	else
	{
		z = z - border_z;
	}

	// Inside inner lattice
	offset += lshift;
	offset += (t * inner.dim_x * inner.dim_y * inner.dim_z) + (x * inner.dim_y * inner.dim_z) + (y * inner.dim_z) + z;
	return (offset / factor);
}

int extended_index_offset(int id, int to, int xo, int yo, int zo)
{
	int t, x, y, z;

	t = site[id].coord[0];
	x = site[id].coord[1];
	y = site[id].coord[2];
	z = site[id].coord[3];

	t = (t + to + outer.dim_t) % outer.dim_t;
	x = (x + xo + outer.dim_x) % outer.dim_x;
	y = (y + yo + outer.dim_y) % outer.dim_y;
	z = (z + zo + outer.dim_z) % outer.dim_z;

	return extended_index(t, x, y, z);
}

int extended_index_global(int t, int x, int y, int z)
{
	int ta, tb;
	int xa, xb;
	int ya, yb;
	int za, zb;

	// t-direction
	ta = cp_t * local.dim_t - border_t;
	tb = ta + outer.dim_t;

	if(t >= ta && t < tb)
	{
		t = t - ta;
	}
	else if(ta < 0 && t >= (ta + global.dim_t))
	{
		t = t + border_t - global.dim_t;
	}
	else if(tb >= global.dim_t && t < (tb - global.dim_t))
	{
		t = t + border_t + local.dim_t;
	}
	else
	{
		return -1;
	}

	// x-direction
	xa = cp_x * local.dim_x - border_x;
	xb = xa + outer.dim_x;

	if(x >= xa && x < xb)
	{
		x = x - xa;
	}
	else if(xa < 0 && x >= (xa + global.dim_x))
	{
		x = x + border_x - global.dim_x;
	}
	else if(xb >= global.dim_x && x < (xb - global.dim_x))
	{
		x = x + border_x + local.dim_x;
	}
	else
	{
		return -1;
	}

	// y-direction
	ya = cp_y * local.dim_y - border_y;
	yb = ya + outer.dim_y;

	if(y >= ya && y < yb)
	{
		y = y - ya;
	}
	else if(ya < 0 && y >= (ya + global.dim_y))
	{
		y = y + border_y - global.dim_y;
	}
	else if(yb >= global.dim_y && y < (yb - global.dim_y))
	{
		y = y + border_y + local.dim_y;
	}
	else
	{
		return -1;
	}

	// z-direction
	za = cp_z * local.dim_z - border_z;
	zb = za + outer.dim_z;

	if(z >= za && z < zb)
	{
		z = z - za;
	}
	else if(za < 0 && z >= (za + global.dim_z))
	{
		z = z + border_z - global.dim_z;
	}
	else if(zb >= global.dim_z && z < (zb - global.dim_z))
	{
		z = z + border_z + local.dim_z;
	}
	else
	{
		return -1;
	}

	// Return local site
	return extended_index(t, x, y, z);
}

int local_index(int t, int x, int y, int z)
{
	t += border_t;
	x += border_x;
	y += border_y;
	z += border_z;
	return extended_index(t, x, y, z);
}

int global_time(int id)
{
	id = site[id].coord[0];
	id = cp_t * local.dim_t - border_t + id;
	id = (id + global.dim_t) % global.dim_t;
	return id;
}

void global_index(int t, int x, int y, int z, int *pid, int *id)
{
	int pt, lt;
	int px, lx;
	int py, ly;
	int pz, lz;

	pt = t / local.dim_t;
	lt = (t % local.dim_t) + border_t;

	px = x / local.dim_x;
	lx = (x % local.dim_x) + border_x;

	py = y / local.dim_y;
	ly = (y % local.dim_y) + border_y;

	pz = z / local.dim_z;
	lz = (z % local.dim_z) + border_z;

	*pid = pid_from_coords(pt, px, py, pz);
	*id = extended_index(lt, lx, ly, lz);
}

void geometry_init(int bdr)
{
	// Border size
	border = bdr;

	// Global dimensions
	global.dim_t = var_int("lat:dim_t");
	global.dim_x = var_int("lat:dim_x");
	global.dim_y = var_int("lat:dim_y");
	global.dim_z = var_int("lat:dim_z");
	global.vol3 = global.dim_x * global.dim_y * global.dim_z;
	global.vol4 = global.dim_t * global.vol3;

	// Parallel dimensions
#ifdef ENABLE_MPI
	np_t = var_int("mp:np_t");
	np_x = var_int("mp:np_x");
	np_y = var_int("mp:np_y");
	np_z = var_int("mp:np_z");
#else
	np_t = 1;
	np_x = 1;
	np_y = 1;
	np_z = 1;
#endif

	// Check grid size
	if((np_t < 1) || (np_x < 1) || (np_y < 1) || (np_z < 1))
	{
		lprintf("GEOMETRY", CRITICAL, "Invalid grid size (must be positive)");
	}

	// Check dimensions of global lattice
	if((global.dim_t % 2) || (global.dim_x % 2) || (global.dim_y % 2) || (global.dim_z % 2))
	{
		lprintf("GEOMETRY", CRITICAL, "Even dimensions required for global lattice");
	}

	// Check for proper parallelization
	if((global.dim_t % np_t) || (global.dim_x % np_x) || (global.dim_y % np_y) || (global.dim_z % np_z))
	{
		lprintf("GEOMETRY", CRITICAL, "The parallelization grid cannot evenly divide the global lattice");
	}

	// Borders
	border_t = (np_t > 1) * border;
	border_x = (np_x > 1) * border;
	border_y = (np_y > 1) * border;
	border_z = (np_z > 1) * border;

	// Local dimensions
	local.dim_t = global.dim_t / np_t;
	local.dim_x = global.dim_x / np_x;
	local.dim_y = global.dim_y / np_y;
	local.dim_z = global.dim_z / np_z;
	local.vol3 = local.dim_x * local.dim_y * local.dim_z;
	local.vol4 = local.dim_t * local.vol3;

	// Outer dimensions
	outer.dim_t = local.dim_t + (2 * border_t);
	outer.dim_x = local.dim_x + (2 * border_x);
	outer.dim_y = local.dim_y + (2 * border_y);
	outer.dim_z = local.dim_z + (2 * border_z);
	outer.vol3 = outer.dim_x * outer.dim_y * outer.dim_z;
	outer.vol4 = outer.dim_t * outer.vol3;

	// Inner dimensions
	inner.dim_t = local.dim_t - (2 * border_t);
	inner.dim_x = local.dim_x - (2 * border_x);
	inner.dim_y = local.dim_y - (2 * border_y);
	inner.dim_z = local.dim_z - (2 * border_z);
	inner.vol3 = inner.dim_x * inner.dim_y * inner.dim_z;
	inner.vol4 = inner.dim_t * inner.vol3;

	// Block size for borders
	sz_outer_border_t = 2 * border_t * outer.dim_x * outer.dim_y * outer.dim_z;
	sz_outer_border_x = 2 * border_x * local.dim_t * outer.dim_y * outer.dim_z;
	sz_outer_border_y = 2 * border_y * local.dim_t * local.dim_x * outer.dim_z;
	sz_outer_border_z = 2 * border_z * local.dim_t * local.dim_x * local.dim_y;

	sz_local_border_t = 2 * border_t * local.dim_x * local.dim_y * local.dim_z;
	sz_local_border_x = 2 * border_x * inner.dim_t * local.dim_y * local.dim_z;
	sz_local_border_y = 2 * border_y * inner.dim_t * inner.dim_x * local.dim_z;
	sz_local_border_z = 2 * border_z * inner.dim_t * inner.dim_x * inner.dim_y;

	// Print log info
	lprintf("GEOMETRY", INFO, "Global size: %dx%dx%dx%d", global.dim_t, global.dim_x, global.dim_y, global.dim_z);
	lprintf("GEOMETRY", INFO, "Grid size: %dx%dx%dx%d", np_t, np_x, np_y, np_z);
	lprintf("GEOMETRY", INFO, "Local size: %dx%dx%dx%d", local.dim_t, local.dim_x, local.dim_y, local.dim_z);
	lprintf("GEOMETRY", INFO, "Extended size: %dx%dx%dx%d", outer.dim_t, outer.dim_x, outer.dim_y, outer.dim_z);
	lprintf("GEOMETRY", INFO, "Border size: %d", border);

#ifdef UPDATE_EO
	lprintf("GEOMETRY", INFO, "Using even-odd preconditioning");
#endif

	// Check MPI world size
	if(mpi_size != np_t * np_z * np_y * np_x)
	{
		lprintf("GEOMETRY", CRITICAL, "The grid size does not match the number of MPI processes");
	}

	// Check dimensions of local lattice
	if((local.dim_t % 2) || (local.dim_x % 2) || (local.dim_y % 2) || (local.dim_z % 2))
	{
		lprintf("GEOMETRY", CRITICAL, "Even dimensions required for local lattice");
	}

	// Check local block size
	if((local.dim_t < 8) || (local.dim_x < 8) || (local.dim_y < 8) || (local.dim_z < 8))
	{
		lprintf("GEOMETRY", CRITICAL, "Local lattice is too small (minimum size is 8x8x8x8)");
	}

	// Allocate space
	site = new site_t[outer.vol4];

	// Adjusted volume
	volume = global.vol4;

	// Coordinates for local process
	cp_z = mpi_rank;
	cp_t = cp_z / (np_x * np_y * np_z);
	cp_z = cp_z % (np_x * np_y * np_z);
	cp_x = cp_z / (np_y * np_z);
	cp_z = cp_z % (np_y * np_z);
	cp_y = cp_z / np_z;
	cp_z = cp_z % np_z;

	// Calculate site neighbors
	for(int t = 0; t < outer.dim_t; t++)
	for(int x = 0; x < outer.dim_x; x++)
	for(int y = 0; y < outer.dim_y; y++)
	for(int z = 0; z < outer.dim_z; z++)
	{
		int id = extended_index(t, x, y, z);

		// Store coordinates
		site[id].coord[0] = t;
		site[id].coord[1] = x;
		site[id].coord[2] = y;
		site[id].coord[3] = z;

		// Store neighbours
		site[id].fw[0] = extended_index_offset(id, +1, 0, 0, 0);
		site[id].fw[1] = extended_index_offset(id, 0, +1, 0, 0);
		site[id].fw[2] = extended_index_offset(id, 0, 0, +1, 0);
		site[id].fw[3] = extended_index_offset(id, 0, 0, 0, +1);
		site[id].bk[0] = extended_index_offset(id, -1, 0, 0, 0);
		site[id].bk[1] = extended_index_offset(id, 0, -1, 0, 0);
		site[id].bk[2] = extended_index_offset(id, 0, 0, -1, 0);
		site[id].bk[3] = extended_index_offset(id, 0, 0, 0, -1);

		// Store parity
		site[id].parity = site_parity(t, x, y, z) ? ODD : EVEN;
	}

	// Communication
	setup_communication();
}
