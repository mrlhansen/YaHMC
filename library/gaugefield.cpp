#include <global.h>
#include <repr.h>
#include <bcs.h>
#include <wilsonloop.h>
#include <memory.h>
#include <cmath>
#include <cstring>

suNg *gfield_fund;
suNg *gfield_copy;
suNf *gfield_repr;
suNa *afield_momenta;
static suNg *gbuf;

void gauge_field_init()
{
	gfield_fund = suNg_allocate(4);
	gfield_copy = suNg_allocate(4);
	gfield_repr = suNf_allocate(4);
	afield_momenta = suNa_allocate(4);

	if((mpi_size > 1) && (mpi_rank == 0))
	{
		gbuf = new suNg[4*global.vol4];
	}
	else
	{
		gbuf = gfield_fund;
	}
}

void gauge_field_backup()
{
	long memsz = 4*outer.vol4*sizeof(suNg);
	memcpy(gfield_copy, gfield_fund, memsz);
}

void gauge_field_restore()
{
	long memsz = 4*outer.vol4*sizeof(suNg);
	memcpy(gfield_fund, gfield_copy, memsz);
}

void gauge_field_save(const char *filename)
{
	header_t info;
	int pid, id;

	// Gather configuration
	if(mpi_size > 1)
	{
		pid = 4 * local.vol4 * sizeof(suNg);
		mp_gather(gfield_fund, gbuf, pid);
	}

	// Information about the configuration
	info.ident = CFG_IDENT;
	info.version = CFG_VERSION;
	info.dim_t = global.dim_t;
	info.dim_x = global.dim_x;
	info.dim_y = global.dim_y;
	info.dim_z = global.dim_z;
	info.nc = NC;
	info.repr = REPR_ID;
	info.num_accept = num_accept;
	info.num_cnfg = num_cnfg;
	info.plaquette = avr_wilson(1,1);

	// Wait until data has been written
	if(mpi_rank)
	{
		mp_barrier();
		return;
	}

	// Write configuration
	FILE *fp = fopen(filename, "wb");

	if(fp != NULL)
	{
		fwrite(&info, sizeof(header_t), 1, fp);

		for(int t = 0; t < global.dim_t; t++)
		for(int x = 0; x < global.dim_x; x++)
		for(int y = 0; y < global.dim_y; y++)
		for(int z = 0; z < global.dim_z; z++)
		{
			global_index(t, x, y, z, &pid, &id);
			id = pid * local.vol4 + id;
			gbuf[4*id+0].write(fp);
			gbuf[4*id+1].write(fp);
			gbuf[4*id+2].write(fp);
			gbuf[4*id+3].write(fp);
		}

		fclose(fp);
		lprintf("GAUGEFIELD", INFO, "Configuration stored: %s", filename);
	}
	else
	{
		lprintf("GAUGEFIELD", CRITICAL, "Unable to open file: %s", filename);
	}

	// Data has been written
	mp_barrier();
}

void gauge_field_load(const char *filename)
{
	header_t info;
	int sz, id, pid;
	FILE *fp;

	fp = fopen(filename, "rb");

	if(fp != NULL)
	{
		sz = fread(&info, sizeof(header_t), 1, fp);
		sz  = (info.ident != CFG_IDENT);
		sz += (info.version != CFG_VERSION);
		sz += (info.dim_t != global.dim_t);
		sz += (info.dim_x != global.dim_x);
		sz += (info.dim_y != global.dim_y);
		sz += (info.dim_z != global.dim_z);
		sz += (info.nc != NC);
		sz += (info.repr != REPR_ID);

		if(sz)
		{
			lprintf("GAUGEFIELD", CRITICAL, "Configuration [%s] is not consistent with specified settings", filename);
		}

		// Information
		num_accept = info.num_accept;
		num_cnfg = info.num_cnfg;

		// Load configuration
		if(mpi_rank == 0)
		{
			for(int t = 0; t < global.dim_t; t++)
			for(int x = 0; x < global.dim_x; x++)
			for(int y = 0; y < global.dim_y; y++)
			for(int z = 0; z < global.dim_z; z++)
			{
				global_index(t, x, y, z, &pid, &id);
				id = pid * local.vol4 + id;
				gbuf[4*id+0].read(fp);
				gbuf[4*id+1].read(fp);
				gbuf[4*id+2].read(fp);
				gbuf[4*id+3].read(fp);
			}
		}

		// Close file
		fclose(fp);

		// Scatter configuration
		if(mpi_size > 1)
		{
			pid = 4 * local.vol4 * sizeof(suNg);
			mp_scatter(gbuf, gfield_fund, pid);
		}

		// Boundary conditions
		apply_bcs_on_gauge_field();

		// Represent gauge field
		represent_gauge_field();

		// Print log info
		info.plaquette -= avr_wilson(1,1);
		lprintf("GAUGEFIELD", INFO, "Configuration loaded: %s", filename);
		lprintf("GAUGEFIELD", INFO, "Plaquette difference: %1.6e", fabs(info.plaquette));
	}
	else
	{
		lprintf("GAUGEFIELD", CRITICAL, "Unable to open file: %s", filename);
	}
}

void gauge_field_set_random()
{
	// Initialize all links
	sites_for(id)
	{
		link(id,0).random();
		link(id,1).random();
		link(id,2).random();
		link(id,3).random();
	}

	// Boundary conditions
	apply_bcs_on_gauge_field();

	// Represent gauge field
	represent_gauge_field();
}

void gauge_field_set_unit()
{
	// Initialize all links
	sites_for(id)
	{
		link(id,0).unit();
		link(id,1).unit();
		link(id,2).unit();
		link(id,3).unit();
	}

	// Boundary conditions
	apply_bcs_on_gauge_field();

	// Represent gauge field
	represent_gauge_field();
}
