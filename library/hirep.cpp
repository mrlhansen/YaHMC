#include <global.h>
#include <repr.h>
#include <bcs.h>
#include <wilsonloop.h>
#include <cstring>

static int init = 0;
static suNg *gbuf;

static void hirep_init()
{
	if(init == 0)
	{
		if((mpi_size > 1) && (mpi_rank == 0))
		{
			gbuf = new suNg[4*global.vol4];
		}
		else
		{
			gbuf = gfield_fund;
		}
		init = 1;
	}
}

static void swap_bytes(void *ptr, int n)
{
	char *p = (char*)ptr;
	int lo = 0;
	int hi = n-1;

	while(hi > lo)
	{
		char tmp = p[lo];
		p[lo++] = p[hi];
		p[hi--] = tmp;
	}
}

static void swap_double(double *ptr, int count)
{
	for(int n = 0; n < count; n++)
	{
		swap_bytes(ptr+n, sizeof(double));
	}
}

static void swap_int(int *ptr, int count)
{
	for(int n = 0; n < count; n++)
	{
		swap_bytes(ptr+n, sizeof(int));
	}
}

void gauge_field_load_hirep(const char *filename)
{
	int sz, id, pid;
	int info[5];
	double plaq;
	FILE *fp;
	char *ptr;
	suNg tmp;

	hirep_init();
	fp = fopen(filename, "rb");

	if(fp != NULL)
	{
		sz = fread(info, sizeof(int), 5, fp);
		sz = fread(&plaq, sizeof(double), 1, fp);
		swap_int(info, 5);
		swap_double(&plaq, 1);

		sz  = (info[0] != NC);
		sz += (info[1] != global.dim_t);
		sz += (info[2] != global.dim_x);
		sz += (info[3] != global.dim_y);
		sz += (info[4] != global.dim_z);

		if(sz)
		{
			lprintf("HIREP", CRITICAL, "Configuration [%s] is not consistent with specified settings", filename);
		}

		// Information
		num_accept = 0;
		num_cnfg = 0;
		ptr = strrchr(filename, 'n');

		if(ptr)
		{
			sscanf(ptr, "n%d", &num_cnfg);
		}

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
				for(int mu = 0; mu < 4; mu++)
				{
					tmp.read(fp);
					swap_double(tmp.re, NC*NC);
					swap_double(tmp.im, NC*NC);
					gbuf[4*id+mu] = tmp;
				}
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
		plaq -= avr_wilson(1,1);
		lprintf("HIREP", INFO, "Configuration loaded: %s", filename);
		lprintf("HIREP", INFO, "Plaquette difference: %1.6e", fabs(plaq));
	}
	else
	{
		lprintf("HIREP", CRITICAL, "Unable to open file: %s", filename);
	}
}


void gauge_field_save_hirep(const char *filename)
{
	int info[5];
	double plaq;
	int pid, id;
	suNg tmp;

	// Allocate auxilliary field
	hirep_init();

	// Gather configuration
	if(mpi_size > 1)
	{
		pid = 4 * local.vol4 * sizeof(suNg);
		mp_gather(gfield_fund, gbuf, pid);
	}

	// Information about the configuration
	info[0] = NC;
	info[1] = global.dim_t;
	info[2] = global.dim_x;
	info[3] = global.dim_y;
	info[4] = global.dim_z;
	plaq = avr_wilson(1,1);

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
		swap_int(info, 5);
		swap_double(&plaq, 1);
		fwrite(info, sizeof(int), 5, fp);
		fwrite(&plaq, sizeof(double), 1, fp);

		for(int t = 0; t < global.dim_t; t++)
		for(int x = 0; x < global.dim_x; x++)
		for(int y = 0; y < global.dim_y; y++)
		for(int z = 0; z < global.dim_z; z++)
		{
			global_index(t, x, y, z, &pid, &id);
			id = pid * local.vol4 + id;
			for(int mu = 0; mu < 4; mu++)
			{
				tmp = gbuf[4*id+mu];
				swap_double(tmp.re, NC*NC);
				swap_double(tmp.im, NC*NC);
				tmp.write(fp);
			}
		}

		fclose(fp);
		lprintf("HIREP", INFO, "Configuration stored: %s", filename);
	}
	else
	{
		lprintf("HIREP", CRITICAL, "Unable to open file: %s", filename);
	}

	// Data has been written
	mp_barrier();
}
