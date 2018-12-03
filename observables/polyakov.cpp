#include <polyakov.h>
#include <global.h>

complex avr_polyakov(int mu)
{
	int np[4] = {np_t, np_x, np_y, np_z};
	int id, size, ndx, norm;
	int ma, mb, mc, x[4];
	suNg *sbuf, *rbuf;
	complex sum;
	suNg loop;

	ma = (mu+1)%4;
	mb = (mu+2)%4;
	mc = (mu+3)%4;

	norm = global.vol4 / global.dim[mu];
	size = local.vol4 / local.dim[mu];
	sbuf = new suNg[size];
	rbuf = sbuf;
	ndx = 0;

	for(x[ma] = 0; x[ma] < local.dim[ma]; x[ma]++)
	for(x[mb] = 0; x[mb] < local.dim[mb]; x[mb]++)
	for(x[mc] = 0; x[mc] < local.dim[mc]; x[mc]++)
	{
		loop.unit();
		for(x[mu] = 0; x[mu] < local.dim[mu]; x[mu]++)
		{
			id = local_index(x[0], x[1], x[2], x[3]);
			loop *= link(id,mu);
		}
		sbuf[ndx++] = loop;
	}

	if(mpi_size > 1)
	{
		if(mpi_rank == 0)
		{
			rbuf = new suNg[mpi_size*size];
		}
		mp_gather(sbuf, rbuf, size * sizeof(suNg));
	}

	if(mpi_rank == 0)
	{
		for(x[ma] = 0; x[ma] < np[ma]; x[ma]++)
		for(x[mb] = 0; x[mb] < np[mb]; x[mb]++)
		for(x[mc] = 0; x[mc] < np[mc]; x[mc]++)
		{
			for(ndx = 0; ndx < size; ndx++)
			{
				loop.unit();
				for(x[mu] = 0; x[mu] < np[mu]; x[mu]++)
				{
					id = pid_from_coords(x[0], x[1], x[2], x[3]);
					loop *= rbuf[id*size+ndx];
				}
				sum += ctrace(loop);
			}
		}
	}

	if((mpi_size > 1) && (mpi_rank == 0))
	{
		delete[] rbuf;
		delete[] sbuf;
	}
	else
	{
		delete[] sbuf;
	}

	mp_global_sum(&sum, 2);
	sum /= NC*norm;
	return sum;
}
