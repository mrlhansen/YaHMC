#include <global.h>
#include <clover.h>
#include <timing.h>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#ifdef ENABLE_OMP
#include <omp.h>
#endif

int mpi_rank;
int mpi_size;

void mp_init(int argc, char *argv[])
{
	#ifdef ENABLE_MPI

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	atexit(mp_finalize);

	#else

	mpi_size = 1;
	mpi_rank = 0;

	#endif
}

void mp_setup()
{
	#ifdef ENABLE_OMP

	int threads = var_int("parallel", "threads");

	// Number of OpenMP threads
	if(threads)
	{
		omp_set_num_threads(threads);
	}
	else
	{
		#pragma omp parallel
		{
			#pragma omp master
			threads = omp_get_num_threads();
		}
	}

	// Print log info
	lprintf("MP", INFO, "OpenMP enabled using %d threads", threads);

	#endif
	#ifdef ENABLE_MPI

	// Print log info
	lprintf("MP", INFO, "MPI enabled using %d processes", mpi_size);

	#endif
}

#ifdef ENABLE_MPI

void mp_finalize()
{
	MPI_Finalize();
}

template <typename Type>
static void sync_send_buffer(gd_t &gd, Type *gf, int cnt)
{
	Type *buf;
	for(auto& blk : gd)
	{
		buf = (Type*)blk.buf_send;
		for(int n = 0; n < blk.size; n++)
		{
			int id = blk.idx_send[n];
			for(int mu = 0; mu < cnt; mu++)
			{
				buf[cnt*n+mu] = gf[cnt*id+mu];
			}
		}
	}
}

template <typename Type>
static void sync_recv_buffer(gd_t &gd, Type *gf, int cnt)
{
	Type *buf;
	for(auto& blk : gd)
	{
		buf = (Type*)blk.buf_recv;
		for(int n = 0; n < blk.size; n++)
		{
			int id = blk.idx_recv[n];
			for(int mu = 0; mu < cnt; mu++)
			{
				gf[cnt*id+mu] = buf[cnt*n+mu];
			}
		}
	}
}

void mp_transfer_spinor_even(SpinorField &spinor)
{
	static MPI_Request *status;
	static int init = 0;
	static gd_t gd;
	int num_status;
	Spinor *sf = spinor.ptr;

	if(init == 0)
	{
		gd = gd_list[2];
		for(auto& blk : gd)
		{
			blk.buf_send = new Spinor[blk.size];
			blk.buf_recv = new Spinor[blk.size];
		}
		status = new MPI_Request[2*gd.size()];
		init = 1;
	}

	num_status = 0;
	timing_start(tm_mpi);
	sync_send_buffer(gd, sf, 1);

	for(auto& blk : gd)
	{
		MPI_Isend(blk.buf_send, blk.size*sizeof(Spinor), MPI_BYTE, blk.pid, blk.tag_send, MPI_COMM_WORLD, &status[num_status++]);
		MPI_Irecv(blk.buf_recv, blk.size*sizeof(Spinor), MPI_BYTE, blk.pid, blk.tag_recv, MPI_COMM_WORLD, &status[num_status++]);
	}

	MPI_Waitall(num_status, status, MPI_STATUSES_IGNORE);
	sync_recv_buffer(gd, sf, 1);
	timing_end(tm_mpi);
}

void mp_transfer_spinor_odd(SpinorField &spinor)
{
	static MPI_Request *status;
	static int init = 0;
	static gd_t gd;
	int num_status;
	Spinor *sf = spinor.ptr;

	if(init == 0)
	{
		gd = gd_list[3];
		for(auto& blk : gd)
		{
			blk.buf_send = new Spinor[blk.size];
			blk.buf_recv = new Spinor[blk.size];
		}
		status = new MPI_Request[2*gd.size()];
		init = 1;
	}

	num_status = 0;
	timing_start(tm_mpi);
	sync_send_buffer(gd, sf, 1);

	for(auto& blk : gd)
	{
		MPI_Isend(blk.buf_send, blk.size*sizeof(Spinor), MPI_BYTE, blk.pid, blk.tag_send, MPI_COMM_WORLD, &status[num_status++]);
		MPI_Irecv(blk.buf_recv, blk.size*sizeof(Spinor), MPI_BYTE, blk.pid, blk.tag_recv, MPI_COMM_WORLD, &status[num_status++]);
	}

	MPI_Waitall(num_status, status, MPI_STATUSES_IGNORE);
	sync_recv_buffer(gd, sf, 1);
	timing_end(tm_mpi);
}

void mp_transfer_spinor(SpinorField &spinor)
{
	static MPI_Request *status;
	static int init = 0;
	static gd_t gd;
	int num_status;
	Spinor *sf = spinor.ptr;

	if(init == 0)
	{
		gd = gd_list[1];
		for(auto& blk : gd)
		{
			blk.buf_send = new Spinor[blk.size];
			blk.buf_recv = new Spinor[blk.size];
		}
		status = new MPI_Request[2*gd.size()];
		init = 1;
	}

	num_status = 0;
	timing_start(tm_mpi);
	sync_send_buffer(gd, sf, 1);

	for(auto& blk : gd)
	{
		MPI_Isend(blk.buf_send, blk.size*sizeof(Spinor), MPI_BYTE, blk.pid, blk.tag_send, MPI_COMM_WORLD, &status[num_status++]);
		MPI_Irecv(blk.buf_recv, blk.size*sizeof(Spinor), MPI_BYTE, blk.pid, blk.tag_recv, MPI_COMM_WORLD, &status[num_status++]);
	}

	MPI_Waitall(num_status, status, MPI_STATUSES_IGNORE);
	sync_recv_buffer(gd, sf, 1);
	timing_end(tm_mpi);
}

void mp_transfer_links()
{
	static MPI_Request *status;
	static int init = 0;
	static gd_t gd;
	int num_status;
	suNg *gf = gfield_fund;

	if(init == 0)
	{
		gd = gd_list[0];
		for(auto& blk : gd)
		{
			blk.buf_send = new suNg[4*blk.size];
			blk.buf_recv = new suNg[4*blk.size];
		}
		status = new MPI_Request[2*gd.size()];
		init = 1;
	}

	num_status = 0;
	timing_start(tm_mpi);
	sync_send_buffer(gd, gf, 4);

	for(auto& blk : gd)
	{
		MPI_Isend(blk.buf_send, 4*blk.size*sizeof(suNg), MPI_BYTE, blk.pid, blk.tag_send, MPI_COMM_WORLD, &status[num_status++]);
		MPI_Irecv(blk.buf_recv, 4*blk.size*sizeof(suNg), MPI_BYTE, blk.pid, blk.tag_recv, MPI_COMM_WORLD, &status[num_status++]);
	}

	MPI_Waitall(num_status, status, MPI_STATUSES_IGNORE);
	sync_recv_buffer(gd, gf, 4);
	timing_end(tm_mpi);
}

void mp_transfer_clover_force()
{
	static MPI_Request *status;
	static int init = 0;
	static gd_t gd;
	int num_status;
	suNf *gf = cl_force;

	if(init == 0)
	{
		gd = gd_list[4];
		for(auto& blk : gd)
		{
			blk.buf_send = new suNf[6*blk.size];
			blk.buf_recv = new suNf[6*blk.size];
		}
		status = new MPI_Request[2*gd.size()];
		init = 1;
	}

	num_status = 0;
	timing_start(tm_mpi);
	sync_send_buffer(gd, gf, 6);

	for(auto& blk : gd)
	{
		MPI_Isend(blk.buf_send, 6*blk.size*sizeof(suNf), MPI_BYTE, blk.pid, blk.tag_send, MPI_COMM_WORLD, &status[num_status++]);
		MPI_Irecv(blk.buf_recv, 6*blk.size*sizeof(suNf), MPI_BYTE, blk.pid, blk.tag_recv, MPI_COMM_WORLD, &status[num_status++]);
	}

	MPI_Waitall(num_status, status, MPI_STATUSES_IGNORE);
	sync_recv_buffer(gd, gf, 6);
	timing_end(tm_mpi);
}

void mp_global_sum(void *buf, int count)
{
	MPI_Allreduce(MPI_IN_PLACE, buf, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void mp_global_max(void *buf, int count)
{
	MPI_Allreduce(MPI_IN_PLACE, buf, count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

void mp_broadcast(int *val)
{
	MPI_Bcast(val, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void mp_barrier()
{
	MPI_Barrier(MPI_COMM_WORLD);
}

void mp_gather(void *sbuf, void *rbuf, int sz)
{
	MPI_Gather(sbuf, sz, MPI_BYTE, rbuf, sz, MPI_BYTE, 0, MPI_COMM_WORLD);
}

void mp_scatter(void *sbuf, void *rbuf, int sz)
{
	MPI_Scatter(sbuf, sz, MPI_BYTE, rbuf, sz, MPI_BYTE, 0, MPI_COMM_WORLD);
}

#endif
