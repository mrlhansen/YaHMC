#include <wilsonflow.h>
#include <linalg.h>
#include <global.h>

static SpinorField *phiA;
static SpinorField *phiB;
static SpinorField *phiC;
static SpinorField *delta_phiA;
static SpinorField *delta_phiB;
static SpinorField *delta_phiC;
static int num_sf;

static void wf_laplacian(SpinorField &dptr, SpinorField &sptr)
{
	mp_transfer_spinor(sptr);

	#pragma omp parallel for
	sites_for(id)
	{
		dptr[id] = -8.0*sptr[id];
		for(int mu = 0; mu < 4; mu++)
		{
			int fw = fw_index(id,mu);
			int bk = bk_index(id,mu);
			dptr[id] += fermion_link(id,mu)*sptr[fw];
			dptr[id] += dagger(fermion_link(bk,mu))*sptr[bk];
		}
	}
}

void wf_fermion_iterate(double epsilon)
{
	// Fermion 0
	represent_gauge_field();
	for(int i = 0; i < num_sf; i++)
	{
		wf_laplacian(delta_phiA[i], phiA[i]);
		spinor_copy(phiB[i], phiA[i]);
		spinor_mulr_add_assign(phiB[i], epsilon/4.0, delta_phiA[i]);
	}

	// Gauge 0
	wf_zeta(epsilon/4.0);
	wf_update(-17.0/9.0);

	// Fermion 1
	represent_gauge_field();
	for(int i = 0; i < num_sf; i++)
	{
		wf_laplacian(delta_phiB[i], phiB[i]);
		spinor_copy(phiC[i], phiA[i]);
		spinor_mulr_sub_assign(phiC[i], 2.0*epsilon/9.0, delta_phiA[i]);
		spinor_mulr_add_assign(phiC[i], 8.0*epsilon/9.0, delta_phiB[i]);
	}

	// Gauge 1
	wf_zeta(8.0*epsilon/9.0);
	wf_update(-1.0);

	// Fermion 2
	represent_gauge_field();
	for(int i = 0; i < num_sf; i++)
	{
		wf_laplacian(delta_phiC[i], phiC[i]);
		spinor_copy(phiA[i], phiB[i]);
		spinor_mulr_add_assign(phiA[i], 3.0*epsilon/4.0, delta_phiC[i]);
	}

	// Gauge 2
	wf_zeta(3.0*epsilon/4.0);
	wf_update(0);
}

void wf_fermion_init(SpinorField *ptr, int nsf)
{
	// Number of fields to be evolved
	num_sf = nsf;

	// Allocate memory
	phiA = ptr;
	phiB = new SpinorField[nsf];
	phiC = new SpinorField[nsf];

	delta_phiA = new SpinorField[nsf];
	delta_phiB = new SpinorField[nsf];
	delta_phiC = new SpinorField[nsf];

	for(int i = 0; i < nsf; i++)
	{
		spinor_allocate(phiB[i]);
		spinor_allocate(phiC[i]);
		spinor_allocate(delta_phiA[i]);
		spinor_allocate(delta_phiB[i]);
		spinor_allocate(delta_phiC[i]);
	}

	// Initialize gauge part
	wf_init();
}

void wf_fermion_free()
{
	delete[] phiB;
	delete[] phiC;
	delete[] delta_phiA;
	delete[] delta_phiB;
	delete[] delta_phiC;
}
