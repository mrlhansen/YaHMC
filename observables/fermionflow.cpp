#include <wilsonflow.h>
#include <linalg.h>
#include <global.h>

static SpinorField *phi0;
static SpinorField *phi1;
static SpinorField *phi2;
static SpinorField dphi;
static int nsf;

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
	for(int i = 0; i < nsf; i++)
	{
		wf_laplacian(dphi, phi0[i]);
		spinor_copy(phi1[i], phi0[i]);
		spinor_mulr_add_assign(phi1[i], epsilon/4.0, dphi);
	}

	// Gauge 0
	wf_zeta(epsilon/4.0);
	wf_update(-17.0/9.0);

	// Fermion 1
	represent_gauge_field();
	for(int i = 0; i < nsf; i++)
	{
		wf_laplacian(dphi, phi1[i]);
		spinor_mulr(phi2[i], -8.0/9.0, phi1[i]);
		spinor_mulr_add_assign(phi2[i], 17.0/9.0, phi0[i]);
		spinor_mulr_add_assign(phi2[i], 8.0*epsilon/9.0, dphi);
	}

	// Gauge 1
	wf_zeta(8.0*epsilon/9.0);
	wf_update(-1.0);

	// Fermion 2
	represent_gauge_field();
	for(int i = 0; i < nsf; i++)
	{
		wf_laplacian(dphi, phi2[i]);
		spinor_copy(phi0[i], phi1[i]);
		spinor_mulr_add_assign(phi0[i], 3.0*epsilon/4.0, dphi);
	}

	// Gauge 2
	wf_zeta(3.0*epsilon/4.0);
	wf_update(0);
}

void wf_fermion_init(SpinorField *ptr, int count)
{
	// Number of fields to be evolved
	nsf = count;

	// Allocate memory
	phi0 = ptr;
	phi1 = new SpinorField[nsf];
	phi2 = new SpinorField[nsf];

	spinor_allocate(dphi);
	for(int i = 0; i < nsf; i++)
	{
		spinor_allocate(phi1[i]);
		spinor_allocate(phi2[i]);
	}

	// Initialize gauge part
	wf_init();
}

void wf_fermion_free()
{
	delete[] phi1;
	delete[] phi2;
}
