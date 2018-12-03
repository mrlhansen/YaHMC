#include <bcs.h>
#include <global.h>

static vector<int> listA;
static vector<int> listB;
static vector<int> listC;
double *plaq_weights;
double *rect_weights;
bc_par bcs;

void bcs_init(double cf)
{
#ifdef BC_PERIODIC
	lprintf("BCS", INFO, "Boundary conditions: periodic");
#endif

#ifdef BC_ANTIPERIODIC
	extended_sites_for(id)
	{
		if(global_time(id) == 0)
		{
			listA.push_back(id);
		}
	}

	lprintf("BCS", INFO, "Boundary conditions: antiperiodic");
#endif

#ifdef BC_OPEN
	plaq_weights = new double[16*outer.vol4];
	rect_weights = plaq_weights;
	volume = (global.dim_t-1)*global.vol3;

	for(int id = 0; id < 16*outer.vol4; id++)
	{
		plaq_weights[id] = 1.0;
	}

	extended_sites_for(id)
	{
		if(global_time(id) == 0)
		{
			listB.push_back(id);
			plaq_weight(id,1,2) = 0.5;
			plaq_weight(id,1,3) = 0.5;
			plaq_weight(id,2,3) = 0.5;
			plaq_weight(id,2,1) = 0.5;
			plaq_weight(id,3,1) = 0.5;
			plaq_weight(id,3,2) = 0.5;
		}

		if(global_time(id) == global.dim_t-1)
		{
			listA.push_back(id);
			listB.push_back(id);
			plaq_weight(id,1,2) = 0.5;
			plaq_weight(id,1,3) = 0.5;
			plaq_weight(id,2,3) = 0.5;
			plaq_weight(id,2,1) = 0.5;
			plaq_weight(id,3,1) = 0.5;
			plaq_weight(id,3,2) = 0.5;
		}
	}

	bcs.cf = cf;
	bcs.Ta = 1;
	bcs.Tb = global.dim_t - 2;

	lprintf("BCS", INFO, "Boundary conditions: open");
	lprintf("BCS", INFO, "Coefficients: cf = %1.6f", cf);
#endif

#ifdef BC_SF
	rect_weights = new double[16*outer.vol4];
	plaq_weights = new double[16*outer.vol4];
	volume = (global.dim_t-2)*global.vol3;

	for(int id = 0; id < 16*outer.vol4; id++)
	{
		rect_weights[id] = 1.0;
		plaq_weights[id] = 1.0;
	}

	extended_sites_for(id)
	{
		if(global_time(id) == 0)
		{
			listA.push_back(id);
			rect_weight(id,1,2) = 0.0;
			rect_weight(id,1,3) = 0.0;
			rect_weight(id,2,3) = 0.0;
			rect_weight(id,0,1) = 1.5;
			rect_weight(id,0,2) = 1.5;
			rect_weight(id,0,3) = 1.5;
			plaq_weight(id,1,2) = 0.5;
			plaq_weight(id,1,3) = 0.5;
			plaq_weight(id,2,3) = 0.5;
			plaq_weight(id,2,1) = 0.5;
			plaq_weight(id,3,1) = 0.5;
			plaq_weight(id,3,2) = 0.5;
		}

		if(global_time(id) == global.dim_t-3)
		{
			rect_weight(id,0,1) = 1.5;
			rect_weight(id,0,2) = 1.5;
			rect_weight(id,0,3) = 1.5;
		}

		if(global_time(id) == global.dim_t-2)
		{
			listB.push_back(id);
			rect_weight(id,1,2) = 0.0;
			rect_weight(id,1,3) = 0.0;
			rect_weight(id,2,3) = 0.0;
			plaq_weight(id,1,2) = 0.5;
			plaq_weight(id,1,3) = 0.5;
			plaq_weight(id,2,3) = 0.5;
			plaq_weight(id,2,1) = 0.5;
			plaq_weight(id,3,1) = 0.5;
			plaq_weight(id,3,2) = 0.5;
		}

		if(global_time(id) == global.dim_t-1)
		{
			listC.push_back(id);
		}
	}

	bcs.cf = cf;
	bcs.Ta = 1;
	bcs.Tb = global.dim_t - 3;

	lprintf("BCS", INFO, "Boundary conditions: sf");
	lprintf("BCS", INFO, "Coefficients: cf = %1.6f", cf);
#endif
}

void apply_bcs_on_spinor_field(SpinorField &sptr)
{
#ifdef BC_OPEN
	Spinor u;
	for(int &id : listB)
	{
		sptr[id] = u;
	}
#endif

#ifdef BC_SF
	Spinor u;
	for(int &id : listA)
	{
		sptr[id] = u;
	}
	for(int &id : listB)
	{
		sptr[id] = u;
	}
	for(int &id : listC)
	{
		sptr[id] = u;
	}
#endif
}

void apply_bcs_on_momentum_field()
{
#ifdef BC_OPEN
	suNa u;
	for(int &id : listA)
	{
		momentum(id,0) = u;
	}
#endif

#ifdef BC_SF
	suNa u;
	for(int &id : listA)
	{
		momentum(id,1) = u;
		momentum(id,2) = u;
		momentum(id,3) = u;
	}
	for(int &id : listB)
	{
		momentum(id,0) = u;
		momentum(id,1) = u;
		momentum(id,2) = u;
		momentum(id,3) = u;
	}
	for(int &id : listC)
	{
		momentum(id,0) = u;
		momentum(id,1) = u;
		momentum(id,2) = u;
		momentum(id,3) = u;
	}
#endif
}

void apply_bcs_on_gauge_field()
{
#ifdef BC_OPEN
	suNg u;
	for(int &id : listA)
	{
		link(id,0) = u;
	}
#endif

#ifdef BC_SF
	suNg c1, c2;
	suNg C1, C2, Z;
	double L;

#if (NC==2)

	double v1[NC] = {-M_PI, M_PI};
	double v2[NC] = {-3.*M_PI, 3.*M_PI};
	L = 4 * global.dim_x;

#elif (NC==3)

	double v1[NC] = {-M_PI, 0, M_PI};
	double v2[NC] = {-5.*M_PI, 2.*M_PI, 3.*M_PI};
	L = 6 * global.dim_x;

#endif

	for(int i = 0; i < NC; i++)
	{
		c1.im[i*NC+i] = v1[i]/L;
		c2.im[i*NC+i] = v2[i]/L;
	}

	C1 = exp(c1);
	C2 = exp(c2);

	for(int &id : listA)
	{
		link(id,1) = C1;
		link(id,2) = C1;
		link(id,3) = C1;
	}

	for(int &id : listB)
	{
		link(id,0) = Z;
		link(id,1) = C2;
		link(id,2) = C2;
		link(id,3) = C2;
	}

	for(int &id : listC)
	{
		link(id,0) = Z;
		link(id,1) = Z;
		link(id,2) = Z;
		link(id,3) = Z;
	}
#endif
}

void apply_bcs_on_represented_gauge_field()
{
#ifdef BC_ANTIPERIODIC
	for(int &id : listA)
	{
		fermion_link(id,0) *= -1;
	}
#endif
}
