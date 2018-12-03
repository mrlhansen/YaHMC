#include <staples.h>
#include <global.h>
#include <wilsonloop.h>
#include <bcs.h>
#include <cmath>

suNg staples_symanzik(int id, int mu)
{
	suNg staple;
	int o1, o2, o3, o4;

	for(int nu = 0; nu < 4; nu++)
	{
		if(nu == mu) continue;

		// *---*---*
		//         |
		// *---*---*
		o1 = fw_index(id,mu);
		o2 = fw_index(o1,nu);
		o3 = fw_index(id,nu);
		o4 = fw_index(o3,nu);
		staple += rect_weight(id,mu,nu) * link(o1,nu) * link(o2,nu) * dagger(link(id,nu) * link(o3,nu) * link(o4,mu));

		// *---*---*
		// |
		// *---*---*
		o1 = bk_index(id,nu);
		o2 = bk_index(o1,nu);
		o3 = fw_index(o2,mu);
		o4 = fw_index(o3,nu);
		staple += rect_weight(o2,mu,nu) * dagger(link(o2,mu) * link(o3,nu) * link(o4,nu)) * link(o2,nu) * link(o1,nu);

		// *---*
		// |   |
		// *   *
		// |
		// *---*
		o1 = bk_index(id,nu);
		o2 = fw_index(o1,mu);
		o3 = fw_index(o2,mu);
		o4 = fw_index(id,mu);
		staple += rect_weight(o1,nu,mu) * link(o4,mu) * dagger(link(o1,mu) * link(o2,mu) * link(o3,nu)) * link(o1,nu);

		// *---*
		// |   |
		// *   *
		//     |
		// *---*
		o1 = fw_index(id,mu);
		o2 = fw_index(o1,mu);
		o3 = fw_index(id,nu);
		o4 = fw_index(o3,mu);
		staple += rect_weight(id,nu,mu) * link(o1,mu) * link(o2,nu) * dagger(link(id,nu) * link(o3,mu) * link(o4,mu));

		// *---*
		// |
		// *   *
		// |   |
		// *---*
		o1 = bk_index(id,mu);
		o2 = bk_index(o1,nu);
		o3 = fw_index(o2,mu);
		o4 = fw_index(o3,mu);
		staple += rect_weight(o2,nu,mu) * dagger(link(o2,mu) * link(o3,mu) * link(o4,nu)) * link(o2,nu) * link(o1,mu);

		// *---*
		//     |
		// *   *
		// |   |
		// *---*
		o1 = bk_index(id,mu);
		o2 = fw_index(o1,nu);
		o3 = fw_index(o2,mu);
		o4 = fw_index(id,mu);
		staple += rect_weight(o1,nu,mu) * link(o4,nu) * dagger(link(o1,nu) * link(o2,mu) * link(o3,mu)) * link(o1,mu);
	}

	return staple;
}

suNg staples_wilson(int id, int mu)
{
	suNg staple;
	int o1, o2;

	for(int nu = 0; nu < 4; nu++)
	{
		if(nu == mu) continue;

		// *---*
		//     |
		// *---*
		o1 = fw_index(id,mu);
		o2 = fw_index(id,nu);
		staple += plaq_weight(id,mu,nu) * link(o1,nu) * dagger(link(id,nu) * link(o2,mu));

		// *---*
		// |
		// *---*
		o1 = bk_index(id,nu);
		o2 = fw_index(o1,mu);
		staple += plaq_weight(o1,nu,mu) * dagger(link(o1,mu) * link(o2,nu)) * link(o1,nu);
	}

	return staple;
}

void staples_test()
{
	suNg staple;
	double s = 0, p = 0;

	#pragma omp parallel for private(staple) reduction(+:s)
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			staple = staples_wilson(id, mu);
			s += trace(link(id,mu)*staple);
		}
	}

	mp_global_sum(&s, 1);
	s = s / (NC*4*6*volume); // NC for scaling, 4 directions, 6 staples per link
	p = avr_wilson(1,1);

	if(fabs(p-s) < 1.0e-12)
	{
		lprintf("STAPLES", INFO, "Staples test PASSED, difference: %1.4e", fabs(p-s));
	}
	else
	{
		lprintf("STAPLES", CRITICAL, "Staples test FAILED, difference: %1.4e", fabs(p-s));
	}
}
