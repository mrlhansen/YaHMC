#include <wilsonloop.h>
#include <global.h>
#include <bcs.h>

suNf clover_repr(int id, int mu, int nu)
{
	suNg gtmp[4];
	suNf ftmp[4], ret;
	int a, b, c;

	a = fw_index(id,mu);
	b = fw_index(id,nu);
	gtmp[0] = link(id,mu) * link(a,nu) * dagger(link(id,nu) * link(b,mu));

	a = bk_index(id,mu);
	b = fw_index(a,nu);
	gtmp[1] = link(id,nu) * dagger(link(a,nu) * link(b,mu)) * link(a,mu);

	a = bk_index(id,mu);
	b = bk_index(id,nu);
	c = bk_index(a,nu);
	gtmp[2] = dagger(link(c,nu) * link(a,mu)) * link(c,mu) * link(b,nu);

	a = bk_index(id,nu);
	b = fw_index(a,mu);
	gtmp[3] = dagger(link(a,nu)) * link(a,mu) * link(b,nu) * dagger(link(id,mu));

	represent_links(ftmp, gtmp, 4);
	ret = ftmp[0] + ftmp[1] + ftmp[2] + ftmp[3];
	return ret;
}

suNg clover_fund(int id, int mu, int nu)
{
	suNg ret;
	int a, b, c;

	a = fw_index(id,mu);
	b = fw_index(id,nu);
	ret = link(id,mu) * link(a,nu) * dagger(link(id,nu) * link(b,mu));

	a = bk_index(id,mu);
	b = fw_index(a,nu);
	ret += link(id,nu) * dagger(link(a,nu) * link(b,mu)) * link(a,mu);

	a = bk_index(id,mu);
	b = bk_index(id,nu);
	c = bk_index(a,nu);
	ret += dagger(link(c,nu) * link(a,mu)) * link(c,mu) * link(b,nu);

	a = bk_index(id,nu);
	b = fw_index(a,mu);
	ret += dagger(link(a,nu)) * link(a,mu) * link(b,nu) * dagger(link(id,mu));

	return ret;
}

double wilson_loop(int id, int mu, int nu, int height, int width)
{
	suNg one = link(id,mu);
	suNg two = link(id,nu);
	int pos = id;
	double tr;

	for(int ndx = 1; ndx < height; ndx++)
	{
		pos = fw_index(pos,mu);
		one *= link(pos,mu);
	}

	pos = fw_index(pos,mu);
	one *= link(pos,nu);

	for(int ndx = 1; ndx < width; ndx++)
	{
		pos = fw_index(pos,nu);
		one *= link(pos,nu);
	}

	pos = id;

	for(int ndx = 1; ndx < width; ndx++)
	{
		pos = fw_index(pos,nu);
		two *= link(pos,nu);
	}

	pos = fw_index(pos,nu);
	two *= link(pos,mu);

	for(int ndx = 1; ndx < height; ndx++)
	{
		pos = fw_index(pos,mu);
		two *= link(pos,mu);
	}

	two = dagger(two);
	tr = trace(one, two);

	if(height == 1 && width == 1)
	{
		tr *= plaq_weight(id,mu,nu);
	}
	else if(height == 1 && width == 2)
	{
		tr *= rect_weight(id,mu,nu);
	}
	else if(height == 2 && width == 1)
	{
		tr *= rect_weight(id,nu,mu);
	}

	return tr;
}

double avr_wilson(int height, int width)
{
	double sum = 0;

	#pragma omp parallel for reduction(+:sum)
	sites_for(id)
	{
		sum += wilson_loop(id, 1, 0, height, width);
		sum += wilson_loop(id, 2, 0, height, width);
		sum += wilson_loop(id, 2, 1, height, width);
		sum += wilson_loop(id, 3, 0, height, width);
		sum += wilson_loop(id, 3, 1, height, width);
		sum += wilson_loop(id, 3, 2, height, width);
	}

	mp_global_sum(&sum, 1);
	return sum/(NC*6*volume);
}
