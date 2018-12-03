#include <global.h>
#include <repr.h>
#include <logger.h>

static int check_homomorphism()
{
	suNg u, v, uv;
	suNf U, V, UV;
	double r1, r2;
	int ret = 0;

	// h(u^dagger) = h(u)^dagger
	u.random();
	v = dagger(u);

	represent_links(&U, &u, 1);
	represent_links(&V, &v, 1);

	U = dagger(U) - V;
	r1 = U.sqnorm();
	r1 = sqrt(r1);

	// Evaluate result
	if(r1 < 1e-14)
	{
		lprintf("CHECK", INFO, "Group homomorphism (part 1): PASSED [diff = %1.6e]", r1);
	}
	else
	{
		lprintf("CHECK", INFO, "Group homomorphism (part 1): FAILED [diff = %1.6e]", r1);
		ret += 1;
	}

	// h(u*v) = h(u)*h(v)
	u.random();
	v.random();
	uv = u*v;

	represent_links(&UV, &uv, 1);
	represent_links(&U, &u, 1);
	represent_links(&V, &v, 1);

	UV -= U*V;
	r2 = UV.sqnorm();
	r2 = sqrt(r2);

	// Evaluate result
	if(r2 < 1e-14)
	{
		lprintf("CHECK", INFO, "Group homomorphism (part 2): PASSED [diff = %1.6e]", r2);
	}
	else
	{
		lprintf("CHECK", INFO, "Group homomorphism (part 2): FAILED [diff = %1.6e]", r2);
		ret += 1;
	}

	// Return result
	return ret;
}

static int check_orthogonality()
{
	double r1 = 0;
	double r2 = 0;

	// Test generators
	for(int i = 0; i < NG; i++)
	{
		for(int j = 0; j < NG; j++)
		{
			r1 += trace(iTfund[i], iTfund[j]);
			r2 += trace(iTrepr[i], iTrepr[j]);
			if(i == j)
			{
				r1 += TF;
				r2 += REPR_TRACE;
			}
		}
	}

	// Average result
	r1 /= NG*NG;
	r2 /= NG*NG;
	r1 = fabs(r1 + r2) / 2.0;

	// Evaluate result
	if(r1 < 1e-14)
	{
		lprintf("CHECK", INFO, "Generator orthogonality: PASSED [diff = %1.6e]", r1);
		return 0;
	}
	else
	{
		lprintf("CHECK", INFO, "Generator orthogonality: FAILED [diff = %1.6e]", r1);
		return 1;
	}
}

int check_group()
{
	int ret = 0;
	ret += check_homomorphism();
	ret += check_orthogonality();
	return ret;
}
