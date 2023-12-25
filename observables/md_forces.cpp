#include <monomials.h>
#include <global.h>

void measure_md_forces(Monomial *mon)
{
	double nsq, avg, var;
	double max = 0;
	double sum = 0;
	double sumsq = 0;
	suNa zero;

	// Set momenta to zero
	#pragma omp parallel for
	sites_for(id)
	{
		momentum(id,0) = zero;
		momentum(id,1) = zero;
		momentum(id,2) = zero;
		momentum(id,3) = zero;
	}

	// Update momenta
	mon->reset();
	begin_fermion_force();
	mon->update(1);
	end_fermion_force(1);

	// Calculate sum of forces
	#pragma omp parallel for private(nsq) reduction(+:sum,sumsq) reduction(max: max)
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			nsq = momentum(id,mu).sqnorm(); // Multiply by TF?
			sum += nsq;
			sumsq += nsq*nsq;
			if(nsq > max) max = nsq;
		}
	}

	mp_global_sum(&sum, 1);
	mp_global_sum(&sumsq, 1);
	mp_global_max(&max, 1);

	// Print information
	avg = sum/(4*global.vol4);
	var = sumsq/(4*global.vol4) - avg*avg;
	lprintf("MDFORCE", INFO, "MD Forces |F|^2: sum = %1.8e, avg = %1.8e, max = %1.8e, var = %1.8e", sum, avg, max, var);
}
