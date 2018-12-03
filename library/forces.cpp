#include <forces.h>
#include <global.h>
#include <linalg.h>
#include <dirac.h>
#include <clover.h>
#include <staples.h>
#include <bcs.h>
#include <timing.h>
#include <cmath>

// Creates an suNf force matrix
// fmat = a_lhs # a_rhs + b_lhs # b_rhs
static suNf fmat_create(Spinor &a_lhs, Spinor &a_rhs, Spinor &b_lhs, Spinor &b_rhs)
{
	suNf fmat;
	for(int i = 0; i < REPR_DIM; i++)
	{
		for(int j = 0; j < REPR_DIM; j++)
		{
			for(int k = 0; k < 4; k++)
			{
				fmat.re[i*REPR_DIM+j] += a_lhs.re[i][k]*a_rhs.re[j][k] + a_lhs.im[i][k]*a_rhs.im[j][k];
				fmat.im[i*REPR_DIM+j] += a_lhs.im[i][k]*a_rhs.re[j][k] - a_lhs.re[i][k]*a_rhs.im[j][k];
				fmat.re[i*REPR_DIM+j] += b_lhs.re[i][k]*b_rhs.re[j][k] + b_lhs.im[i][k]*b_rhs.im[j][k];
				fmat.im[i*REPR_DIM+j] += b_lhs.im[i][k]*b_rhs.re[j][k] - b_lhs.re[i][k]*b_rhs.im[j][k];
			}
		}
	}
	return fmat;
}

// Core function for the clover force
static void clover_force_core(double dt)
{
	double coeff = (get_csw() * dt) / 8.0;
	suNf fmat, W[9];
	int num, sign;

	// Communicate force
	mp_transfer_clover_force();

	// Timing
	timing_start(tm_cforce);

	// Perform momentum update
	#pragma omp parallel for private(fmat,num,sign,W)
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			for(int nu = 0; nu < 4; nu++)
			{
				if(mu == nu) continue;

				// Coordinates
				int o1 = fw_index(id,mu); // x + mu
				int o2 = fw_index(id,nu); // x + nu
				int o3 = bk_index(id,nu); // x - nu
				int o4 = fw_index(o3,mu); // x + mu - nu
				int o5 = fw_index(o2,mu); // x + mu + nu

				if(mu < nu)
				{
					num = nu*(nu-1)/2+mu;
					sign = +1;
				}
				else
				{
					num = mu*(mu-1)/2+nu;
					sign = -1;
				}

				// Force matrices
				suNf &Z0 = clover_force(id,num);
				suNf &Z1 = clover_force(o1,num);
				suNf &Z2 = clover_force(o3,num);
				suNf &Z3 = clover_force(o4,num);
				suNf &Z4 = clover_force(o5,num);
				suNf &Z5 = clover_force(o2,num);

				// Construct links
				W[0] = dagger(fermion_link(o3,mu));
				W[1] = fermion_link(o3,nu);
				W[2] = fermion_link(o1,nu);
				W[3] = dagger(fermion_link(o2,mu));
				W[4] = dagger(fermion_link(id,nu));
				W[5] = dagger(fermion_link(o4,nu));
				W[6] = W[0]*W[1];
				W[7] = W[2]*W[3];
				W[8] = W[7]*W[4]-W[5]*W[6];

				// Calculate sum of forces
				fmat  = W[8]*Z0+Z1*W[8];
				fmat -= W[5]*(W[0]*Z2*W[1]+Z3*W[6]);
				fmat += (W[2]*Z4*W[3]+W[7]*Z5)*W[4];
				fmat  = fermion_link(id,mu)*fmat;

				// Project on force matrix
				momentum(id,mu) += algebra_project(fmat, sign*coeff);
			}
		}
	}

	// Timing
	timing_end(tm_cforce);
}

// Calculates the force: log(det(D_oo))
void clover_logdet_force(double mass, double residue)
{
	// Timing
	timing_start(tm_cforce);

	// Compute force matrices
	compute_clover_force(mass, residue);

	// Timing
	timing_end(tm_cforce);
}

// Calculates the force: lhs^\dagger*\dot{Q}_xx*rhs
void clover_fermion_force(double residue, SpinorField &lhs, SpinorField &rhs)
{
	Spinor tmp_lhs, tmp_rhs;

	// Timing
	timing_start(tm_cforce);

	// Construct force matrices
	#pragma omp parallel for private(tmp_lhs,tmp_rhs)
	sites_for(id)
	{
		tmp_rhs = rhs[id];
		tmp_rhs.g5_sigma(0,1);
		tmp_lhs = lhs[id];
		tmp_lhs.g5_sigma(0,1);
		clover_force(id,0) += residue * fmat_create(tmp_lhs, rhs[id], tmp_rhs, lhs[id]);

		tmp_rhs = rhs[id];
		tmp_rhs.g5_sigma(0,2);
		tmp_lhs = lhs[id];
		tmp_lhs.g5_sigma(0,2);
		clover_force(id,1) += residue * fmat_create(tmp_lhs, rhs[id], tmp_rhs, lhs[id]);

		tmp_rhs = rhs[id];
		tmp_rhs.g5_sigma(1,2);
		tmp_lhs = lhs[id];
		tmp_lhs.g5_sigma(1,2);
		clover_force(id,2) += residue * fmat_create(tmp_lhs, rhs[id], tmp_rhs, lhs[id]);

		tmp_rhs = rhs[id];
		tmp_rhs.g5_sigma(0,3);
		tmp_lhs = lhs[id];
		tmp_lhs.g5_sigma(0,3);
		clover_force(id,3) += residue * fmat_create(tmp_lhs, rhs[id], tmp_rhs, lhs[id]);

		tmp_rhs = rhs[id];
		tmp_rhs.g5_sigma(1,3);
		tmp_lhs = lhs[id];
		tmp_lhs.g5_sigma(1,3);
		clover_force(id,4) += residue * fmat_create(tmp_lhs, rhs[id], tmp_rhs, lhs[id]);

		tmp_rhs = rhs[id];
		tmp_rhs.g5_sigma(2,3);
		tmp_lhs = lhs[id];
		tmp_lhs.g5_sigma(2,3);
		clover_force(id,5) += residue * fmat_create(tmp_lhs, rhs[id], tmp_rhs, lhs[id]);
	}

	// Timing
	timing_end(tm_cforce);
}

// Calculates the force: lhs^\dagger*\dot{Q}*rhs
void wilson_fermion_force(double dt, double mass, double residue, SpinorField &lhs, SpinorField &rhs)
{
	double coeff = dt * residue;
	Spinor tmp_lhs, tmp_rhs;
	suNf fmat;
	int fw;

#ifdef UPDATE_EO
	Dphi_oe(rhs, rhs);
	Dphi_oe(lhs, lhs);
#ifdef CLOVER_TERM
	Dphi_oo_inv(mass, rhs, rhs);
	Dphi_oo_inv(mass, lhs, lhs);
#endif
	coeff = -coeff;
#endif

#ifdef CLOVER_TERM
	clover_fermion_force(residue, lhs, rhs);
#endif

	// Communicate fields
	mp_transfer_spinor(lhs);
	mp_transfer_spinor(rhs);

	// Timing
	timing_start(tm_wforce);

	// Perform momentum update
	#pragma omp parallel for private(fmat,tmp_lhs,tmp_rhs,fw)
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			fw = fw_index(id,mu);

			// Construct force matrix
			tmp_rhs = rhs[fw];
			tmp_rhs.g5_one_minus_gamma(mu);
			tmp_lhs = lhs[fw];
			tmp_lhs.g5_one_minus_gamma(mu);

			fmat = fmat_create(tmp_lhs, rhs[id], tmp_rhs, lhs[id]);
			fmat = fermion_link(id,mu)*fmat;

			// Project on force matrix
			momentum(id,mu) += algebra_project(fmat, coeff);
		}
	}

	// Timing
	timing_end(tm_wforce);
}

// Calculates the force: rhs^\dagger*Q*\dot{Q}*rhs
void wilson_fermion_force_q(double dt, double mass, double residue, SpinorField &rhs)
{
	static SpinorField tmp;
	spinor_allocate(tmp);
	Hphi(mass, tmp, rhs);
	wilson_fermion_force(dt, mass, residue, tmp, rhs);
}

// Start fermion force evaluation
void begin_fermion_force()
{
#ifdef CLOVER_TERM
	suNf zero;
	sites_for(id)
	{
		for(int mu = 0; mu < 6; mu++)
		{
			clover_force(id,mu) = zero;
		}
	}
#endif
}

// Finalize fermion force evaluation
void end_fermion_force(double dt)
{
#ifdef CLOVER_TERM
	clover_force_core(dt);
#endif
}

// Calculates the force of the Wilson plaquette action
void wilson_gauge_force(double dt, double beta)
{
	double coeff = dt * (beta / NC);
	suNg force;

	// Timing
	timing_start(tm_gforce);

	// Perform momentum update
	#pragma omp parallel for private(force)
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			force = link(id,mu)*staples_wilson(id, mu);
			momentum(id,mu) -= algebra_project(force, coeff);
		}
	}

	// Timing
	timing_end(tm_gforce);
}

// Calculates the force of the Symanzik gauge action
void improved_gauge_force(double dt, double beta, double c0, double c1)
{
	double coeff = dt * (beta / NC);
	suNg part1, part2, force;

	// Timing
	timing_start(tm_gforce);

	// Perform momentum update
	#pragma omp parallel for private(part1,part2,force)
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			part1 = link(id,mu)*staples_wilson(id, mu);
			part2 = link(id,mu)*staples_symanzik(id, mu);
			force = c0*part1 + c1*part2;
			momentum(id,mu) -= algebra_project(force, coeff);
		}
	}

	// Timing
	timing_end(tm_gforce);
}
