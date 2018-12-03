#ifndef MONOMIALS_H
#define MONOMIALS_H

#include <forces.h>
#include <rational.h>
#include <mre.h>

#include <vector>
using namespace std;

/* Monomial base class */
class Monomial {
	public:
		virtual void reset() = 0;
		virtual double action() = 0;
		virtual void update(double) = 0;
};

/* HMC monomial */
class MonomialHMC : public Monomial {
	private:
		double omass;
		double smass;
		double precision;
		double la;
		SpinorField pf;
		SpinorField atmp;
		SpinorField btmp;
		mre_par mpar1;
		mre_par mpar2;
	public:
		MonomialHMC(double, double, double, int);
		void reset();
		void update(double);
		double action();
};

/* RHMC monomial */
class MonomialRHMC : public Monomial {
	private:
		double mass;
		int n_pf, d_pf;
		int n_fc, d_fc;
		int order;
		double inv_precision;
		double rap_precision;
		double la;
		rational_app rpf;
		rational_app rfc;
		SpinorField pf;
		SpinorField tmp;
		SpinorField *out;
	public:
		MonomialRHMC(double, double, double);
		void reset();
		void update(double);
		double action();
};

/* Hasenbusch monomial */
class MonomialHB : public Monomial {
	private:
		double omass;
		double smass;
		double dm;
		double precision;
		double la;
		SpinorField pf;
		SpinorField atmp;
		SpinorField btmp;
		SpinorField ctmp;
	public:
		MonomialHB(double, double, double);
		void reset();
		void update(double);
		double action();
};

/* Gauge monomial */
class MonomialGauge : public Monomial {
	private:
		double beta;
		double c0;
		double c1;
		int improved;
	public:
		MonomialGauge(double, double);
		void reset();
		void update(double);
		double action();
};

/* MonomialList class */
typedef struct {
	Monomial *ptr;
	int level;
} MonomialItem;

class MonomialList {
	private:
		vector<MonomialItem> list;
	public:
		void add(Monomial*, int);
		void update(double, int);
		void reset();
		double action();
};

#endif
