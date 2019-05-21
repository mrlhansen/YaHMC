#include <complex.h>
#include <logger.h>

void lu(int N, complex *A, complex *x, complex *b)
{
	double big;
	int row, itmp;
	complex ctmp;
	int mutate[N];

	// Setup mutate
	for(int i = 0; i < N; i++)
	{
		mutate[i] = i;
	}

	// LU factorization
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			for(int k = 0; k < j; k++)
			{
				A[j*N+i] -= A[j*N+k] * A[k*N+i];
			}
		}

		big = cabs(A[i*N+i]);
		row = i;

		for(int j = i+1; j < N; j++)
		{
			for(int k = 0; k < i; k++)
			{
				A[j*N+i] -= A[j*N+k] * A[k*N+i];
			}

			if(cabs(A[j*N+i]) > big)
			{
				big = cabs(A[j*N+i]);
				row = j;
			}
		}

		if(big < 1.0e-14)
		{
			lprintf("LU", CRITICAL, "LU decomposition failed: matrix is singular");
		}

		if(row != i)
		{
			for(int k = 0; k < N; k++)
			{
				ctmp = A[row*N+k];
				A[row*N+k] = A[i*N+k];
				A[i*N+k] = ctmp;
			}

			itmp = mutate[row];
			mutate[row] = mutate[i];
			mutate[i] = itmp;
		}

		for(int k = i+1; k < N; k++)
		{
			A[k*N+i] /= A[i*N+i];
		}
	}

	// Should we continue?
	if(x == 0 || b == 0)
	{
		return;
	}

	// Forward substitution
	for(int i = 0; i < N; i++)
	{
		x[i] = b[mutate[i]];
		for(int k = 0; k < i; k++)
		{
			x[i] -= A[i*N+k] * x[k];
		}
	}

	// Backward substitution
	for(int i = N-1; i >= 0; i--)
	{
		for(int k = i+1; k < N; k++)
		{
			x[i] -= A[i*N+k] * x[k];
		}
		x[i] /= A[i*N+i];
	}
}

void ldl(int N, complex *A, complex *x, complex *b)
{
	// LDL factorization
	for(int i = 0; i < N; i++)
	{
		for(int k = 0; k < i; k++)
		{
			A[i*N+i].re -= cabs2(A[i*N+k]) * A[k*N+k].re;
		}

		for(int j = i+1; j < N; j++)
		{
			for(int k = 0; k < i; k++)
			{
				A[j*N+i] -= A[j*N+k] * conj(A[i*N+k]) * A[k*N+k].re;
			}
			A[j*N+i] /= A[i*N+i].re;
		}
	}

	// Should we continue?
	if(x == 0 || b == 0)
	{
		return;
	}

	// Forward substitution
	for(int i = 0; i < N; i++)
	{
		x[i] = b[i];
		for(int k = 0; k < i; k++)
		{
			x[i] -= A[i*N+k] * x[k];
		}
	}

	// Backward substitution
	for(int i = N-1; i >= 0; i--)
	{
		x[i] /= A[i*N+i].re;
		for(int k = i+1; k < N; k++)
		{
			x[i] -= conj(A[k*N+i]) * x[k];
		}
	}
}
