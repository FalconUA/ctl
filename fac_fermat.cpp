#include "ctlfac.hpp"
#include "typespec.hpp"

namespace mp {

	#define __TRIVIAL_FACTORIZATION__ \
		T N = n;                      \
		if (N == 1){                  \
			r1 = 1, r2 = 1;           \
			return;                   \
		}                             \
		if (N % 2 == 0){              \
			N = N/2; r1 = 2; r2 = N;  \
			return ;                  \
		}                             \
		if (!(algo.isprime == NULL))  \
			if (algo.isprime(N)){     \
				r1 = N, r2 = 1;       \
				return;               \
			}

	template <typename T>
	void fac_fermat(T const& n, T& r1, T& r2, TypeSpecifications<T> const& algo)
	{
		__TRIVIAL_FACTORIZATION__
	
		auto isfullsqr = [&algo](T x) -> bool {
				T t = algo.sqrt(x);
				return (t*t == x);
			};
	
		T a = algo.sqrt(N);
		if (a*a < N) a = a + 1;
		T b = a*a - N;
		while (!isfullsqr(b)){
			b = b + 2*a + 1;	
			a = a + 1;
		}
		r1 = a - algo.sqrt(b);
		r2 = a + algo.sqrt(b);
		return ;
	}

	#undef __TRIVIAL_FACTORIZATION__

}
