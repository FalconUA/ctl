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
	void fac_pollard(T const& n, T& r1, T& r2, TypeSpecifications<T> const& algo)
	{
		__TRIVIAL_FACTORIZATION__
		
	    T t, d;
		t = algo.sqrt(n);

	    for(T i = 1, a = 2, b = 2; i <= t; i = i + 1){
	        a = (a*a + 1)%n, b = (b*b + 1)%n, b = (b*b + 1)%n;
	        d = algo.gcd(a-b,n);
	        if((d > 1)&&(d < n)){
				r1 = d;
				r2 = n / d;
	            return ;
	        }
		}
		r1 = n;
		r2 = 1;
	}

	#undef __TRIVIAL_FACTORIZATION__

}
