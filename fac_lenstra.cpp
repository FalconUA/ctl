#include "typespec.hpp"
#include "ctlfac.hpp"

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
	class eliptic_curve{
	public:
	    T x, y, a, n;
	    eliptic_curve(): a(0), n(0), x(0), y(0) {};
	    eliptic_curve(T N, TypeSpecifications<T> const& algo){
        	T b, g;
	        n=N;
	        do {
	            x = algo.random(n);
	            y = algo.random(n);
	            a = algo.random(n);
            	b = (y*y-x*x*x-a*x) % n;
	            T temp=(4*a*a*a+27*b*b) % n;
	            g = algo.gcd(n,temp);
	        } while(g == n);
	    }
	};

	template <typename T>	
	void fac_lenstra(T const& n, T& r1, T& r2, TypeSpecifications<T> const& algo)
	{
		__TRIVIAL_FACTORIZATION__

	    T B1 = 10000;
		eliptic_curve<T> P;
	    while(true){
		    P = eliptic_curve<T>(n, algo);
			T b;
			b = (P.y*P.y - P.x*P.x*P.x - P.a*P.x) % n;
			T temp=(4*P.a*P.a*P.a + 27*b*b) % n;
			T g = algo.gcd(n,temp);
			if(g == n) continue;
	        if(g<n && g>1){
				r1 = g;
				r2 = n/g;
	            return;
	        }
	        break;
		}

	    T B2 = B1;

	    while(true){
		    for(T p = 2; p<B1; p=p+1){
	            if(algo.isprime(p)){
	                T r = 0, c_p = 1;
	                for(T t; c_p < B1; ){
	                    t = c_p*c_p;
	                    if(t > 1)
	                        if(t < B1)
	                            c_p = t, r = r*2;
	                        else
	                            while(c_p < B1)
									c_p = c_p * p, r = r+1;
	                    else
	                        c_p = c_p * p, r = r+1;
	                }
	                r = r - 1;
	                for(T j = 0; j < r; j = j+1)
	                    P.x = P.x * p, P.y = P.y * p;
		            T t = algo.rev(2*P.y,n);
	                T d = algo.gcd(n,t);
	                if(d>1){
	                    r1 = d;
	                    r2 = n/d;
	                    return;
	                }
	            }
	        }
	    }
	}


	#undef __TRIVIAL_FACTORIZATION__
}
