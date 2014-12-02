#include "ctllog.hpp"
#include "typespec.hpp"

namespace mp {

	template <typename T> 
	T log_pollard(T const& a, T const& b, T const& m, TypeSpecifications<T> const& algo)
	{
	init: 
		T a0 = 0, b0 = 0, x0 = 1, i = 1, lower_bound = 0, r, x;
		auto f = [&a, &b, &m](T x) -> T {
				return (x % 3 == 0)? (x*x):(
					   (x % 3 == 1)? (x*a):(
					   (x % 3 == 2)? (x*b):(0))) % (m+1);
			};
		auto g = [&m](T x, T n) -> T {
				return (x % 3 == 0)? (2*n):(
					   (x % 3 == 1)? (n+1):(
					   (x % 3 == 2)? (n  ):(0))) % (m);
			};
		auto h = [&m](T x, T n) -> T {
				return (x % 3 == 0)? (2*n):(
					   (x % 3 == 1)? (n  ):(
					   (x % 3 == 2)? (n+1):(0))) % (m);			
			};

	step1:
		T xi, xi_1 = x0, x2i, x2i_2 = x0;
		T ai, ai_1 = a0, a2i, a2i_2 = a0;
		T bi, bi_1 = b0, b2i, b2i_2 = b0;
	step2: // calculating xi
		xi = f(xi_1);
		ai = g(xi_1, ai_1);
		bi = h(xi_1, bi_1);
	step3: // and calculate x2i at the same time
		x2i = f(f(x2i_2));
		a2i = g(f(x2i_2), g(x2i_2, a2i_2));	
		b2i = h(f(x2i_2), h(x2i_2, b2i_2));
	step4:
		if (xi == x2i){
			r = (bi - b2i) % m;
			if (r == 0){
			//	a0 ++lower_bound;
				return -1;
			}
			x = algo.rev(r, m)*(ai - a2i) % m;	
			return (x % m + m) % m;
		} else {
			// now we have to go to next i
			i++;
			xi_1 = xi,   ai_1 = ai,   bi_1 = bi  ;
			x2i_2 = x2i, a2i_2 = a2i, b2i_2 = b2i;
			goto step2;
		}						
	}

}
