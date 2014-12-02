#include "ctllog.hpp"
#include "typespec.hpp"

namespace mp {

	template <typename T> 
	T __log_pollard(T const& a, T const& b, T const& m, TypeSpecifications<T> const& algo)
	{
	init: 
		T a0 = 0, b0 = 0, x0 = 1, i = 1, lower_bound = 0, r, x;
		auto f = [&a, &b, &m](T x) -> T {
				return (x % 3 == 0)? (x*x):(
					   (x % 3 == 1)? (x*a):(
					   (x % 3 == 2)? (x*b):(0))) % (m);
			};
		auto g = [&m](T x, T n) -> T {
				return (x % 3 == 0)? (2*n):(
					   (x % 3 == 1)? (n+1):(
					   (x % 3 == 2)? (n  ):(0))) % (m-1);
			};
		auto h = [&m](T x, T n) -> T {
				return (x % 3 == 0)? (2*n):(
					   (x % 3 == 1)? (n  ):(
					   (x % 3 == 2)? (n+1):(0))) % (m-1);			
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
		return -1;
	}

template <typename Integer> static Integer modInverse(Integer n, Integer mod){
    Integer t(0), newt(1);
    Integer r(mod), newr(n);
    Integer quot, temp;
    while (newr != Integer(0)){
        quot = r / newr;

        //(t, newt) = (newt, t - quotient * newt)
        temp = newt;
        newt = t - quot*newt;
        t = temp;

        //(r, newr) := (newr, r - quotient * newr)
        temp = newr;
        newr = r - quot*newr;
        r = temp;
    }
    if (r > Integer(1)) return -1;
    if (t < Integer(0)) t = t+mod;
    return t;
}



template<typename T> 
T log_pollard(T const& A, T const& B, T const& M, TypeSpecifications<T> const& algo){
	T a=A, b=B, p=M;
    T x(1), apow(0), bpow(0), X(1), APOW(0), BPOW(0);
    T S1 = p/T(3); T S2 = T(2)*S1;

    auto next_triplet = [S1, S2, a, b, p](T& x, T& apow, T& bpow){
        if (x <= S1){
            apow = (apow+T(1))%(p-T(1));                                           x = (b*x)%p;
        } else if (x <= S2){
            apow = (T(2)*apow)%(p-T(1)); bpow = (T(2)*bpow)%(p-T(1));  x = (x*x)%p;
        } else {
                                                     bpow = (bpow+T(1))%(p-T(1));  x = (a*x)%p;
        }
    };

    while(true){
        next_triplet(x, apow, bpow);
        next_triplet(X, APOW, BPOW);
        next_triplet(X, APOW, BPOW);

        if (x == X){
            T P = p-T(1);
            T L = (APOW-apow) % P; while (L < T(0)) L += P;
            T R = (bpow-BPOW) % P; while (R < T(0)) R += P;
            T D = GCD(L, P);

            if ((R/D)*D != R) continue;
            T K = (R/D * modInverse(L/D, P/D))%P; while (K < T(0)) K += P;
            T K_ = K;
            K = ((K+P/D)%P);
            while (modPow(a, K, p) != b && K_ != K ) K = ((K+P/D)%P);
            if (modPow(a, K, p) == b) return K;
        }
    }
    return -1;
}

} // end of namespace mp
