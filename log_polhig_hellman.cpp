#include "ctllog.hpp"
#include "ctlfac.hpp"
#include "typespec.hpp"
#include <vector>
#include <iostream>

template <typename T>
static void factorize(T& t, std::vector<T>& res, mp::TypeSpecifications<T> const& algo){
	T r1, r2;
	fac_pollard(t, r1, r2, algo);
	
	if ((r1 == t)&&(r2==1)) { res.push_back(r1); return ; }
	if ((r2 == t)&&(r1==1)) { res.push_back(r2); return ; }

	if (r1 > 1)	factorize(r1, res, algo);
	if (r2 > 1)	factorize(r2, res, algo);	
}

template <typename T>
static T phi (T n) {
	T result = n;
	for (T i=2; i*i<=n; ++i)
	if (n % i == 0) {
			while (n % i == 0)
				n /= i;
			result -= result / i;
		}
	if (n > 1)
		result -= result / n;
	return result;
}

template <typename T> 
static T exp_mod (T a, T b, T p) {
	T res = 1;
	while (b)
		if (b % 2 == 1)
			res = T (res * 1ll  * a % p),  --b;
	else
		a = T (a * 1ll * a % p),  b = b/2;
	return res;
}

namespace mp {

	template <typename T>	
	T log_polhig_hellman(T const& G, T const& B, T const& N, TypeSpecifications<T> const& algo){
		//1.
		T g=G, b=B, n=N;

		T m = phi(n);
	    std::vector<T> base;
		std::vector<unsigned> exp_base;
	    std::cout << "m=" << m << std::endl;
		factorize(m, base, algo);
		for(unsigned i = 0; i < base.size(); i++)
			std::cout << base[i] << ", ";

	    std::vector<T> base_tmp = base;
		base.clear();
	    while(!base_tmp.empty()) {
	        unsigned new_size = 0, old_size = base_tmp.size(), deg = 0;
	        T p = base_tmp[0];
	        for(unsigned i = 0; i < old_size; i++) {
	            if(p != base_tmp[i])
	                base_tmp[new_size++] = p;
	            else
	                deg++;
	        }
	        base_tmp.resize(new_size);
	        base.push_back(p);
	        exp_base.push_back(deg);
	    }

	    //debug start
	    std::cout << "prime: ";
	    for(unsigned i = 0; i < base.size(); i++)
	        std::cout << base[i] << ", ";
	    std::cout << "\nexp: ";
	    for(unsigned i = 0; i < base.size(); i++)
	        std::cout << exp_base[i] << ", ";
	    system("pause");
	    //debug end
	
	    //2.
	    unsigned size = base.size();
	    //Table start
	    std::vector< std::vector<T> > table;
	    table.resize(size);
	    T t;
	    for(unsigned i = 0; i < size; i++) {
	        t = m / base[i];
	        for(unsigned j = 0; j < base[i]; j++)
	            table[i].push_back(exp_mod(g,T(t*j),n));
	    }
	    //table end
	    T *x = new T[size];
	    T *v1 = new T[size], *v2 = new T[size];
	    for(unsigned i = 0; i < size; i++) {
	        static T w,q;
	        q = 1;
	        for(unsigned j = 0; j < exp_base[i]; j++) {
	            t = 0;
	            w = 1;
	            for(unsigned k = 0; k < j; k++) {
	                t = t + x[k] * w;
	                w = w * base[i];
	            }
	            w = 1;
	            for(unsigned k = 0; k <= j; k++)
	                w = w * base[i];
	            T gt = 1;
	            for(unsigned k = 0; k < t; k++)
	                gt = gt * g;
	            t = (exp_mod(b,m/w,n) * algo.rev(exp_mod(gt,m/w,n),n)) % n;
	            for(unsigned k = 0, s = table[i].size(); k < s; k++)
	                if(t == table[i][k])
	                    x[j] = (int)k;
	            v1[i] = v1[i] + x[j] * q;
	            T p_exp = 1;
	            for(unsigned k = 0; k < exp_base[i]; k++)
	                p_exp = p_exp * base[i];
	            v1[i] = v1[i] % p_exp;
	            q = q * base[i];
	        }
	    }
	    for(unsigned i = 0; i < size; i++) {
	        T p_exp = 1;
	        for(unsigned j = 0; j < exp_base[i]; j++)
	            p_exp = p_exp * base[i];
	        base[i] = p_exp;
	    }
	
	    T s1,s2;
	    for(unsigned i = 0; i < size; i++) {
	        v2[i] = v1[i];
	        for(unsigned j = 0; j < i; j++) {
	            v2[i] = (v2[i] - v2[j]) * algo.rev(base[j],base[i]);
	            s2 = v2[i];
	            v2[i] = v2[i] % base[i];
				if(v2[i] < 0)
                v2[i] = v2[i] + base[i];
			}
		}
		s1 = v2[0];
		t = 1;
		for(unsigned i = 1; i < size; i++) {
			t = t * base[i-1];
			s1 = s1 + t * v2[i];
		}

		return s1;

	}



}
