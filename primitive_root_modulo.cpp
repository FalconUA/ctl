#include "ctlfac.hpp"
#include "ctlmod.hpp"
#include "typespec.hpp"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <algorithm>

namespace mp {

	template <typename T>
	static void factorize(T const& t, std::vector<T>& res, TypeSpecifications<T> const& algo){
		T r1, r2;
		fac_pollard(t, r1, r2, algo);

		if ((r1 == t)&&(r2==1)) { res.push_back(r1); return ; }
		if ((r2 == t)&&(r1==1)) { res.push_back(r2); return ; }
	
		if (r1 > 1)	factorize(r1, res, algo);
		if (r2 > 1)	factorize(r2, res, algo);	
	}

	template <typename T> 
	static T powmod (T a, T b, T p) {
		T res = 1;
		while (b)
			if (b % 2 == 1)
				res = T (res * 1ll  * a % p),  --b;
		else
			a = T (a * 1ll * a % p),  b = b/2;
		return res;
	}

	template <typename T>
	static T phi (T n, TypeSpecifications<T> const& algo) {
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
	T primitive_root_modulo(T const& n, TypeSpecifications<T> const& algo){

	    if(n == 2) return 1;
	    if(n == 4) return 3;
	    if(n % 4 == 0) return -1;
	    bool ev_numb_n = (n%2 == 0);
	    T t = (ev_numb_n)?(n/2):(n);
	    T p = 0;
		
		{
			std::vector<T> fact;
			factorize(t, fact, algo);
			for (int i=0; i<fact.size(); i++)
				if (fact[0] != fact[i])
					return -1;
			p = fact[0];
		}
	
	    t = phi(n, algo);

	    std::vector<T> rank_G_factors;
	    factorize(t, rank_G_factors, algo);
		
	    int m = (int)rank_G_factors.size();

	    std::vector<T> Group;

	    for(T i = 1, k = ev_numb_n ? 2 : 1 ; i < n; i = i+k)
		    if(i % p != 0)
			    Group.push_back(i);

	    T g,b;
		int group_size = Group.size(), index = 0;

		std::random_shuffle(Group.begin(), Group.end());
	
		
	    for (int index = 0; index < group_size; index++){
	        //1
	        g = Group[index];
	        //2
	        for(int i = 0; i < m; i++){
	            //2.1
	            b = powmod(g, t/rank_G_factors[i], n);
	            //2.2
	            if(b == 1)
	                break;
	        }
	        //3
	        if(b != 1)
	            return g;
	    }
	}

}
