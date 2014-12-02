#include "ctllog.hpp"
#include "typespec.hpp"

namespace mp {

	template <typename T> 
	T log_primitive(T const& a, T const& b, T const& m, TypeSpecifications<T> const& algo)
	{
		T x = 0, t = 1, z = b % m;
		if (b==1) return 0;
		if (a==b) return 1;
		for (; ++x < m; ){
			t = (t * a) % m;
			if (t == z) return x; 
		}
		return -1;
	}


}
