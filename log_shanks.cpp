#include "ctllog.hpp"
#include "typespec.hpp"
#include <map>

namespace mp {

	template <typename T> 
	T log_shanks(T const& a, T const& b, T const& m, TypeSpecifications<T> const& algo)
	{
		T n = algo.sqrt(m) + 1;
 
		T an = 1;
		for (T i=0; i<n; ++i)
			an = (an * a) % m;
	 
		std::map<T, T> vals;
		for (T i=1, cur=an; i<=n; ++i) {
			if (!vals.count(cur))
				vals[cur] = i;
			cur = (cur * an) % m;
		}
	 
		for (T i=0, cur=b; i<=n; ++i) {
			if (vals.count(cur)) {
				T ans = vals[cur] * n - i;
				if (ans < m)
					return ans;
			}
			cur = (cur * a) % m;
		}
		return -1;
	}

}
