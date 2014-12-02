#ifndef CTLLOG_HPP
#define CTLLOG_HPP

#include "typespec.hpp"
namespace mp {

	#define INITLOGALG(x) \
		template <typename T> T x(T const& a, T const& b, T const& m, TypeSpecifications<T> const& algo);

	INITLOGALG(log_primitive)
	INITLOGALG(log_shanks)
	INITLOGALG(log_pollard)
	INITLOGALG(log_index)
	INITLOGALG(log_polhig_hellman)
	INITLOGALG(log_adleman)

	#undef INITLOGALG		
}


#endif
