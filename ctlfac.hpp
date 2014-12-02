#ifndef CTLLOG_HPP
#define CTLLOG_HPP

#include "typespec.hpp"
namespace mp {

	#define INITFACALG(x) \
		template <typename T> void x(T const& n, T& r1, T& r2, TypeSpecifications<T> const& algo);

	INITFACALG(fac_fermat)
	INITFACALG(fac_pollard)
	INITFACALG(fac_qsieve)
	INITFACALG(fac_rgenerator)
	INITFACALG(fac_lenstra)

	#undef INITFACALG
}

#endif
