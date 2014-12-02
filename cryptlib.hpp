#ifndef CRYPTLIB_HPP
#define CRYPTLIB_HPP

#include "ctlfac.hpp"
#include "ctllog.hpp"
#include "ctlmod.hpp"
#include "typespec.hpp"

// костыли костыли
	#include "fac_fermat.cpp"
	#include "fac_pollard.cpp"
	#include "fac_qsieve.cpp"
	#include "fac_lenstra.cpp"

	#include "log_primitive.cpp"
	#include "log_shanks.cpp"
	#include "log_pollard.cpp"
	#include "log_index.cpp"
	#include "log_polhig_hellman.cpp"
	#include "log_adleman.cpp"

	#include "primitive_root_modulo.cpp"

#endif
