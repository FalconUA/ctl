#ifndef CTLMOD_HPP
#define CTLMOD_HPP

#include "typespec.hpp"

namespace mp{
	
	template <typename T> 
	T primitive_root_modulo(T const& n, TypeSpecifications<T> const& algo);
}

#endif
