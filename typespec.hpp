#ifndef TYPESPEC_HPP
#define TYPESPEC_HPP
#include <cstddef>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdlib.h>

template <typename T> static inline T		mpstdgcd(T, T);
template <typename T> static inline T		mpstdabs(T);
template <typename T> static inline T		mpstdrev(T, T);
template <typename T> static inline size_t	mpstdnumlen(T);
template <typename T> static inline bool    mpstdisprime(T);

namespace mp {	

	template <typename T>
	class TypeSpecifications
	{
	public:
		T		(*sqrt)		(T n);
		T		(*gcd)		(T lhs, T rhs);
		T		(*abs)		(T n);
		T		(*rev)		(T a, T m);
		T		(*pow)		(T n, T a);
		T		(*random)	(T range); 
		bool	(*isprime)	(T a);
		size_t  (*numlen)	(T a);

		TypeSpecifications(): 
			sqrt		(NULL), 
			gcd			(&mpstdgcd<T>), 
			rev			(&mpstdrev<T>),
			abs			(&mpstdabs<T>),
			pow			(NULL),
			random		(NULL),
			isprime		(&mpstdisprime<T>),
			numlen		(&mpstdnumlen<T>){};
	};
};

template <typename T>
static inline T mpstdgcd(T a, T b){
	T c; 
	while (a!=0){c=a; a=b%a; b=c;}
	return b;
}

template <typename T>
static inline T mpstdabs(T x){return (x < 0)?(-x):(x);}

template <typename T> 
static inline void gcdex(T const& a, T const& b, T& d, T& x, T& y){
	if (b==0) {d=a; x=1; y=0; return;}
	T d1, x1, y1;
	gcdex(b, a % b, d1, x1, y1);
	d = d1, x = y1, y = x1 - (a/b)*y1;
}

template <typename T>
static inline T mpstdrev(T a, T m){
	T x, y, d;
	gcdex(a, m, d, x, y);
	if (d != 1) return -1;
	else 	
		return (x % m + m) % m;
}

template <typename T>
static inline size_t mpstdnumlen(T a){
	size_t ans = 0;
	if (a==0) return 1;
	while (a>0){
		ans++;
		a = a/10;
	}
	return ans;	
}
template <typename T>
static inline bool mpstdisprime(T a){
	std::ifstream primesfile;
	primesfile.open("primes.dat");
	if (!primesfile.is_open()){
		primesfile.close();
		return false;
	}
	size_t n;
	primesfile >> n; 
	T p = 0;
	for (size_t i=0; (i<n)&&(p < a); i++){		
		primesfile >> p;
		if (p == a){
			primesfile.close();
			return true;
		}
	}
	
	primesfile.close();
	return false;
}


#endif
