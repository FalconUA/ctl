#include "cryptlib.hpp"
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

inline int isqrt(int x){return (int)std::sqrt(x);} 
inline int ipow(int n, int a){return std::pow(n, a);}
inline int irand(int range){return rand() % range;}

int main(int argc, char* argv[]){

	srand(0);

	mp::TypeSpecifications<int> typespec;
	typespec.sqrt	 = &isqrt;
	typespec.pow	 = &ipow;
	typespec.random  = &irand;

	int i = 0;
	if (argc <= 1) return 0;
	while (i < argc){
		if (i==1) std::cout << "\n";	
		if (std::string(argv[i]) == "-r"){
			int n = atoi(argv[++i]);
			std::cout << "Finding primitive root modulo "<<n<<"\n";
			int g = mp::primitive_root_modulo(n, typespec);
			std::cout << "    primitive root modulo: "<< g<< "\n\n";
			int r=1;
			for (int j=0; j<50; j++){
				r = r*g % n;
				std::cout << r << " ";
			}
			std::cout << "\n\n";
		} else
		if (std::string(argv[i]) == "-f"){
			int n = atoi(argv[++i]);
			std::cout << "Factorizing " << n << ":\n";
			int r1, r2;
			mp::fac_fermat(n, r1, r2, typespec); 
			std::cout << "    fermat:    " << r1 << " " << r2 << "\n";
			mp::fac_pollard(n, r1, r2, typespec);
			std::cout << "    pollard:   " << r1 << " " << r2 << "\n";
			mp::fac_qsieve(n, r1, r2, typespec);
			std::cout << "    qsieve:    " << r1 << " " << r2 << "\n";
			mp::fac_lenstra(n, r1, r2, typespec);
			std::cout << "    lenstra:   " << r1 << " " << r2 << "\n\n";			
		} else
		if (std::string(argv[i]) == "-l"){
			int a = atoi(argv[++i]);
			int b = atoi(argv[++i]);
			int m = atoi(argv[++i]);
			std::cout << "Finding solution of equation "<<a<<"^x == "<<b
													<<" (mod "<<m<<")\n";

			std::cout << "    primitive:     " << mp::log_primitive(a, b, m, typespec) << "\n";
			std::cout << "    shanks:        " << mp::log_shanks(a, b, m, typespec) << "\n";
			std::cout << "    pollard:       " << mp::log_pollard(a, b, m, typespec) << "\n";
			std::cout << "    index:         " << mp::log_index(a, b, m, typespec) << "\n";
			std::cout << "    polhg hellman: " << mp::log_polhig_hellman(a, b, m, typespec) << "\n";
			std::cout << "    adleman:       " << mp::log_adleman(a, b, m, typespec) << "\n\n";

		}
		i++;		
	}
	return 0;
}
