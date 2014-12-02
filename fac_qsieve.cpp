#include <cmath>
#include <cstdlib>
#include <bitset>
#include <set>
#include <vector>
#include <algorithm>
#include <map>
#include <utility>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

#include "ctllog.hpp"
#include "typespec.hpp"

template <typename Integer> static Integer pow(Integer val, Integer deg){
    if (deg == Integer(0)) return Integer(1);
    if (deg%2 == 0) return pow(val, deg/2) * pow(val, deg/2);
    else return pow(val, deg/2)*pow(val, deg/2)*val;
}

template <typename Integer> static Integer modPow(Integer val, Integer deg, Integer mod){
    if (deg == Integer(0)) return Integer(1);
    else if (deg%2 == 0) return (modPow(val, deg/Integer(2), mod) * modPow(val, deg/Integer(2), mod)) % mod;
    else return (((modPow(val, deg/Integer(2), mod)*modPow(val, deg/Integer(2), mod))%mod)*val) % mod;
}

template<typename Integer> Integer GCD (Integer a, Integer b) { return b ? GCD (b, a % b) : a;}
 
template <typename Integer> static bool isPrime(Integer N){
    Integer M = sqrt(N)+1;
    for(Integer i = 2; i < M; i = i+Integer(1)){
        if (N%i == 0) return false;
    }
    return true;
}

template <typename Integer> static Integer newPrime(bool reset = false){
    static Integer prime = 1;
    if (reset) prime = 1;
    do{
        prime = prime+Integer(1);
    } while (!isPrime(prime));
    return prime;
}

namespace mp {

template <typename Integer> Integer factorQuadratic(Integer N){
    const int FCount = 6; // YAY, MAGIC NUMBERS!
    const int TCount = 6; // YAY, AWESOME CONSTANTS
    Integer F[FCount] = {-1};
    Integer prime = newPrime<Integer>(true);
    for(int i = 1; i < FCount;){
        if (modPow(N, (prime-Integer(1))/Integer(2), prime) == Integer(1)) F[i++] = prime;
        prime = newPrime<Integer>();
    }

    Integer M = sqrt(N);
    std::vector< std::pair< std::pair<Integer, Integer>, std::pair< std::bitset<FCount>, std::vector<Integer> > > > pairSet;
    Integer X = 0;
    auto Q = [M, N](Integer X){return (X+M)*(X+M)-N;};
    auto nextX = [](Integer X){if (X < 0) return (-X)+1; else if (X > 0) return -X; else return Integer(1);};
    while (pairSet.size() < TCount+1){
        Integer B = Q(X);
        Integer B_ = B;
        std::bitset<FCount> Bfactormod2;
        std::vector<Integer> Bfactor(FCount);
        if (B < Integer(0)) Bfactormod2.flip(0), B = -B, Bfactor[0] = 1;
        for(int i = 1; i < FCount && B != 1; i++){
            while (B % F[i] == 0) Bfactormod2.flip(i), B/=F[i], Bfactor[i]++;
        }
        if (B == 1){
            pairSet.push_back(std::make_pair(std::make_pair(X+M, B_), std::make_pair(Bfactormod2, Bfactor)));
        }
        X = nextX(X);
    }

    for(int i = 1; i < 1<<(TCount+1); i++){
        std::bitset<FCount> res;
        std::bitset<TCount+1> mask(i);
        for(int j = 0; j < TCount+1; j++){
            if (mask[j]) res ^= pairSet[j].second.first;
        }
        if (res == 0){
            Integer x(1);
            for(int j = 0; j < TCount+1; j++) if (mask[j]) x = (x * pairSet[j].first.first) % N;
            std::vector<Integer> l(FCount);
            for(int i = 0; i < FCount; i++){
                for(int j = 0; j < TCount+1; j++){
                    if (mask[j]) l[i] += pairSet[j].second.second[i];
                }
                l[i] /= 2;
            }

            Integer y = 1;
            for(int i = 0; i < FCount; i++) y = (y * modPow(F[i], l[i], N))%N;

            if (x != y && x != -y){
                Integer T = x-y;
                Integer res = GCD(x-y, N);
                if (res == Integer(1) || res == N || res == Integer(-1)) continue;
                return res;
            }
        }
    }

    return Integer(1);
}

template <typename T>
void fac_qsieve(T const& n, T& r1, T& r2, TypeSpecifications<T> const& algo){
	T N = n;
	r1 = factorQuadratic(N);
	r2 = N / r1;
}

} // end of namespace mp
