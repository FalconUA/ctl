#include "ctllog.hpp"
#include "ctlfac.hpp"
#include "typespec.hpp"
#include <vector>
#include <iostream>

template <typename Integer> static Integer modInverse(Integer n, Integer mod){
    Integer t(0), newt(1);
    Integer r(mod), newr(n);
    Integer quot, temp;
    while (newr != Integer(0)){
        quot = r / newr;

        //(t, newt) = (newt, t - quotient * newt)
        temp = newt;
        newt = t - quot*newt;
        t = temp;

        //(r, newr) := (newr, r - quotient * newr)
        temp = newr;
        newr = r - quot*newr;
        r = temp;
    }
    if (r > Integer(1)) return -1;
    if (t < Integer(0)) t = t+mod;
    return t;
}


template<typename Integer> static Integer solveChineseRemainders(std::vector< std::pair<Integer, Integer> > values){
    std::vector<Integer> X(values.size());
    for(size_t i = 0; i < values.size(); i++){
        X[i] = values[i].first;
        for(size_t j = 0; j < i; j++){
            X[i] = modInverse(values[j].second, values[i].second) * (X[i] - X[j]);
            X[i] = X[i] % values[i].second;
            if (X[i] < Integer(0)) X[i] += values[i].second;
        }
    }

    Integer res = 0;
    for(int i = int(values.size())-1; i >= 0; i--){
        res *= values[i].second;
        res += X[i];
    }

    return res;
}

template<typename Integer> static std::vector < std::pair<Integer, size_t> > factorize(Integer A){
    std::vector< std::pair<Integer, size_t> > res;
    for(Integer t(2); t <= A; t = t + Integer(1)){
        auto temp = std::make_pair(t, 0);
        for(;(A%t == 0); A/=t, temp.second++);
        if (temp.second) res.push_back(temp);
    }
    return res;
}




namespace mp {

template <typename T> T log_polhig_hellman(T const& a, T const& b, T const& p, TypeSpecifications<T> const& algo){
    T P = p-T(1);
    std::vector< std::pair<T, size_t> > factor = factorize(P);
    std::vector< std::vector<T> > R(factor.size()+1);
    std::vector< std::pair<T, T> > mods;

    for(int i = 1; i < R.size(); i++)
        for(T j(0); j < factor[i-1].first; j+=T(1))
            R[i].push_back(modPow(a, j*(p-1)/factor[i-1].first, p));

    for(int i = 1; i < R.size(); i++){
        T X(0);
        T Q = factor[i-1].first;
        T Q_(1);
        T P = (p-1);
        T A;

        for(size_t j(0); j < factor[i-1].second; j++){
            P /= Q;
            A = modPow(b*modInverse(modPow(a, X, p) ,p), P, p);
            for(int z = 0; z < R[i].size(); z++){
                if (R[i][z] == A){
                    X += T(z)*Q_;
                    Q_ *= Q;
                    break;
                }
                if (z+1 == R[i].size()) return 0;
            }
        }

        mods.push_back(std::make_pair(X%Q_, Q_));
    }

    return solveChineseRemainders(mods);
}

} // end of namespace mp

