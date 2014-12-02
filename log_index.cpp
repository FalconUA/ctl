#include "ctllog.hpp"
#include "typespec.hpp"
#include <algorithm>
#include <iostream>
using std::sqrt;

class MathError{
    virtual void debug_print(){std::cerr<<"Math error happend\n";}
};

class InversionError:public MathError{
    virtual void debug_print(){std::cerr<<"Integer is not invertible\n";}
};

class LogarithmError:public MathError{
    virtual void debug_print(){std::cerr<<"Logarithm error\n";}
};

template <typename Integer> static bool newisPrime(Integer N){
    Integer M = sqrt(N)+1;
    for(Integer i = 2; i < M; i = i+Integer(1)){
        if (N%i == 0) return false;
    }
    return true;
}

template <typename Integer> static Integer newnewPrime(bool reset = false){
    static Integer prime = 1;
    if (reset) prime = 1;
    do{
        prime = prime+Integer(1);
    } while (!newisPrime(prime));
    return prime;
}



namespace mp {

template <typename Integer> Integer log_index(Integer const& A, Integer const& B, Integer const& P, TypeSpecifications<Integer> const& algo){
	Integer a=A, b=B, p=P;
    const size_t factorBaseLength = (std::log2(p)/std::log2(10))+1;
    std::vector<Integer> factorBase;
    std::vector<Integer> factorBaseLog;

    factorBase.push_back(newnewPrime<Integer>(true));
    while (factorBase.size() != factorBaseLength) factorBase.push_back(newnewPrime<Integer>());

    std::vector< std::vector<Integer> > relations_list, rt;
    Integer Q(a), C(1);

    // Vector normalization: (0, 0, ..., 0, xi, ..., xn) -> (0, 0, ..., 0, 1, ..., modInverse(xi, P)*xn mod P)
    auto normalize = [](std::vector<Integer> &V, const Integer& P){
        auto rep = std::find_if(V.begin(), V.end(), std::bind2nd(std::not_equal_to<Integer>(), Integer(0)));
        if (rep == V.end()) return;
        try{
            Integer T = modInverse(*rep, P);
            std::for_each(rep, V.end(), [&T, &P](Integer& a){a = (a*T)%P;});
        } catch (InversionError){
            V.assign(V.size(), Integer(0));
        }
    };

    // WHY DO I THROW?..
    auto factorize = [](Integer Q, std::vector<Integer>& base){
        std::vector<Integer> factorization(base.size());
        for(size_t i = 0; i < base.size() && Q != 1; i++){
            while (Q != 1 && Q%base[i] == Integer(0)){
                factorization[i]++;
                Q /= base[i];
            }
        }
        if (Q == Integer(1)) throw factorization;
    };

    // Push back vector if it's independent
    auto push_back_independent = [p, normalize](std::vector< std::vector<Integer> >& V, std::vector<Integer> E){
        Integer P = p - Integer(1);
        for(auto it = V.begin(); it != V.end(); it++){
            auto rep = std::find_if(it->begin(), it->end(), std::bind2nd(std::not_equal_to<Integer>(), Integer(0)));
            auto pep = *(E.begin() + (rep-it->begin()));
            std::transform(E.begin(), E.end(), it->begin(), E.begin(), [&pep, &P](Integer first, Integer& second){
                first = (first - pep*second)%(P);
                if (first < Integer(0)) first += P;
                return first;
            });
        }
        normalize(E, P);
        if (std::none_of(E.begin(), E.end()-1, std::bind2nd(std::not_equal_to<Integer>(), Integer(0)))) return;
        V.push_back(E);
    };

    // ...BECAUSE I CAN!
    while (relations_list.size() < factorBaseLength){
        try{
            factorize(Q, factorBase);
        }
        catch (std::vector<Integer> factor){
            factor.push_back(C);
            push_back_independent(relations_list, factor);
        }
        Q = (Q*a)%p;
        C = C+Integer(1);
    }
    // IT'S MAGIC

    std::sort(relations_list.begin(), relations_list.end());
    for(auto it = relations_list.begin(); it != relations_list.end(); it++) push_back_independent(rt, *it);
    relations_list.assign(rt.rbegin(), rt.rend());


    // Relations list debug output
//    std::for_each(relations_list.begin(), relations_list.end(), [](std::vector<Integer> t){
//        std::for_each(t.begin(), t.end(), [](Integer t){std::cerr<<t<<" ";});
//        std::cerr<<std::endl;
//    });

    for(size_t i = 0; i < factorBaseLength; i++) factorBaseLog.push_back(relations_list.at(i).back());

    Q = a*b, C = Integer(1);

    for(;;){
        try{
            factorize(Q, factorBase);
        }
        catch (std::vector<Integer> factor){
            return (std::inner_product(factor.begin(), factor.end(), factorBaseLog.begin(), Integer(0)) - C)%(p-Integer(1));
        }
        Q = (Q*a)%p;
        C = C + Integer(1);
    }

    throw LogarithmError();
}



}
