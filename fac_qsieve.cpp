#include "ctlfac.hpp"
#include "typespec.hpp"
#include <vector>
#include <iostream>


template <typename T> 
static T powmod (T a, T b, T p) {
	T res = 1;
	while (b)
		if (b % 2 == 1)
			res = T (res * 1ll  * a % p),  --b;
	else
		a = T (a * 1ll * a % p),  b = b/2;
	return res;
}

template <typename T>
static inline T LegendreSymbol(T a, T p){ return powmod<T>(a, (p-1)/2, p); }

template <typename T>
static inline T JakobiSymbol(T a, T b, mp::TypeSpecifications<T> const& algo){

	T r, t, c;
step1:
	if (algo.gcd(a, b) != 1) return 0;

step2:
	r = 1;

step3:
	if (a<0){
		a = -a;
		if (b % 4 == 3) r = -r;
	}

step4:
	t = 0;
	while (a % 2 == 0){
		t = t + 1;
		a = a / 2;
	}
	if (t % 2 != 0)
		if ((b % 8 == 3)||(b % 8 == 5))
			r = -r;
	if ((a % 4 == b % 4) && (b % 4 == 3)) 
		r = -r;

step5:
    c = a; a = b % c; b = c;

step6:
	if (a != 0) goto step4;
	return r;	
}

static const int PRIME_N = 60;
static const int PRIME[PRIME_N] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                            31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                            73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
                            127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
                            179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
                            233, 239, 241, 251, 257, 263, 269, 271, 277, 281};


namespace mp{


static bool check_subset(bool *ch_func, int t, int **v){
    int *v_s = new int[t];
    for(int i = 0; i < t; i++)
        v_s[i] = 0;
    for(int i = 0; i <= t; i++)
        if(ch_func[i])
            for(int j = 0; j < t; j++)
                v_s[j] ^= v[i][j];
    for(int j = 0; j < t; j++)
        if(v_s[j] != 0){
            v_s[0] = 1;
            break;
        }
    if(v_s[0] == 0){
        delete[] v_s;
        return true;
    }
    delete[] v_s;
    return false;
}

static bool find_next_subset(bool* last_ch_func, std::vector<int> &subset, unsigned int t, int **v){
    //brute force all subsets
    int ch_j = -1;
    bool skip_first;
    for(int i = t; i >= 0; i--)
        if(last_ch_func[i] == 1){
            ch_j = i;
            skip_first = true;
            break;
        }
    if(ch_j == -1)
        last_ch_func[0] = true, ch_j = 0, skip_first = false;
    int j0, N;
    while(true){
        j0 = ch_j;

        if(!skip_first && check_subset(last_ch_func, t, v)){
            ch_j = -1;//flag, that we found it;
            break;
        }

   //     std::cout << "ch_func: ";
   //     for(int i = 0; i <= t; i++)
   //         std::cout << this->last_ch_func[i];
   //     std::cout << std::endl;

        while(ch_j < (int)t){
            last_ch_func[ch_j++] = false;
            last_ch_func[ch_j] = true;
       //     std::cout << "ch_func: ";
       //     for(int i = 0; i <= t; i++)
       //         std::cout << this->last_ch_func[i];
       //     std::cout << std::endl;
            if(check_subset(last_ch_func,t,v)){
                ch_j = -1;//flag, that we found it;
                break;
            }
        }
        if(ch_j == -1)
            break;

        //int N;
        for(N = 0; (last_ch_func[ch_j] == 1)&&(ch_j > 0); ch_j--)
            last_ch_func[ch_j] = false, N++;

       // if(N > 10)
       //     std::cout << "N = " << N << "\n";
        if((N == (int)t) && last_ch_func[ch_j])
            break;
        else{
            for(ch_j = j0; (last_ch_func[ch_j] != true)&&(ch_j > 0); ch_j--);
            if(last_ch_func[ch_j] != true){
                last_ch_func[ch_j] = true;
                while(ch_j < N)
                    last_ch_func[++ch_j] = true;
            }
            else{
                last_ch_func[ch_j] = false;
                for(int i = 0; i <= N; i++)
                    last_ch_func[++ch_j] = true;
            }
        }
    }
    if(ch_j != -1){//if we didn't found the subset T: Sum[v[i],i<-T] = 0 (in Z_2);
        //res.push_back(*this);
        return false;
    }
    subset.clear();
    //std::vector<int> T;
    //unsigned int k_T = 0;
    for(unsigned int i = 0; i <= t; i++)
        if(last_ch_func[i])
            /*Tau*/subset.push_back(i);//, k_T++;
    return true;//T;
    //delete[] ch_func;
}


static bool find_first_subset(bool* last_ch_func, std::vector<int> &subset, unsigned int t, int **v){
    if (last_ch_func != NULL) delete[] last_ch_func;
    last_ch_func = new bool [t+1];
    for(unsigned int i = 0; i <= t; i++)
        last_ch_func[i] = false;
    return find_next_subset(last_ch_func, subset, t, v);
}

	#define __TRIVIAL_FACTORIZATION__ \
		T N = n;                      \
		if (N == 1){                  \
			r1 = 1, r2 = 1;           \
			return;                   \
		}                             \
		if (N % 2 == 0){              \
			N = N/2; r1 = 2; r2 = N;  \
			return ;                  \
		}                             \
		if (!(algo.isprime == NULL))  \
			if (algo.isprime(N)){     \
				r1 = N, r2 = 1;       \
				return;               \
			}


	template <typename T>
	void fac_qsieve(T const& n, T& r1, T& r2, TypeSpecifications<T> const& algo){

		__TRIVIAL_FACTORIZATION__

		{
			mp::fac_pollard(n, r1, r2, algo);
			return ;
		}

	    //1
	    std::vector<int> prime_base;
	    size_t t = 1, up_lim = algo.numlen(n);
	    prime_base.push_back(-1);
	    for(unsigned int i = 0; i < PRIME_N && t < up_lim; i++)
	        if( JakobiSymbol(n, PRIME[i], algo) == 1)
	            prime_base.push_back(PRIME[i]), t++;
    //2
    T m = algo.sqrt(n);
    //3
    T x = 0, b, tmp;

    T *A = new T[t+1], *B = new T[t+1];
    int** e = new int*[t+1];
    for(unsigned int i = 0; i <= t; i++)
        e[i] = new int[t];
    int** v = new int*[t+1];
    for(unsigned int i = 0; i <= t; i++)
        v[i] = new int[t];

    //int prime_degree[PRIME_N];
    for(unsigned int i = 0, j; i <= t; i++){
        //3.1
        while(true){
            tmp = x + m;
            b = tmp*tmp - (n);
            for(j = 0; j < t; j++)
                e[i][j] = 0;
                //prime_degree[i] = 0;
            if(b < 0)
                e[i][0] = 1;
            for(j = 1, tmp = b; j < t; j++){
                while(tmp % prime_base[j] == 0){
                    tmp = tmp / prime_base[j];
                    e[i][j]++;
                    //prime_degree[j]++;
                }
            }
            if((tmp != 1)&&(tmp != -1)){
                if(x > 0)
                    x = -x;
                else{
                    x = -x;
                    x = x + 1;
                }
                continue;
            }
            break;
        }
        //3.2
        A[i] = x + m;
        B[i] = b;
        for(j = 0; j < t; j++)
            v[i][j] = e[i][j] % 2;
        //3.3
        if(x > 0)
            x = -x;
        else{
            x = -x;
            x = x + 1;
        }
    }

    //4
    //brute force all subsets
    bool* ch_func = NULL;

    std::vector<int> Tau;
    if(find_first_subset(ch_func, Tau,t,v) == false){//if we didn't found the subset T: Sum[v[i],i<-T] = 0 (in Z_2);
        r1 = n;
		r2 = n/r1;
        return ;
    }
    unsigned int k_T = Tau.size();

    T X, Y;
    do{
        //5
        X = 1;
        for(unsigned int i = 0; i < k_T; i++)
            X = (X * A[Tau[i]]) % n;
        //6
        T *l = new T[t];
        for(unsigned int j = 0; j < t; j++){
            T sum;
            for(unsigned int i = 0; i < k_T; i++)
                sum = sum + e[Tau[i]][j];
            l[j] = sum/2;
        }
        //7
        Y = 1;
        for(unsigned int j = 0; j < t; j++){
            T pow_p = 1;
            for(T i = 0; i < l[j]; i = i + 1){
                pow_p = pow_p * prime_base[j];
                while(pow_p < 0)
                    pow_p = pow_p + n;
                pow_p = pow_p % n;
            }
            Y = (Y*pow_p) % n;
        }
        //8
        if((X != Y)&&(X != n - Y))
            break;
        if(find_next_subset(ch_func, Tau ,t ,v) == false){//if we didn't found the subset T: Sum[v[i],i<-T] = 0 (in Z_2);
            //res.push_back(*this);
            r1 = n;
			r2 = n/r1;
			std::cout << "find next subset failed\n";
            return ;
        }
        k_T = Tau.size();
    } while(true);
    if(X < Y){
        T W;
        W = X;
        X = Y;
        Y = W;
    }
    T D = algo.gcd(X-Y, n);
    if((D > 1)&&(D < n)){
        r1 = n / D;
        r2 = D;
		std::cout << "yaaazzz\n";
        return ;
    }
	
    r1 = n;
	r2 = n/r1;
    return ;
}

#undef __TRIVIAL_FACTORIZATION__

}
