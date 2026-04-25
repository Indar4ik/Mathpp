#pragma once
#include <vector>
#include <cmath>
#include <algorithm>

inline std::vector<unsigned> sieve(unsigned n) noexcept{
    [[assume(n>0)]];
    std::vector<unsigned> lp(n + 1, 0), primes;
    primes.reserve(n / log(n));

    for (unsigned i = 2; i <= n; ++i){
        if (lp[i] == 0){
            lp[i] = i;
            primes.push_back(i);
        }

        unsigned max_p = std::min(n / i, lp[i]);
        for (unsigned p : primes){
            if (p > max_p) break;
            lp[p * i] = p;
        }
    }

    return primes;
}

inline std::vector<unsigned> firstNprimes(unsigned n) noexcept{
    [[assume(n>0)]];
    if (n == 1) return {2};
    if (n == 2) return {2, 3};
    if (n == 3) return {2, 3, 5};

    size_t k = 0;
    unsigned limit = (n < 6) ? 13 : (unsigned)(n * (log(n) + log(log(n))));
    std::vector<unsigned> lp(limit + 1, 0), primes(n);

    primes[0] = 2; primes[1] = 3; primes[2] = 5;
    
    static const unsigned increments[] = {6, 4, 2, 4, 2, 4, 6, 2};
    unsigned i = 7, idx = 0;

    while (i <= limit){
        if (lp[i] == 0){
            lp[i] = i;
            primes[k++] = i;
            if (k == n) break;
        }

        unsigned max_p = std::min(limit / i, lp[i]);
        for (unsigned j = 0; j < k; ++j){
            if (primes[j] > max_p) break;
            lp[primes[j] * j] = primes[j];
        }

        i += increments[idx];
        idx = (idx + 1) & 7;
    }
    
    return primes;
}