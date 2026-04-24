#pragma once
#include <vector>
#include <cmath>

inline std::vector<unsigned long> sieve(unsigned long n) noexcept{
    std::vector<unsigned long> lp(n + 1, 0), primes;
    primes.reserve(n / log(n));

    for (unsigned long i = 2; i <= n; ++i){
        if (lp[i] == 0){
            lp[i] = i;
            primes.push_back(i);
        }
        for (const unsigned long& p : primes){
            if (p > lp[i] || p > n / i) break;
            lp[p * i] = p;
        }
    }

    return primes;
}