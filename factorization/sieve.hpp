#pragma once
#include <vector>
#include <cmath>

inline std::vector<int> sieve(int n) noexcept{
    std::vector<int> lp(n + 1, 0), primes;
    primes.reserve(n / log(n));

    for (int i = 2; i <= n; ++i){
        if (lp[i] == 0){
            lp[i] = i;
            primes.push_back(i);
        }
        for (const int& p : primes){
            if (p > lp[i] || p > n / i) break;
            lp[p * i] = p;
        }
    }

    return lp;
}