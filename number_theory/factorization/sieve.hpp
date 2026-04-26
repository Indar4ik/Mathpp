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
    [[assume(n > 0)]];
    if (n == 1) return {2};

    size_t k = 1;
    unsigned limit = (n < 6) ? 13 : (unsigned)(n * (log(n) + log(log(n))));
    std::vector<unsigned> lp(limit / 2+ 1, 0), primes(n);

    primes[0] = 2;

    for (unsigned i = 3; i <= limit; i += 2){
        unsigned idx = i / 2;
        if (lp[idx] == 0){
            lp[idx] = i;
            primes[k++] = i;
            if (k == n) break;
        }

        unsigned max_p = std::min(limit / i, lp[idx]);
        for (unsigned j = 1; j < k; ++j){
            if (primes[j] > max_p) break;
            lp[primes[j] * i / 2] = primes[j];
        }
    }

    return primes;
}