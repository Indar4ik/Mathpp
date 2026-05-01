#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <ranges>

inline std::vector<uint32_t> sieve(uint32_t n) noexcept {
    [[assume(n > 0)]];
    std::vector<uint32_t> lp((n >> 1) + 1, 0), primes = {2};
    primes.reserve(n / log(n));

    for (uint32_t i = 3; i <= n; i += 2){
        uint32_t idx = i >> 1;
        if (lp[idx] == 0){
            lp[idx] = i;
            primes.push_back(i);
        }

        uint32_t max_p = std::min(n / i, lp[idx]);
        for (const uint32_t& p : primes | std::views::drop(1)){
            if (p > max_p) break;
            lp[p * i >> 1] = p;
        }
    }

    return primes;
}

inline std::vector<uint32_t> firstNprimes(uint32_t n) noexcept {
    [[assume(n > 0)]];
    if (n == 1) return {2};

    size_t k = 1;
    uint32_t limit = (n < 6) ? 13 : (uint32_t)(n * (log(n) + log(log(n))));
    std::vector<uint32_t> lp((limit >> 1) + 1, 0), primes(n);

    primes[0] = 2;

    for (uint32_t i = 3; i <= limit; i += 2){
        uint32_t idx = i >> 1;
        if (lp[idx] == 0){
            lp[idx] = i;
            primes[k++] = i;
            if (k == n) break;
        }

        uint32_t max_p = std::min(limit / i, lp[idx]);
        for (uint32_t j = 1; j < k; ++j){
            if (primes[j] > max_p) break;
            lp[primes[j] * i >> 1] = primes[j];
        }
    }

    return primes;
}