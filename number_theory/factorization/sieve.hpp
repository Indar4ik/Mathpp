#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <ranges>

inline std::vector<uint32_t> sieve(uint32_t n) noexcept {
    [[assume(n > 0)]];
    std::vector<uint32_t> lp((n >> 1) + 1, 0), primes = {2, 3};
    primes.reserve(n / log(n));

    uint32_t i = 5;
    const uint32_t max_idx = ((n - 1) & ~1u) / 3;
    for (uint32_t idx = 1; idx <= max_idx; ++idx){
        if (lp[idx] == 0){
            lp[idx] = i;
            primes.push_back(i);
        }

        uint32_t max_p = std::min(n / i, lp[idx]);
        for (const uint32_t& p : primes | std::views::drop(2)){
            if (p > max_p) break;
            lp[p * i / 3] = p;
        }

        i += 6 - ((i % 3) << 1);
    }

    return primes;
}

inline std::vector<uint32_t> firstNprimes(uint32_t n) noexcept {
    [[assume(n > 0)]];
    if (n == 1) return {2};

    size_t k = 2;
    uint32_t limit = (n < 6) ? 13 : (uint32_t)(n * (log(n) + log(log(n))));
    std::vector<uint32_t> lp((limit >> 1) + 1, 0), primes(n);

    primes[0] = 2; primes[1] = 3;

    uint32_t i = 5;
    const uint32_t max_idx = ((limit - 1) & ~1u) / 3;
    for (uint32_t idx = 1; idx <= max_idx; ++idx){
        if (lp[idx] == 0){
            lp[idx] = i;
            primes[k++] = i;
            if (k == n) break;
        }

        uint32_t max_p = std::min(limit / i, lp[idx]);
        for (const uint32_t& p : primes | std::views::drop(2) | std::views::take(k - 2)){
            if (p > max_p) break;
            lp[p * i / 3] = p;
        }

        i += 6 - ((i % 3) << 1);
    }

    return primes;
}