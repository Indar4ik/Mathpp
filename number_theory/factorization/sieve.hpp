#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <ranges>

inline uint32_t idx_to_x (uint32_t idx) noexcept {
    return (0x62642424U >> ((idx & 7) << 2)) & 15;
}

inline uint32_t x_to_idx (uint32_t x) noexcept {
    uint32_t x_offset = x - 7;

    uint32_t k = x_offset / 30;
    uint32_t rem = x_offset % 30;

    uint32_t r = (0x777655544322100ULL >> ((rem >> 1) << 2)) & 15;
    return (k << 3) | r;
}

inline std::vector<uint32_t> sieve(uint32_t n) noexcept {
    [[assume(n > 0)]];
    if (n < 2) return {};
    if (n < 3) return {2};
    if (n < 5) return {2, 3};
    if (n < 7) return {2, 3, 5};

    const uint32_t max_idx = x_to_idx(n);
    std::vector<uint32_t> lp(max_idx + 1, 0), primes = {2, 3, 5};
    primes.reserve(n / log(n));

    uint32_t i = 7;
    for (uint32_t idx = 0; idx <= max_idx; ++idx){
        if (lp[idx] == 0){
            lp[idx] = i;
            primes.push_back(i);
        }

        uint32_t max_p = std::min(n / i, lp[idx]);
        for (const uint32_t& p : primes | std::views::drop(3)){
            if (p > max_p) break;
            lp[x_to_idx(p * i)] = p;
        }

        i += idx_to_x(idx);
    }

    return primes;
}

inline std::vector<uint32_t> firstNprimes(uint32_t n) noexcept {
    if (n == 0) return {};
    if (n == 1) return {2};
    if (n == 2) return {2, 3};
    if (n == 3) return {2, 3, 5};

    size_t k = 3;
    uint32_t limit = (n < 6) ? 13 : (uint32_t)(n * (log(n) + log(log(n))));
    const uint32_t max_idx = x_to_idx(limit);

    std::vector<uint32_t> lp(max_idx + 1, 0), primes(n);

    primes[0] = 2; primes[1] = 3; primes[2] = 5;

    uint32_t i = 7;
    for (uint32_t idx = 0; idx <= max_idx; ++idx){
        if (lp[idx] == 0){
            lp[idx] = i;
            primes[k++] = i;
            if (k == n) break;
        }

        uint32_t max_p = std::min(limit / i, lp[idx]);
        for (const uint32_t& p : primes | std::views::drop(3) | std::views::take(k - 3)){
            if (p > max_p) break;
            lp[x_to_idx(i * p)] = p;
        }

        i += idx_to_x(idx);
    }

    return primes;
}