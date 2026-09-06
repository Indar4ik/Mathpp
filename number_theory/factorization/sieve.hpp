#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <ranges>

// idx: value
// 0: 4
// 1: 2
// 2: 4
// 3: 2
// 4: 4
// 5: 6
// 6: 2
// 7: 6

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

inline std::vector<uint32_t> sieve_lp(uint32_t n) noexcept {
    if (n < 2) return {};
    
    std::vector<uint32_t> lp(n - 1, 0u), primes;
    primes.reserve(n / log(n));

    for (uint32_t i = 2; i <= n; ++i){
        if (lp[i - 2] == 0){
            lp[i - 2] = i;
            primes.push_back(i);
        }

        // uint32_t max_p = std::min(n / i, lp[i - 2]);
        for (const uint32_t& p : primes){
            if (p > lp[i - 2] || static_cast<uint64_t>(p) * i > static_cast<uint64_t>(n)) break;
            lp[i * p - 2] = p;
        }
    }

    return lp;
}

inline std::vector<uint32_t> EulerFunc(uint32_t n, const std::vector<uint32_t>& lp) noexcept {
    if (n == 0u) return {};
    if (n == 1u) return {0u};

    std::vector<uint32_t> EulerFunc = std::views::iota(1u, n + 1u) | std::ranges::to<std::vector>();
    EulerFunc[0] = 0u;

    for (auto& phi : EulerFunc){
        uint32_t x = phi;
        while (x > 1u){
            uint32_t p = lp[x - 2];
            while (x % p == 0) x /= p;
            phi = phi / p * (p - 1);
        }
    }

    return EulerFunc;
}

inline std::vector<int32_t> MobiusFunc(uint32_t n, const std::vector<uint32_t>& lp) noexcept {
    if (n == 0u) return {};
    if (n == 1u) return {1u};

    std::vector<int32_t> MobiusFunc(n, 1);

    for (uint32_t x = 2u; x <= n; ++x){
        uint32_t mu = x;
        while (mu > 1u){
            uint32_t p = lp[mu - 2];
            mu /= p;
            if (mu % p == 0) {
                MobiusFunc[x - 1] = 0;
                break;
            }
            MobiusFunc[x - 1] ^= 0xfffffffe;
        }
    }

    return MobiusFunc;
}