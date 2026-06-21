#pragma once
#include <bit>
#include <cstdint>
#include <vector>

inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) noexcept {
    [[assume(m > 0)]];
    return static_cast<__uint128_t>(a) * b % m;
}

inline uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) noexcept {
    [[assume(m > 0)]];
    uint64_t r = 1ull;
    while (e){
        if (e & 1ull) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

inline bool is_prime(uint64_t n) noexcept {
    if (n < 2) return false;
    if (n == 2) return true;
    if ((n & 1ull) == 0ull) return false;
    for (const uint64_t a : {3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull}){
        if (n == a) return true;
        if (n % a == 0ull) return false;
    }
    uint64_t d = n - 1;
    const uint64_t s = std::countr_zero(d);
    d >>= s;
    for (const uint64_t a : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull}){
        uint64_t l = powmod(a, d, n);
        if (l == 1ull || l == (n - 1)) continue;
        for (uint64_t i = 0; i < s; ++i){
            bool flag = (l == (n - 1));
            l = mulmod(l, l, n);
            if (l == 1ull && flag) break;
            if (l == 1ull && !flag) return false;
        }
        if (l != 1ull) return false;
    }
    return true;
}