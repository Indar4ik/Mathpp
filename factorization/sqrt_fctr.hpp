#pragma once
#include <vector>
#include <bit>
#include <cstdint>

inline void sqrt_fctr(unsigned long n, std::vector<unsigned long>& factors) noexcept{
    factors.clear();
    if (n <= 1) return;

    int twos = std::countr_zero(n);
    n >>= twos;
    for (int i = 0; i < twos; ++i) factors.push_back(2UL);
    // factors.reserve(64);

    for (unsigned long i = 3; i * i <= n; i += 2){
        // while (n % i == 0){
        //     factors.push_back(i);
        //     n /= i;
        // }
        while (true){
            unsigned long q = n / i;
            unsigned long r = n % i;
            if (r != 0) break;
            factors.push_back(i);
            n = q;
        }
    }
    if (n > 1) factors.push_back(n);
}