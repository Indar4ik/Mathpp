#pragma once
#include <vector>
#include <bit>

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

inline void wheel_sqrt_fctr(unsigned long n, std::vector<unsigned long>& factors) noexcept{
    factors.clear();
    if (n <= 1) return;

    int twos = std::countr_zero(n);
    n >>= twos;
    for (int i = 0; i < twos; ++i) factors.push_back(2);
    while (n % 3 == 0) { factors.push_back(3); n /= 3; }
    while (n % 5 == 0) { factors.push_back(5); n /= 5; }
    while (n % 7 == 0) { factors.push_back(7); n /= 7; }
    // factors.reserve(64);

    static const unsigned long incremenst[] = {};
    unsigned long i = 11;
    int idx = 0;

    while (i * i < n){
        while (true){
            unsigned long q = n / i;
            unsigned long r = n % i;
            if (r != 0) break;
            factors.push_back(i);
            n = q;
        }
        i += incremenst[idx];
        idx = (idx + 1) & 1;
    }
    if (n > 1) factors.push_back(n);
}