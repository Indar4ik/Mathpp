#pragma once
#include <algorithm>
#include <vector>
#include <cmath>
#include <bit>

inline std::vector<unsigned long> sqrt_fctr(unsigned long n) noexcept{
    std::vector<unsigned long> factors;
    factors.reserve(n / log(n));

    std::fill_n(factors.data(), std::countr_zero(n), 2L);
    for (int i = 3; i * i <= n; i += 2){
        while (n % i == 0){
            factors.push_back(i);
            n /= i;
        }
    }
    if (n > 1) factors.push_back(n);
    
    return factors;
}