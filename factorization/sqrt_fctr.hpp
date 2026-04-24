#pragma once
#include <vector>
#include <bit>

inline std::vector<unsigned long> sqrt_fctr(unsigned long n) noexcept{
    if (n == 0) return {};
    [[assume(n > 0)]];
    int twos = std::countr_zero(n);
    n >>= twos;
    std::vector<unsigned long> factors(twos, 2UL);
    // factors.reserve(64);

    for (unsigned long i = 3; i * i <= n; i += 2){
        while (n % i == 0){
            factors.push_back(i);
            n /= i;
        }
    }
    if (n > 1) factors.push_back(n);
    
    return factors;
}