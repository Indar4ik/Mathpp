#include "sieve.hpp"
#include <iostream>
#include <print>

int main(){
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    uint32_t n;
    std::cin >> n;

    for (const uint32_t& p : sieve(n)){
        std::print("{}, ", p);
    }

    return 0;
}