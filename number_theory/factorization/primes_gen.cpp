#include "sieve.hpp"
#include <iostream>
#include <print>

int main(){
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    unsigned n;
    std::cin >> n;

    for (const unsigned& p : sieve(n)){
        std::print("{}, ", p);
    }

    return 0;
}