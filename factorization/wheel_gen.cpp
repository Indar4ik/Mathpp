#include "sieve.hpp"
#include <iostream>
#include <print>

int main(){
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    unsigned int n;
    unsigned long p = 1;
    std::cin >> n; // Количество чисел в колесе

    for (const unsigned& prime : firstNprimes(n)){
        p *= prime;
    }

    std::vector<unsigned> primes_to_p = sieve(p + 1);

    std::print("{{");
    for (unsigned i = n + 1; i < primes_to_p.size(); ++i){
        std::print("{}, ", primes_to_p[i] - primes_to_p[i - 1]);
    }
    std::println("{}}}", p + primes_to_p[n] - primes_to_p.back());

    std::println("i = {}", primes_to_p[n]);
    std::println("idx = (idx + 1) % {}", primes_to_p.size() - n);
    std::println("p={}", p);

    return 0;
}