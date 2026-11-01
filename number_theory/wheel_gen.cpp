#include "sieve.hpp"
#include <iostream>
#include <numeric>
#include <print>
#include <vector>

int main(){
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    uint32_t n;
    std::cin >> n;

    // Колесо с базисом длины n будет пропускать (1 - phi(p_n#)/p_n#) * 100% чисел, где phi - функция Эйлера, p_n# - примориал n-ого простого числа
    // Выражение стремится к 100% (потому что на бесконечности плотность простых бесконечно мала), но размер массива офсетов растёт грубо как (nlnn)^(n-1)
    uint64_t primorial = 1;
    double checked_fraction = 1.0;
    for (uint32_t prime : firstNprimes(n)) {
        primorial *= prime;
        checked_fraction *= (1.0 - 1.0 / static_cast<double>(prime));
    }

    uint64_t start = n;
    while (std::gcd(start, primorial) != 1) ++start;

    std::vector<uint64_t> coprimes;
    for (uint64_t i = start; i < start + primorial; ++i) if (std::gcd(i, primorial) == 1) coprimes.push_back(i);

    std::print("{{");
    for (size_t i = 1; i < coprimes.size(); ++i) std::print("{}, ", coprimes[i] - coprimes[i - 1]);
    std::println("{}}}", start + primorial - coprimes.back());

    std::println("i = {}", start);
    std::println("idx = (idx + 1) % {}", coprimes.size());
    std::println("primorial = {}", primorial);
    std::println("wheel size = phi(primorial) = {}", coprimes.size());
    std::println("checked fraction = {:.5f}", checked_fraction);

    return 0;
}
