#pragma once
#include <initializer_list>

#include "redc.hpp"

namespace hpc {

// Один раунд Миллера-Рабина в Montgomery-форме: false -> a свидетель составности.
// M построен на нечётном n; mone == (n-1) в форме; вход a в обычном виде.
inline bool mr_witness(uint64_t a, const mont& M, uint64_t minus_one, uint64_t d, int s) noexcept {
    uint64_t x = M.pow(M.to(a), d);                // a^d в форме
    if (x == M.one || x == minus_one) return true;      // a^d == 1 или -1
    for (int r = 1; r < s; ++r) {
        x = M.mul(x, x);
        if (x == minus_one) return true;                // встретили -1 -> раунд пройден
    }
    return false;
}


inline bool mr_witness_lazy(uint64_t a, const mont& M, uint64_t minus_one, uint64_t d, int s) noexcept {
    const uint64_t minus_one_lazy = minus_one + M.m;
    uint64_t x = M.pow_lazy(M.to(a), d);                    // a^d в форме
    if (x == M.one || x == M.one + M.m || x == minus_one || x == minus_one_lazy) return true;      // a^d == 1 или -1
    for (int r = 1; r < s; ++r) {
        x = M.mul_lazy(x, x);
        if (x == minus_one || x == minus_one_lazy) return true;  // встретили -1 -> раунд пройден
    }
    return false;
}

// Детерминированный тест простоты для всего диапазона uint64_t.
inline bool is_prime(uint64_t n) noexcept {
    if (n < 2) return false;
    // Снимаем все делители свидетелей для корректности (иначе n | a -> a^d == 0 mod n)
    // Плюс засекает ~77% составных без умножений
    #pragma clang loop unroll(full)
    for (uint64_t p : {2ull, 3ull, 5ull, 7ull, 13ull, 19ull, 73ull, 193ull, 407521ull, 299210837ull}) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    // n нечётное и не делится на малые простые: раскладываем n - 1 = d * 2^s один раз
    uint64_t d = n - 1ull;
    int s = std::countr_zero(d);
    d >>= s;
    // Montgomery-контекст и -1 в форме считаем один раз на n.
    mont M(n);
    uint64_t minus_one = n - M.one;                     // (n-1) в Montgomery-форме = -1

    const auto all_pass = [&](std::initializer_list<uint64_t> bases) noexcept {
        for (uint64_t a : bases) if (!mr_witness(a, M, minus_one, d, s)) return false;
        return true;
    };
    const auto all_pass_lazy = [&](std::initializer_list<uint64_t> bases) noexcept {
        for (uint64_t a : bases) if (!mr_witness_lazy(a, M, minus_one, d, s)) return false;
        return true;
    };
    // Набор свидетелей по величине n: полный семёрочный нужен только выше 3.2e9,
    // ниже хватает коротких наборов Яшке. Это x2.3 на простых до 2e7 и x1.4 на
    // диапазоне 1e9..3e9; на 64-битных входах ярусы не срабатывают и не мешают.
    // Два младших яруса остаются СТРОГИМИ намеренно: там показатель d короткий, раунд
    // и так заканчивается за десяток умножений, а ленивый вариант на таких n стабильно
    // теряет ~3.5% (замерено раздельными бинарями, не раскладка кода). С третьего яруса
    // ленивый выигрывает: x1.03 на 2.5e7..3.2e9 и x1.06-1.07 на простых ~2^61.
    if (n < 1373653ull)    return all_pass({2, 3});
    if (n < 25326001ull)   return all_pass({2, 3, 5});
    if (n <= LAZY_MAX) {
        if (n < 3215031751ull) return all_pass_lazy({2, 3, 5, 7});
        return all_pass_lazy({2, 325, 9375, 28178, 450775, 9780504, 1795265022});
    }
    if (n < 3215031751ull) return all_pass({2, 3, 5, 7});
    return all_pass({2, 325, 9375, 28178, 450775, 9780504, 1795265022});
}

} // namespace hpc