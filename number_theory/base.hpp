#pragma once
#include <bit>
#include <concepts>
#include <cstdint>

namespace hpc {

// x^{-1} mod 2^W для нечётного x, метод Ньютона.
// Старт y = (3*x)⨁2 точен на 5 бита , каждый шаг
// удваивает: 5 -> 10 -> 20 -> 40 -> 80. Значит 3 шага хватает для 32 бит, 4 для 64.
template <std::unsigned_integral T>
[[nodiscard]] constexpr T inv_mod2(T x) noexcept {
    T y = (T{3} * x) ^ T{2};
    constexpr int itr = std::countr_zero(sizeof(T));
    for (int i = 0; i <= itr; ++i) y *= T{2} - x * y;
    return y;
}

/**
 * @brief Вычисляет a^e методом бинарного возведения в степень.
 * @param a Основание степени (тип T)
 * @param e Показатель степени (любой беззнаковый тип E)
 */
template <std::integral T, std::unsigned_integral E>
[[nodiscard]] constexpr T pow(T a, E e) noexcept {
    T res{1};

    while (e > 0) {
        if (e & 1) res *= a;
        e >>= 1;
        if (e == 0) [[unlikely]] break;
        a *= a;
    }

    return res;
}

namespace detail {

template <std::unsigned_integral T>
consteval auto get_wider_type() {
    if constexpr (sizeof(T) == 1) return std::uint16_t{};
    else if constexpr (sizeof(T) == 2) return std::uint32_t{};
    else if constexpr (sizeof(T) == 4) return std::uint64_t{};
    else if constexpr (sizeof(T) == 8) return static_cast<__uint128_t>(0);
    else static_assert(sizeof(T) <= 8, "Unsupported type size for wider_uint.");
}

template <std::unsigned_integral T>
using wider_uint_t = decltype(get_wider_type<T>());

} // namespace detail

/**
 * @brief Вычисляет (a * b) % m.
 * @note Гарантируется отсутствие переполнения для типов до 64 бит включительно
 */
template <std::unsigned_integral T>
[[nodiscard]] constexpr T mulmod(T a, T b, T m) noexcept {
    [[assume(m != 0)]];

    using W = detail::wider_uint_t<T>;
    return static_cast<T>((static_cast<W>(a) * static_cast<W>(b)) % m);
}

/**
 * @brief Вычисляет (a ^ e) % m методом бинарного возведения в степень.
 * @param a Основание степени (тип T)
 * @param e Показатель степени (любой беззнаковый тип E)
 * @param m Модуль (тип T)
 */
template <std::unsigned_integral T, std::unsigned_integral E>
[[nodiscard]] constexpr T powmod(T a, E e, T m) noexcept {
    [[assume(m != 0)]];
    if (m == 1) [[unlikely]] return 0;

    T res{1};
    a %= m;

    while (e > 0) {
        if (e & 1) res = mulmod(res, a, m);
        e >>= 1;
        if (e == 0) [[unlikely]] break;
        a = mulmod(a, a, m);
    }

    return res;
}

} // namespace hpc