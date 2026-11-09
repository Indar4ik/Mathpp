#pragma once
#include <cassert>
#include <span>
#include <utility>

#include "redc.hpp"

namespace hpc {

/**
 * @brief Вычисляет a^(-1) mod m для нечетного m и известного m_inv = -m^(-1) mod 2^64.
 *        Основное вычислительное ядро без повторного вычисления m_inv.
 */
template <std::unsigned_integral W>
[[nodiscard]] constexpr W inv_mod_odd(W a, W m, W m_inv) noexcept {
    [[assume((m & W{1}) == W{1})]];
    [[assume(m > W{1})]];
    using D = detail::wider_uint_t<W>;

    if (a >= m) a %= m;
    if (a == 0) return 0; // Обратного не существует
    if (a == 1) return 1;

    W u = m;
    W v = a;
    W A = 0;
    W C = 1;

    // Убираем младшие нули v
    int tz = std::countr_zero(v);
    if (tz > 0) {
        v >>= tz;
        W q = m_inv & ((W{1} << tz) - W{1});
        C = static_cast<W>((static_cast<D>(q) * m + D{1}) >> tz);
        if (v == 1) return C;
    }

    // Главный цикл: u и v всегда нечётные
    while (u != v) {
        if (u < v) {
            std::swap(u, v);
            std::swap(A, C);
        }

        u -= v;
        A = (A >= C) ? (A - C) : (A + m - C);

        // u стал чётным, сдвигаем все нули
        int t = std::countr_zero(u);
        u >>= t;

        W q = (A * m_inv) & ((W{1} << t) - W{1});
        A = static_cast<W>((static_cast<D>(q) * m + A) >> t);

        if (u == 1) return A;
    }

    // Если сошлись к числу > 1, значит gcd(a, m) > 1, обратного нет
    return (u == 1) ? A : 0;
}

template <std::unsigned_integral W>
[[nodiscard]] constexpr W inv_mod_odd(W a, W m) noexcept {
    [[assume((m & W{1}) == W{1})]];
    if (m == 1) return 0; // по модулю 1 единственный вычет - 0
    
    const W m_inv = -inv_mod2(m);
    return inv_mod_odd(a, m, m_inv);
}

/**
 * @brief Обращает массив вычетов в Montgomery-форме одним обращением (трюк Монтгомери).
 * @param in  Входы в Montgomery-форме: in[i] = a_i * R mod m.
 * @param out Выход в Montgomery-форме: out[i] = a_i^(-1) * R mod m.
 * @note Все a_i обязаны быть обратимы: один необратимый элемент обнуляет общее
 *       произведение, и тогда ВСЕ выходы становятся нулями.
 * @note in и out не должны пересекаться: прямой проход затирает out префиксами,
 *       а обратному нужны исходные in[i].
 */
template <std::unsigned_integral W>
void batch_inverse(std::span<const W> in, std::span<W> out, const montgomery<W> &M) {
    const size_t n = in.size();
    if (n == 0) return;
    assert(out.size() >= n);

    // 1. Прямой проход: префиксные произведения P[i] = in[0] * ... * in[i] mod m
    // Префиксы храним in-place прямо в out
    out[0] = in[0];
    for (size_t i = 1; i < n; ++i) {
        out[i] = M.mul(out[i - 1], in[i]);
    }

    // 2. Единственное обращение общего произведения. M.m_inv хранит +m^(-1), а ядру
    //    нужен -m^(-1) - берём готовый вместо повторного метода Ньютона.
    W inv = M.to(inv_mod_odd(M.from(out[n - 1]), M.m, static_cast<W>(W{0} - M.m_inv)));

    // 3. Обратный проход
    for (size_t i = n - 1; i > 0; --i) {
        W prefix = out[i - 1];            // in[0] * ... * in[i-1]
        out[i] = M.mul(inv, prefix); // in[i]^(-1) = inv * prefix
        inv = M.mul(inv, in[i]);     // отсекаем in[i], оставляем inv для in[0...i-1]
    }
    out[0] = inv;
}

} // namespace hpc
