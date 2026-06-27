#pragma once
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
// #include <mkl.h>
#include <numbers>
#include <utility>
#include <vector>

namespace fft {

using cd = std::complex<double>;

inline void zero_pad(std::vector<cd>& arr) noexcept {
    if (arr.empty()) {
        arr.resize(1, cd(0.0, 0.0));
        return;
    }
    size_t n = arr.size();
    [[assume(n > 0)]];

    if ((n & (n - 1)) == 0) return;

    size_t new_size = 1ULL << (64 - std::countl_zero(n - 1));
    arr.resize(new_size, cd(0.0, 0.0));
}

// Передискретизация замкнутого контура до степени двойки БЕЗ искажения формы.
// В отличие от zero_pad (который дописывает нули и заставляет эпициклы "танцевать"
// вокруг нуля, а заодно расстраивает соответствие бин<->гармоника), здесь мы линейно
// интерполируем ту же кривую в N = bit_ceil точек, равномерно по параметру. Кривая та же,
// меняется только число отсчётов, поэтому бин i по-прежнему отвечает гармонике i.
inline void resample_pow2(std::vector<cd>& c, size_t &n) noexcept {
    const size_t P = c.size();
    if (P < 2) return;

    // Если контур замкнут (последняя точка совпадает с первой), отбрасываем дубликат:
    // период образуют Q различных вершин c[0..Q-1] с переходом c[Q-1] -> c[0].
    const size_t Q = (c.front() == c.back()) ? P - 1 : P;
    const size_t N = std::bit_ceil(Q);
    n = N; // выставляем ВСЕГДА — иначе при раннем выходе caller получит старый размер
    if (N == Q) { c.resize(Q); return; }

    std::vector<cd> out(N);
    const double scale = static_cast<double>(Q) / static_cast<double>(N);
    for (size_t j = 0; j < N; ++j) {
        const double u = static_cast<double>(j) * scale; // позиция в [0, Q)
        const size_t i0 = static_cast<size_t>(u);
        size_t i1 = i0 + 1;
        if (i1 >= Q) i1 -= Q;                            // замыкание по кругу
        const double frac = u - static_cast<double>(i0);
        out[j] = c[i0] + (c[i1] - c[i0]) * frac;
    }
    c.swap(out);
}

// Предвычисленный план: бит-обратные индексы и твиддлы для заданного n.
struct Plan {
    size_t n;
    std::vector<size_t> revers_idx;
    std::vector<cd> twiddles;

    explicit Plan(size_t N) noexcept : n(N){
        [[assume((n & (n - 1)) == 0)]];
        [[assume(n >= 2)]];
        revers_idx.resize(n, 0);
        twiddles.resize(n / 2);

        size_t* __restrict rev = revers_idx.data();
        for (size_t i = 1, j = 0; i < n; ++i) {
            size_t bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            rev[i] = j;
        }

        cd* __restrict tw = twiddles.data();
        const double ang_step = 2.0 * std::numbers::pi / static_cast<double>(n);
        for (size_t i = 0; i < n / 2; ++i) {
            double ang = ang_step * static_cast<double>(i);
            tw[i] = {std::cos(ang), -std::sin(ang)};
        }
    }
};

// Итеративный Кули-Тьюки O(n log n). arr.size() должен совпадать с plan.n.
inline void transform(std::vector<cd>& arr, const Plan& plan) noexcept {
    const size_t n = plan.n;
    if (n <= 1) return;
    [[assume((n & (n - 1)) == 0)]];
    [[assume(n >= 2)]];

    cd* __restrict a = arr.data();
    const size_t* __restrict rev = plan.revers_idx.data();
    const cd* __restrict tw = plan.twiddles.data();

    // Бит-реверсная перестановка
    for (size_t i = 0; i < n; ++i) if (i < rev[i]) std::swap(a[i], a[rev[i]]);

    // Вырожденный случай - бабочки длины 2
    for (size_t i = 0; i < n; i += 2) {
        cd u = a[i];
        cd v = a[i + 1];
        a[i] = u + v; a[i + 1] = u - v;
    }

    // Вырожденный случай - бабочки длины 4
    for (size_t i = 0; i < n; i += 4) {
        cd u1 = a[i];
        cd u2 = a[i + 1];
        cd v1 = a[i + 2];
        cd v2 = {a[i + 3].imag(), -a[i + 3].real()};
        a[i] = u1 + v1; a[i + 1] = u2 + v2;
        a[i + 2] = u1 - v1; a[i + 3] = u2 - v2;
    }

    // Бабочки по каждой паре с увеличивающейся длиной 8 -> 16 -> 32...
    for (size_t len = 8; len <= n; len <<= 1) {
        const size_t half = len >> 1;
        const size_t step = n / len;
        for (size_t i = 0; i < n; i += len) {
            cd u = a[i], v = a[i + half];
            a[i] = u + v;
            a[i + half] = u - v;
            #pragma clang loop vectorize(enable) interleave(enable)
            for (size_t j = 1; j < half; ++j) {
                cd w = tw[j * step];
                cd u = a[i + j], v = a[i + j + half] * w;
                a[i + j] = u + v;
                a[i + j + half] = u - v;
            }
        }
    }
}

// Удобная обёртка для разового преобразования (строит план на месте).
inline void transform(std::vector<cd>& arr) noexcept {
    const size_t n = arr.size();
    if (n <= 1) return;
    transform(arr, Plan(n));
}

} // namespace fft
