#pragma once
#include "../../number_theory/factorization/fctr.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <omp.h>
#include <utility>
#include <vector>

namespace fft {

using cd = std::complex<double>;

inline void zero_pad(std::vector<cd> &arr) noexcept {
    if (arr.empty()) {
        arr.resize(1, cd(0.0, 0.0));
        return;
    }
    size_t n = arr.size();
    [[assume(n > 0)]];

    if ((n & (n - 1)) == 0)
        return;

    size_t new_size = 1ULL << (64 - std::countl_zero(n - 1));
    arr.resize(new_size, cd(0.0, 0.0));
}

// Передискретизация замкнутого контура до степени двойки без искажения формы.
// В отличие от zero_pad (который дописывает нули и заставляет эпициклы
// танцевать вокруг нуля), здесь линейно интерполируется та же кривая в N =
// bit_ceil точек, равномерно по параметру. Кривая та же, меняется только число
// отсчётов, поэтому бин i по-прежнему отвечает гармонике i.
inline void resample_pow2(std::vector<cd> &c, size_t &n) noexcept {
    const size_t P = c.size();
    if (P < 2)
        return;

    // Если контур замкнут (последняя точка совпадает с первой), отбрасываем
    // дубликат: период образуют Q различных вершин c[0..Q-1] с переходом c[Q-1]
    // -> c[0].
    const size_t Q = (c.front() == c.back()) ? P - 1 : P;
    const size_t N = std::bit_ceil(Q);
    n = N;
    if (N == Q) {
        c.resize(Q);
        return;
    }

    std::vector<cd> out(N);
    const double scale = static_cast<double>(Q) / static_cast<double>(N);
    for (size_t j = 0; j < N; ++j) {
        const double u = static_cast<double>(j) * scale; // позиция в [0, Q)
        const size_t i0 = static_cast<size_t>(u);
        size_t i1 = i0 + 1;
        if (i1 >= Q)
            i1 -= Q; // замыкание по кругу
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

    explicit Plan(size_t N) noexcept : n(N) {
        [[assume((n & (n - 1)) == 0)]];
        [[assume(n >= 2)]];
        revers_idx.resize(n, 0);
        twiddles.resize(n / 2);

        size_t *__restrict rev = revers_idx.data();
        for (size_t i = 1, j = 0; i < n; ++i) {
            size_t bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            rev[i] = j;
        }

        cd *__restrict tw = twiddles.data();
        const double ang_step = 2.0 * std::numbers::pi / static_cast<double>(n);
        for (size_t i = 0; i < n / 2; ++i) {
            double ang = ang_step * static_cast<double>(i);
            tw[i] = {std::cos(ang), -std::sin(ang)};
        }
    }
};

// Итеративный Кули-Тьюки O(n log n). arr.size() должен совпадать с plan.n.
inline void fft(std::vector<cd> &arr, const Plan &plan) noexcept {
    const size_t n = plan.n;
    if (n <= 1)
        return;
    [[assume((n & (n - 1)) == 0)]];
    [[assume(n >= 2)]];

    cd *__restrict a = arr.data();
    const size_t *__restrict rev = plan.revers_idx.data();
    const cd *__restrict tw = plan.twiddles.data();

    // Бит-реверсная перестановка
    for (size_t i = 0; i < n; ++i)
        if (i < rev[i])
            std::swap(a[i], a[rev[i]]);

    // Вырожденный случай - бабочки длины 2
    for (size_t i = 0; i < n; i += 2) {
        cd u = a[i];
        cd v = a[i + 1];
        a[i] = u + v;
        a[i + 1] = u - v;
    }

    // Вырожденный случай - бабочки длины 4
    for (size_t i = 0; i < n; i += 4) {
        cd u1 = a[i];
        cd u2 = a[i + 1];
        cd v1 = a[i + 2];
        cd v2 = {a[i + 3].imag(), -a[i + 3].real()};
        a[i] = u1 + v1;
        a[i + 1] = u2 + v2;
        a[i + 2] = u1 - v1;
        a[i + 3] = u2 - v2;
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

// Обратный итеративный Кули-Тьюки O(n log n). arr.size() должен совпадать с
// plan.n.
inline void ifft(std::vector<cd> &arr, const Plan &plan) noexcept {
    const size_t n = plan.n;
    if (n <= 1)
        return;
    [[assume((n & (n - 1)) == 0)]];
    [[assume(n >= 2)]];

    cd *__restrict a = arr.data();
    const size_t *__restrict rev = plan.revers_idx.data();
    const cd *__restrict tw = plan.twiddles.data();

    // Бит-реверсная перестановка
    for (size_t i = 0; i < n; ++i)
        if (i < rev[i])
            std::swap(a[i], a[rev[i]]);

    // Вырожденный случай - бабочки длины 2
    for (size_t i = 0; i < n; i += 2) {
        cd u = a[i];
        cd v = a[i + 1];
        a[i] = u + v;
        a[i + 1] = u - v;
    }

    // Вырожденный случай - бабочки длины 4
    for (size_t i = 0; i < n; i += 4) {
        cd u1 = a[i];
        cd u2 = a[i + 1];
        cd v1 = a[i + 2];
        cd v2 = {-a[i + 3].imag(), a[i + 3].real()};
        a[i] = u1 + v1;
        a[i + 1] = u2 + v2;
        a[i + 2] = u1 - v1;
        a[i + 3] = u2 - v2;
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
                cd w = std::conj(tw[j * step]);
                cd u = a[i + j], v = a[i + j + half] * w;
                a[i + j] = u + v;
                a[i + j + half] = u - v;
            }
        }
    }

    // Единственная нормировка на 1/n вместо умножения на 0.5 на каждом этапе.
    const double inv_n = 1.0 / static_cast<double>(n);
    for (cd &z : arr)
        z *= inv_n;
}

// Многопоточный прямой БПФ (OpenMP). Один parallel-регион на весь вызов, чтобы
// не платить fork/join на каждой из log2(n) стадий. Ранние стадии (много
// блоков) распараллеливаются по блокам, поздние (мало блоков, но большой half)
// - по j. Для мелких n накладные не окупаются -> делегируем в серийную fft().
inline void fft_mt(std::vector<cd> &arr, const Plan &plan) noexcept {
    const size_t n = plan.n;
    if (n < (1u << 15)) {
        fft(arr, plan);
        return;
    }
    [[assume((n & (n - 1)) == 0)]];
    [[assume(n >= 2)]];

    cd *__restrict a = arr.data();
    const size_t *__restrict rev = plan.revers_idx.data();
    const cd *__restrict tw = plan.twiddles.data();

#ifdef _OPENMP
    const size_t nthreads = static_cast<size_t>(omp_get_max_threads());
#else
    const size_t nthreads = 1;
#endif

#pragma omp parallel
    {
// Бит-реверсная перестановка: пары {i, rev[i]} не пересекаются -> без гонок.
        #pragma omp for schedule(static)
        for (size_t i = 0; i < n; ++i)
            if (i < rev[i])
                std::swap(a[i], a[rev[i]]);

// Вырожденный случай - бабочки длины 2
        #pragma omp for schedule(static)
        for (size_t i = 0; i < n; i += 2) {
            cd u = a[i];
            cd v = a[i + 1];
            a[i] = u + v;
            a[i + 1] = u - v;
        }

// Вырожденный случай - бабочки длины 4
        #pragma omp for schedule(static)
        for (size_t i = 0; i < n; i += 4) {
            cd u1 = a[i];
            cd u2 = a[i + 1];
            cd v1 = a[i + 2];
            cd v2 = {a[i + 3].imag(), -a[i + 3].real()};
            a[i] = u1 + v1;
            a[i + 1] = u2 + v2;
            a[i + 2] = u1 - v1;
            a[i + 3] = u2 - v2;
        }

        // Бабочки по каждой паре с увеличивающейся длиной 8 -> 16 -> 32...
        for (size_t len = 8; len <= n; len <<= 1) {
            const size_t half = len >> 1;
            const size_t step = n / len;
            const size_t blocks = n / len;

            if (blocks >= nthreads) {
            // Достаточно блоков -> параллель по внешнему циклу (блокам)
            #pragma omp for schedule(static)
                for (size_t b = 0; b < blocks; ++b) {
                    const size_t i = b * len;
                    #pragma clang loop vectorize(enable) interleave(enable)
                    for (size_t j = 0; j < half; ++j) {
                        cd w = tw[j * step];
                        cd u = a[i + j], v = a[i + j + half] * w;
                        a[i + j] = u + v;
                        a[i + j + half] = u - v;
                    }
                }
            } else {
                // Мало блоков, но большой half -> параллель по j внутри блока
                for (size_t i = 0; i < n; i += len) {
                    #pragma omp for schedule(static)
                    for (size_t j = 0; j < half; ++j) {
                        cd w = tw[j * step];
                        cd u = a[i + j], v = a[i + j + half] * w;
                        a[i + j] = u + v;
                        a[i + j + half] = u - v;
                    }
                }
            }
        }
    }
}

// Удобные обёртки для разового преобразования (строят план на месте).
inline void fft(std::vector<cd> &arr) noexcept {
    const size_t n = arr.size();
    if (n <= 1)
        return;
    fft(arr, Plan(n));
}

inline void ifft(std::vector<cd> &arr) noexcept {
    const size_t n = arr.size();
    if (n <= 1)
        return;
    ifft(arr, Plan(n));
}

inline void fft_mt(std::vector<cd> &arr) noexcept {
    const size_t n = arr.size();
    if (n <= 1)
        return;
    fft_mt(arr, Plan(n));
}

// ===========================================================================
// NTT - та же схема Кули-Тьюки, но над Z_p вместо C.
//
// Про "знак" корня. В Z_p нет направления вращения: если w имеет порядок n, то
// w^{-1} тоже примитивный корень n-й степени, и никакого признака, который
// делал бы один из них аналогом e^{+2pi i/n}, а другой - e^{-2pi i/n}, в поле
// нет. Обе матрицы DFT одинаково законны и взаимно обратны, так что выбор
// прямой - чистое соглашение. Существенна только согласованность: прямое идёт
// по w, обратное - по w^{-1} и делит на n. Для свёртки выбор вообще
// ненаблюдаем.
// ===========================================================================

// NTT-дружественные простые p = c*2^k + 1: длина преобразования ограничена 2^k.
inline constexpr uint64_t NTT_P30 = 469762049ull; // 7*2^26 + 1,    n <= 2^26
inline constexpr uint64_t NTT_P62 =
        4179340454199820289ull; // 29*2^57 + 1, n <= 2^57

// Максимальная длина NTT по модулю p: наибольшая степень двойки, делящая p - 1.
inline constexpr size_t ntt_max_len(uint64_t p) noexcept {
    return size_t{1} << std::countr_zero(p - 1ull);
}

// Корень степени ровно n из единицы, в Montgomery-форме.
// Кандидат r даёт g = r^{(p-1)/n}; порядок g заведомо делит n, а раз n -
// степень двойки, порядок равен n тогда и только тогда, когда g^{n/2} == -1.
// Проверка избавляет от необходимости знать примитивный корень p или
// факторизовать p - 1. Возвращает 0, если n не делит p - 1 (или p не простое).
inline uint64_t ntt_root(const mont &M, size_t n) noexcept {
    [[assume(n >= 2)]];
    const uint64_t p = M.m;
    if ((p - 1ull) % static_cast<uint64_t>(n) != 0ull)
        return 0ull;

    const uint64_t minus_one = p - M.one;
    const uint64_t e = (p - 1ull) / static_cast<uint64_t>(n);
    for (uint64_t r = 2ull; r < 1000ull; ++r) {
        const uint64_t g = M.pow(M.to(r), e);
        if (M.pow(g, n / 2) == minus_one)
            return g;
    }
    return 0ull;
}

// Предвычисленный план: бит-обратные индексы и твиддлы (прямые и обратные) для
// n. Требуется: n - степень двойки, n >= 2, n <= ntt_max_len(mod), mod -
// нечётное простое.
struct NTTPlan {
    size_t n;
    mont M;
    uint64_t n_inv; // n^{-1} в Montgomery-форме (нормировка intt)
    std::vector<size_t> revers_idx;
    std::vector<uint64_t> twiddles;         // w^i,        i < n/2, в Montgomery-форме
    std::vector<uint64_t> twiddles_inv; // w^{-i}, i < n/2, в Montgomery-форме

    explicit NTTPlan(size_t N, uint64_t mod = NTT_P62) noexcept : n(N), M(mod) {
        [[assume((n & (n - 1)) == 0)]];
        [[assume(n >= 2)]];
        revers_idx.resize(n, 0);
        twiddles.resize(n / 2);
        twiddles_inv.resize(n / 2);

        size_t *__restrict rev = revers_idx.data();
        for (size_t i = 1, j = 0; i < n; ++i) {
            size_t bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            rev[i] = j;
        }

        // const uint64_t w = ntt_root(M, n);
        const uint64_t w = M.pow(M.add(M.one, M.add(M.one, M.one)), (mod - 1ull) / n);
        const size_t H = n / 2;
        uint64_t *__restrict tw = twiddles.data();

        // Бегущее произведение, но K независимыми цепочками: tw[i] = tw[i-K] * w^K.
        // Одна цепочка упёрлась бы в латентность REDC (~13 тактов на элемент), а
        // pow(w, i) на каждый элемент - это ~2*log2(i) таких же последовательных
        // умножений. Дистанция зависимости K разводит умножения по конвейеру и
        // упирает цикл в пропускную способность умножителя (~3 такта на элемент).
        // Арифметика точная, поэтому перегруппировка даёт побитово тот же результат
        // (в отличие от cos/sin в вещественном Plan, где рекуррентность копила бы
        // дрейф).
        constexpr size_t K = 8;
        const size_t head = std::min(K, H);
        uint64_t wk = M.one;
        for (size_t i = 0; i < head; ++i) {
            tw[i] = wk;
            wk = M.mul(wk, w);
        }
        for (size_t i = K; i < H; ++i)
            tw[i] = M.mul(tw[i - K], wk); // здесь wk == w^K

        // w^{-i} = w^{n-i} = w^{n/2} * w^{n/2-i} = -w^{n/2-i}, так как w^{n/2} =
        // -1. Обратная таблица выходит одним вычитанием на элемент, без второго
        // прохода умножений. Все твиддлы - обратимые элементы, поэтому нулей здесь
        // не бывает.
        uint64_t *__restrict twi = twiddles_inv.data();
        twi[0] = M.one;
        for (size_t i = 1; i < H; ++i)
            twi[i] = mod - tw[H - i];

        n_inv = M.pow(M.to(static_cast<uint64_t>(n)), mod - 2ull); // n^{-1} в форме
    }
};

namespace detail {

// Тело преобразования: направление задаётся только таблицей твиддлов.
// Данные - обычные вычеты в [0, p), твиддлы - в Montgomery-форме. Тогда
// M.mul(данные, твиддл) = данные * твиддл снова в обычном виде: множители R и
// R^{-1} сокращаются, и переводить массив в форму и обратно не нужно.
inline void ntt_impl(uint64_t *__restrict a, const NTTPlan &plan, const uint64_t *__restrict tw) noexcept {
    const size_t n = plan.n;
    [[assume((n & (n - 1)) == 0)]];
    [[assume(n >= 2)]];
    const mont M = plan.M;
    const size_t *__restrict rev = plan.revers_idx.data();

    // Бит-реверсная перестановка
    for (size_t i = 0; i < n; ++i)
        if (i < rev[i])
            std::swap(a[i], a[rev[i]]);

    // Вырожденный случай - бабочки длины 2 (твиддл равен 1, умножения нет)
    for (size_t i = 0; i < n; i += 2) {
        uint64_t u = a[i], v = a[i + 1];
        a[i] = M.add(u, v);
        a[i + 1] = M.sub(u, v);
    }

    for (size_t len = 4; len <= n; len <<= 1) {
        const size_t half = len >> 1;
        const size_t step = n / len;
        for (size_t i = 0; i < n; i += len) {
            uint64_t u = a[i], v = a[i + half];
            a[i] = M.add(u, v);
            a[i + half] = M.sub(u, v);
            for (size_t j = 1; j < half; ++j) {
                const uint64_t w = tw[j * step];
                const uint64_t u = a[i + j], v = M.mul(a[i + j + half], w);
                a[i + j] = M.add(u, v);
                a[i + j + half] = M.sub(u, v);
            }
        }
    }
}

} // namespace detail

// Прямое NTT. arr.size() должен совпадать с plan.n, значения - вычеты в [0, p).
inline void ntt(std::vector<uint64_t> &arr, const NTTPlan &plan) noexcept {
    if (plan.n <= 1) return;
    detail::ntt_impl(arr.data(), plan, plan.twiddles.data());
}

// Обратное NTT: та же схема по обратным твиддлам плюс единственная нормировка
// на 1/n.
inline void intt(std::vector<uint64_t> &arr, const NTTPlan &plan) noexcept {
    if (plan.n <= 1) return;
    detail::ntt_impl(arr.data(), plan, plan.twiddles_inv.data());

    const mont M = plan.M;
    const uint64_t n_inv = plan.n_inv; // в Montgomery-форме -> mul снимает лишний R
    for (uint64_t &x : arr) x = M.mul(x, n_inv);
}

// Удобные обёртки для разового преобразования (строят план на месте).
inline void ntt(std::vector<uint64_t> &arr, uint64_t mod = NTT_P62) noexcept {
    if (arr.size() <= 1) return;
    ntt(arr, NTTPlan(arr.size(), mod));
}

inline void intt(std::vector<uint64_t> &arr, uint64_t mod = NTT_P62) noexcept {
    if (arr.size() <= 1) return;
    intt(arr, NTTPlan(arr.size(), mod));
}

inline std::vector<cd> polynom_mul(std::vector<cd> a, std::vector<cd> b) {
    if (a.empty() || b.empty())
        return {};

    const size_t k = a.size() + b.size() - 1ull;
    const size_t n = std::bit_ceil(k);

    a.resize(n, cd(0.0, 0.0));
    b.resize(n, cd(0.0, 0.0));

    std::vector<cd> c(n);
    Plan plan(n);
    fft(a, plan);
    fft(b, plan);
    for (size_t i = 0; i < n; ++i)
        c[i] = a[i] * b[i];
    ifft(c, plan);

    c.resize(k);
    return c;
}

inline std::vector<double> polynom_mul(std::vector<double> a, std::vector<double> b) {
    if (a.empty() || b.empty())
        return {};

    const size_t n = a.size(), m = b.size();
    const size_t k = n + m - 1ull;
    const size_t new_size = std::bit_ceil(k);

    std::vector<cd> ab(new_size, cd(0.0, 0.0));
    std::vector<double> c(k);
    Plan plan(new_size);

    // Упаковка двух вещественных многочленов в один комплексный: P = a + i*b.
    for (size_t i = 0; i < n; ++i)
        ab[i].real(a[i]);
    for (size_t i = 0; i < m; ++i)
        ab[i].imag(b[i]);

    fft(ab, plan);
    // P = a + i*b    =>    P^2 = a^2 - b^2 + 2i*(a*b) -> a*b=Im(P^2)*0.5
    // IFFT(FFT(a)*FFT(b))= a * b свёртка
    for (std::complex<double> &z : ab)
        z = z * z;
    ifft(ab, plan);

    for (size_t i = 0; i < k; ++i)
        c[i] = ab[i].imag() * 0.5;
    return c;
}

// Точная свёртка целочисленных последовательностей.
// Входные коэффициенты должны лежать в [0, p), а истинные коэффициенты
// произведения - не превосходить p - 1 (для входов < C и длины перекрытия L это
// значит L*C^2 < p), иначе результат завернётся по модулю и точность будет
// потеряна молча.
inline std::vector<uint64_t> polynom_mul(std::vector<uint64_t> a, std::vector<uint64_t> b, uint64_t mod = NTT_P62) {
    if (a.empty() || b.empty())
        return {};

    const size_t k = a.size() + b.size() - 1ull;
    const size_t n = std::max<size_t>(2ull, std::bit_ceil(k));

    NTTPlan plan(n, mod);
    const mont M = plan.M;

    // a переводим в Montgomery-форму (преобразование линейно, форма переживает
    // его). Тогда поточечное M.mul(A, B) сокращает лишний R, и результат снова
    // обычные вычеты.
    for (uint64_t &x : a) x = M.to(x);
    a.resize(n, 0ull);
    b.resize(n, 0ull);

    ntt(a, plan);
    ntt(b, plan);
    for (size_t i = 0; i < n; ++i) a[i] = M.mul(a[i], b[i]);
    intt(a, plan);

    a.resize(k);
    return a;
}
} // namespace fft
