#pragma once
#include <omp.h>
#include <utility>

#include "../../number_theory/fctr.hpp"

namespace ntt {

using namespace hpc;

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
inline constexpr uint64_t NTT_P62 = 4179340454199820289ull; // 29*2^57 + 1, n <= 2^57

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
    if ((p - 1ull) % static_cast<uint64_t>(n) != 0ull) return 0ull;

    const uint64_t minus_one = p - M.one;
    const uint64_t e = (p - 1ull) / static_cast<uint64_t>(n);
    for (uint64_t r = 2ull; r < 1000ull; ++r) {
        const uint64_t g = M.pow(M.to(r), e);
        if (M.pow(g, n / 2) == minus_one) return g;
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
            for (; j & bit; bit >>= 1) j ^= bit;
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
        for (size_t i = K; i < H; ++i) tw[i] = M.mul(tw[i - K], wk); // здесь wk == w^K

        // w^{-i} = w^{n-i} = w^{n/2} * w^{n/2-i} = -w^{n/2-i}, так как w^{n/2} =
        // -1. Обратная таблица выходит одним вычитанием на элемент, без второго
        // прохода умножений. Все твиддлы - обратимые элементы, поэтому нулей здесь
        // не бывает.
        uint64_t *__restrict twi = twiddles_inv.data();
        twi[0] = M.one;
        for (size_t i = 1; i < H; ++i) twi[i] = mod - tw[H - i];

        n_inv = M.pow(M.to(static_cast<uint64_t>(n)), mod - 2ull); // n^{-1} в форме
    }
};

// Порог, до которого выгоден Шоуп. Его таблица твиддлов вдвое толще (пары
// (w, w') - 16 байт против 8), и как только она перестаёт помещаться в L2,
// выигрыш в арифметике съедается лишним трафиком. Замерено на i5-13400F
// (L2 = 1.25 МБ на P-ядро): до 2^16 Шоуп даёт ~1.20x к строгому, на 2^18 уже
// 0.87x, на 2^22 - 0.64x. Выше порога берём ленивый Montgomery.
inline constexpr size_t NTT_SHOUP_MAX = size_t{1} << 16;

// Таблицы для умножения Шоупа: пары (w, w') подряд, чтобы обе половины пары
// приходили одной кэш-линией. w здесь ОБЫЧНЫЙ вычет (не Montgomery-форма),
// w' = floor(w*2^64/p). Держим отдельно от плана: Шоуп стоит вдвое больше памяти
// на твиддлы, и платить за это должен только тот, кто им пользуется.
struct ShoupTables {
    std::vector<uint64_t> fwd;  // [w0, w0', w1, w1', ...]
    std::vector<uint64_t> inv;

    explicit ShoupTables(const NTTPlan& plan) {
        const size_t H = plan.n / 2;
        const mont& M = plan.M;
        fwd.resize(2 * H);
        inv.resize(2 * H);
        for (size_t i = 0; i < H; ++i) {
            const uint64_t wm = plan.twiddles[i];       // Montgomery-форма
            fwd[2 * i] = M.from(wm);                    // обычный вычет
            fwd[2 * i + 1] = M.shoup(wm);               // floor(w*2^64/p), без деления
            const uint64_t wim = plan.twiddles_inv[i];
            inv[2 * i] = M.from(wim);
            inv[2 * i + 1] = M.shoup(wim);
        }
    }
};

namespace detail {

// --- Ленивые бабочки -------------------------------------------------------
// Схема Харви: значения гуляют в [0, 4p), полностью приводим один раз в конце.
// Экономия в том, что modular add/sub (каждый - вычитание + cmov + сложение)
// вырождаются в голое сложение и в "вычесть и добавить 2p", а REDC теряет
// финальное условное вычитание. Требуется p < 2^62: тогда вход бабочки < 4p
// даёт произведение < 4p^2 < p*2^64, что укладывается в контракт REDC.
//
// Инвариант по стадиям:
//   u := a[i+j], приводим до [0, 2p) одним сравнением;
//   v := ленивое умножение, само по себе уже в [0, 2p);
//   вверх u + v            -> [0, 4p)
//   вниз  u - v + 2p       -> [0, 4p)

// Финальное приведение [0, 4p) -> [0, p): два безусловных вычитания.
inline void reduce_full(uint64_t* __restrict a, size_t n, uint64_t p) noexcept {
    const uint64_t p2 = 2 * p;
    for (size_t i = 0; i < n; ++i) {
        uint64_t x = a[i];
        if (x >= p2) x -= p2;
        if (x >= p) x -= p;
        a[i] = x;
    }
}

// Ленивое NTT на Montgomery. tw - твиддлы в Montgomery-форме (как в NTTPlan).
inline void ntt_lazy_impl(uint64_t* __restrict a, const NTTPlan& plan, const uint64_t* __restrict tw) noexcept {
    const size_t n = plan.n;
    [[assume((n & (n - 1)) == 0)]];
    [[assume(n >= 2)]];
    const mont M = plan.M;
    const uint64_t p = M.m, p2 = 2 * p;
    const size_t* __restrict rev = plan.revers_idx.data();

    for (size_t i = 0; i < n; ++i) if (i < rev[i]) std::swap(a[i], a[rev[i]]);

    // Бабочки длины 2: вход ещё в [0, p), поэтому приводить нечего, а выход в [0, 2p).
    for (size_t i = 0; i < n; i += 2) {
        const uint64_t u = a[i], v = a[i + 1];
        a[i] = u + v;
        a[i + 1] = u - v + p;
    }

    for (size_t len = 4; len <= n; len <<= 1) {
        const size_t half = len >> 1;
        const size_t step = n / len;
        for (size_t i = 0; i < n; i += len) {
            // j = 0: твиддл равен 1, но привести надо обе половины, иначе
            // сумма двух значений из [0, 4p) вылезет за инвариант.
            uint64_t u = a[i], v = a[i + half];
            if (u >= p2) u -= p2;
            if (v >= p2) v -= p2;
            a[i] = u + v;
            a[i + half] = u - v + p2;
            // AVX2 не умеет 64x64->128, и векторизатор собирает такое умножение из
            // 32-битных кусков, гоняя значения между GPR и векторными регистрами.
            // Для модулярной арифметики это чистая потеря - запрещаем явно.
            #pragma clang loop vectorize(disable)
            for (size_t j = 1; j < half; ++j) {
                uint64_t u = a[i + j];
                if (u >= p2) u -= p2;
                const uint64_t v = M.mul_lazy(a[i + j + half], tw[j * step]);
                a[i + j + half] = u - v + p2;
            }
        }
    }
    reduce_full(a, n, p);
}

// Ленивое NTT на умножении Шоупа. t - таблица пар (w, w') из ShoupTables.
inline void ntt_shoup_impl(uint64_t* __restrict a, const NTTPlan& plan, const uint64_t* __restrict t) noexcept {
    const size_t n = plan.n;
    [[assume((n & (n - 1)) == 0)]];
    [[assume(n >= 2)]];
    const uint64_t p = plan.M.m, p2 = 2 * p;
    const size_t* __restrict rev = plan.revers_idx.data();

    for (size_t i = 0; i < n; ++i) if (i < rev[i]) std::swap(a[i], a[rev[i]]);

    for (size_t i = 0; i < n; i += 2) {
        const uint64_t u = a[i], v = a[i + 1];
        a[i] = u + v;
        a[i + 1] = u - v + p;
    }

    for (size_t len = 4; len <= n; len <<= 1) {
        const size_t half = len >> 1;
        const size_t step = n / len;
        for (size_t i = 0; i < n; i += len) {
            uint64_t u = a[i], v = a[i + half];
            if (u >= p2) u -= p2;
            if (v >= p2) v -= p2;
            a[i] = u + v;
            a[i + half] = u - v + p2;
            #pragma clang loop vectorize(disable)
            for (size_t j = 1; j < half; ++j) {
                uint64_t u = a[i + j];
                if (u >= p2) u -= p2;
                const size_t k = 2 * (j * step);
                const uint64_t v = mont::shoup_mul_lazy(a[i + j + half], t[k], t[k + 1], p);
                a[i + j] = u + v;
                a[i + j + half] = u - v + p2;
            }
        }
    }
    reduce_full(a, n, p);
}

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
    for (size_t i = 0; i < n; ++i) if (i < rev[i]) std::swap(a[i], a[rev[i]]);

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
            #pragma clang loop vectorize(disable)
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

// Ленивые версии: значения гуляют в [0, 4p), приведение одно на весь проход.
// Требуют p < 2^62; контракт снаружи тот же - вход и выход вычеты в [0, p).
inline void ntt_lazy(std::vector<uint64_t> &arr, const NTTPlan &plan) noexcept {
    if (plan.n <= 1) return;
    detail::ntt_lazy_impl(arr.data(), plan, plan.twiddles.data());
}

inline void intt_lazy(std::vector<uint64_t> &arr, const NTTPlan &plan) noexcept {
    if (plan.n <= 1) return;
    detail::ntt_lazy_impl(arr.data(), plan, plan.twiddles_inv.data());

    const mont M = plan.M;
    const uint64_t n_inv = plan.n_inv;
    for (uint64_t &x : arr) x = M.mul(x, n_inv);
}

// То же, но умножение на твиддл - по Шоупу (нужны ShoupTables).
inline void ntt_shoup(std::vector<uint64_t> &arr, const NTTPlan &plan, const ShoupTables &st) noexcept {
    if (plan.n <= 1) return;
    detail::ntt_shoup_impl(arr.data(), plan, st.fwd.data());
}

inline void intt_shoup(std::vector<uint64_t> &arr, const NTTPlan &plan, const ShoupTables &st) noexcept {
    if (plan.n <= 1) return;
    detail::ntt_shoup_impl(arr.data(), plan, st.inv.data());

    const mont M = plan.M;
    const uint64_t n_inv = plan.n_inv;
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

// Точная свёртка целочисленных последовательностей.
// Входные коэффициенты должны лежать в [0, p), а истинные коэффициенты
// произведения - не превосходить p - 1 (для входов < C и длины перекрытия L это
// значит L*C^2 < p), иначе результат завернётся по модулю и точность будет
// потеряна молча.
inline std::vector<uint64_t> polynom_mul(std::vector<uint64_t> a, std::vector<uint64_t> b, uint64_t mod = NTT_P62) {
    if (a.empty() || b.empty()) return {};

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

    // До NTT_SHOUP_MAX таблица пар (w, w') ещё держится в кэше и Шоуп выигрывает
    // ~1.20x; выше - платим за неё трафиком, и лучше ленивый Montgomery.
    if (n <= NTT_SHOUP_MAX) {
        const ShoupTables st(plan);
        ntt_shoup(a, plan, st);
        ntt_shoup(b, plan, st);
        for (size_t i = 0; i < n; ++i) a[i] = M.mul(a[i], b[i]);
        intt_shoup(a, plan, st);
    } else {
        ntt_lazy(a, plan);
        ntt_lazy(b, plan);
        for (size_t i = 0; i < n; ++i) a[i] = M.mul(a[i], b[i]);
        intt_lazy(a, plan);
    }

    a.resize(k);
    return a;
}

} // namespace ntt