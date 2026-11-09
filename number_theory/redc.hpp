#pragma once
#include <concepts>

#include "base.hpp"

namespace hpc {

// Montgomery-арифметика по нечётному модулю m. Параметры шаблона - машинное
// слово W (оно же R = 2^bits) и двойное слово D для произведений; вся структура
// 32- и 64-битной версий совпадает с точностью до этих двух типов.
// Все значения внутри хранятся в форме a_bar = a * R mod m; константы считаются
// раз на m. Всё constexpr: при фиксированном модуле контекст строится компилятором.
template <std::unsigned_integral W>
struct montgomery {
    static constexpr int BITS = static_cast<int>(sizeof(W) * 8);
    using D = detail::wider_uint_t<W>;

    W m;      // модуль (нечётный)
    W m_inv;  // m^{-1} mod R
    W r2;     // R^2 mod m   (для перевода в форму)
    W one;    // R   mod m   (Montgomery-форма единицы)

    constexpr explicit montgomery(W mod) noexcept : m(mod) {
        [[assume((m & W(1)) == W(1))]];
        m_inv = inv_mod2<W>(m);                           // m_inv * m == 1 (mod R)
        one = static_cast<W>(-m) % m;                       // R mod m = (R - m mod m) = (-m) mod m
        r2 = static_cast<W>(static_cast<D>(one) * one % m); // (R mod m)^2 mod m = R^2 mod m
    }

    // REDC: по x < m*R возвращает x * R^{-1} mod m в [0, m).
    // Вычитательная форма: t = lo(x) * m^{-1}, тогда t*m == x (mod R), то есть
    // младшие лимбы x и t*m совпадают ТОЧНО -> при вычитании они гасятся без заёма,
    // а частное (x - t*m)/R равно разности старших лимбов и лежит в (-m, m).
    // Переполнения двойного слова здесь не бывает по построению (в отличие от формы
    // x + t*m, где при m ~ R приходится отдельно ловить перенос), поэтому хвост -
    // ровно sub + cmov + add вместо add/adc + cmp/sbb + двух cmov.
    constexpr W redc(D x) const noexcept {
        W t = static_cast<W>(x) * m_inv;
        W th = static_cast<W>((static_cast<D>(t) * m) >> BITS);
        W r;
        bool borrow = __builtin_sub_overflow(static_cast<W>(x >> BITS), th, &r);
        return borrow ? static_cast<W>(r + m) : r;  // ушли ниже нуля -> возвращаем в [0, m)
    }

    // «Ленивый» REDC: без финального условного вычитания, результат в [0, 2m).
    // Здесь нужна именно АДДИТИВНАЯ форма (x + t*m): её частное неотрицательно по
    // построению, тогда как вычитательная даёт результат в (-m, m) и без коррекции
    // не обходится. Отрицательный обратный считаем на месте - это один neg, и он
    // выносится из цикла как инвариант.
    // Контракт: x < m*R, результат в [0, 2m); требуется m < R/2 (иначе 2m не влезет
    // в слово, а сумма x + t*m - в двойное).
    constexpr W redc_lazy(D x) const noexcept {
        W t = static_cast<W>(x) * static_cast<W>(W(0) - m_inv);
        return static_cast<W>((x + static_cast<D>(t) * m) >> BITS);
    }

    constexpr W mul(W a, W b) const noexcept { return redc(static_cast<D>(a) * b); }
    constexpr W mul_lazy(W a, W b) const noexcept { return redc_lazy(static_cast<D>(a) * b); }
    constexpr W to(W a) const noexcept { return mul(a, r2); }   // a -> a*R mod m
    constexpr W from(W a) const noexcept { return redc(a); }    // a*R -> a mod m

    // Предвычисление для умножения Шоупа на ФИКСИРОВАННЫЙ множитель w: w' = floor(w*R/m).
    // Деления не нужно: пусть r = w*R mod m (это ровно Montgomery-форма w). Тогда
    // w'*m = w*R - r, и по модулю R получаем w' == -r * m^{-1}. Так как w < m,
    // частное w' < R, поэтому сравнение по модулю задаёт его точно.
    // Вход - w в Montgomery-форме (то есть уже r), выход - w' для shoup_mul_lazy.
    constexpr W shoup(W w_mont) const noexcept { return static_cast<W>(W(0) - w_mont) * m_inv; }

    // Умножение Шоупа: b * w mod m, где w - ОБЫЧНЫЙ вычет, w_sh = shoup(to(w)).
    // Умножений столько же, сколько у Montgomery (три), но цепочка зависимостей
    // короче: mulhi и mullo стартуют от b одновременно. Результат в [0, 2m).
    // Контракт: w < m, m < R/2, b - любое слово.
    static constexpr W shoup_mul_lazy(W b, W w, W w_sh, W m) noexcept {
        W q = static_cast<W>((static_cast<D>(w_sh) * b) >> BITS);
        return static_cast<W>(w * b - q * m);
    }

    // Сложение и вычитание по модулю: входы в [0, m), выход в [0, m).
    // Montgomery-форма линейна, поэтому одни и те же add/sub годятся и для обычных
    // вычетов, и для значений в форме.
    constexpr W add(W a, W b) const noexcept {
        // a + b - m считаем как a - (m - b): так не переполняем слово даже при m ~ R.
        W d;
        bool borrow = __builtin_sub_overflow(a, static_cast<W>(m - b), &d);
        return borrow ? static_cast<W>(d + m) : d;  // a + b < m -> вернуть a + b
    }

    constexpr W sub(W a, W b) const noexcept {
        W d;
        bool borrow = __builtin_sub_overflow(a, b, &d);
        return borrow ? static_cast<W>(d + m) : d;
    }

    // a^e mod m, вход/выход в Montgomery-форме.
    constexpr W pow(W a, W e) const noexcept {
        W r = one;
        while (e) {
            if (e & W(1)) r = mul(r, a);
            a = mul(a, a);
            e >>= 1;
        }
        return r;
    }

    // Контракт: 4m < R
    constexpr W pow_lazy(W a, W e) const noexcept {
        W r = one;
        while (e) {
            if (e & W(1)) r = mul_lazy(r, a);
            a = mul_lazy(a, a);
            e >>= 1;
        }
        return r;
    }
};

using mont = montgomery<uint64_t>;
using mont32 = montgomery<uint32_t>;

// Общий потолок ленивой арифметики: держим величины в [0, 2m), тогда произведение
// двух таких меньше 4m^2, а контракт redc_lazy требует вход < mR — отсюда 4m <= R.
inline constexpr uint64_t LAZY_MAX = ~0ull / 4;

template <typename W>
struct barrett {
    static constexpr int BITS = static_cast<int>(sizeof(W) * 8);
    using D = detail::wider_uint_t<W>;

    W m;     // любой модуль
    W m_inv; // floor(R / m)

    // Скобки обязательны: '/' связывает сильнее '<<', и без них выходило D(1) << (BITS / m),
    // то есть m_inv = 1 при любом m > BITS. Контракт: m >= 2 (иначе R / m не влезает в W).
    constexpr explicit barrett(W mod) noexcept : m(mod) {
        m_inv = static_cast<W>((D{1} << BITS) / m);
    }

    // REDC: x < R, возвращает x mod m
    constexpr W redc(W x) const noexcept {
        W q = static_cast<W>(((static_cast<D>(x) * m_inv) >> BITS));
        W r_bar = x - m * q;
        W r;
        bool borrow = __builtin_sub_overflow(r_bar, m, &r);
        return borrow ? r + m : r;
    }

    // «Ленивый» REDC: без финального условного вычитания.
    // Контракт: x < R (любое слово), результат в [0, 2m) и при этом не больше x.
    constexpr W redc_lazy(W x) const noexcept {
        W q = static_cast<W>(((static_cast<D>(x) * m_inv) >> BITS));
        return x - m * q;
    }
};

} // namespace hpc