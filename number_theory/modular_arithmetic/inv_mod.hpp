#pragma once
#include <bit>
#include <cstdint>
#include <immintrin.h>

// ax≡1 (mod m)
// m - нечётное!
inline uint64_t inv_mod(uint64_t a, uint64_t m) noexcept {
    if (m <= 1) return 0;
    
    uint64_t u = a % m;
    if (u == 0) return 0; // Если a кратно m, обратного не существует

    // 1. Вычисляем -m^{-1} mod 2^64 с помощью метода Ньютона-Рафсона.
    // Это полностью разворачивается компилятором (без циклов, только imul).
    uint64_t inv = m;
    inv *= 2 - m * inv;
    inv *= 2 - m * inv;
    inv *= 2 - m * inv;
    inv *= 2 - m * inv;
    inv *= 2 - m * inv; // 5 итераций достаточно для 64-битной точности
    uint64_t neg_inv_m = -inv;

    uint64_t v = m;
    uint64_t x1 = 1, x2 = 0;

    // 2. Предварительный шаг: если u изначально четное, сдвигаем его.
    uint64_t s = _tzcnt_u64(u);
    if (s) {
        u >>= s;
        // Трюк Монтгомери: вычисляем сдвиг на s бит без цикла
        uint64_t q = _bzhi_u64(x1 * neg_inv_m, s); 
        x1 = (uint64_t)(((unsigned __int128)x1 + (unsigned __int128)q * m) >> s);
    }

    // 3. Основной цикл: u и v всегда нечетные на входе в итерацию
    while (u != v) {
        // Устраняем ветвление: используем булеву переменную
        // Компилятор применит инструкции CMP и CMOV
        bool u_gt_v = (u > v);
        
        // Считаем разности для обоих случаев
        uint64_t diff = u_gt_v ? (u - v) : (v - u);

        uint64_t x_diff1 = x1 - x2;
        if (x1 < x2) x_diff1 += m; // Компилируется в branchless (sub + sbb + add / cmov)

        uint64_t x_diff2 = x2 - x1;
        if (x2 < x1) x_diff2 += m;

        uint64_t x_diff = u_gt_v ? x_diff1 : x_diff2;

        // 4. Трюк Монтгомери для разности (избавляемся от цикла while(s--))
        s = _tzcnt_u64(diff);
        diff >>= s;
        
        // _bzhi_u64 обнуляет все биты старше s. Быстрая инструкция из BMI2.
        uint64_t q = _bzhi_u64(x_diff * neg_inv_m, s);
        
        // icpx транслирует это в пару инструкций (mul, add, adc, shrd). Это мега-быстро.
        uint64_t x_shifted = (uint64_t)(((unsigned __int128)x_diff + (unsigned __int128)q * m) >> s);

        // 5. Обновляем u, v, x1, x2 безусловно через CMOV
        if (u_gt_v) {
            u = diff;
            x1 = x_shifted;
        } else {
            v = diff;
            x2 = x_shifted;
        }
    }

    return (v == 1) ? x2 : 0;
}