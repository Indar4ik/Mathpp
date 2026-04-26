#pragma once
#include <bit>
#include <cstdint>
#include <immintrin.h>

// inline uint64_t add_mod(uint64_t a, uint64_t m){
//     uint64_t res, sum;
//     uint8_t carry = __builtin_add_overflow(a, m, &sum);
//     return (sum >> 1) | (carry ? 0x8000000000000000ULL : 0);
// }

// inline uint64_t inv_mod1(uint64_t a, uint64_t m) noexcept{
//     if (m == 1) return 0;
    
//     uint64_t u = a % m; // u = x1*a+y1*m
//     uint64_t v = m; // v = x2*a+y2*m
//     uint64_t x1 = 1, x2 = 0;
    
//     while (u != 0){
//         int shift = std::countr_zero(u);
//         u >>= shift;
//         for (int i = 0; i < shift; i++){
//             uint64_t mask = -(x1 & 1ULL);
//             x1 = add_mod(x1, m & mask);
//         }

//         shift = std::countr_zero(v);
//         v >>= shift;
//         for (int i = 0; i < shift; i++){
//             uint64_t mask = -(x2 & 1ULL);
//             x2 = add_mod(x2, m & mask);
//         }

//         if (u >= v){
//             u -= v;
//             if (x1 < x2) x1 += m;
//             x1 -= x2;
//         }else{
//             v -= u;
//             if (x2 < x1) x2 += m;
//             x2 -= x1;
//         }
//     }

//     return (v == 1) ? x2 : 0;
// }

// ax≡1 (mod m)
inline uint64_t inv_mod(uint64_t a, uint64_t m) {
    if (m <= 1) return 0;
    uint64_t u = a % m;
    uint64_t v = m;
    uint64_t x1 = 1, x2 = 0;

    while (u != 0) {
        // 1. Быстрый сдвиг u и x1
        uint64_t s = _tzcnt_u64(u);
        u >>= s;
        while (s--) {
            uint64_t mask = -(x1 & 1);
            unsigned long long sum;
            uint8_t carry = _addcarry_u64(0, x1, m & mask, &sum);
            x1 = (sum >> 1) | (carry ? 0x8000000000000000ULL : 0);
        }

        // 2. Быстрый сдвиг v и x2
        s = _tzcnt_u64(v);
        v >>= s;
        while (s--) {
            uint64_t mask = -(x2 & 1);
            unsigned long long sum;
            uint8_t carry = _addcarry_u64(0, x2, m & mask, &sum);
            x2 = (sum >> 1) | (carry ? 0x8000000000000000ULL : 0);
        }

        // 3. Branchless вычитание через тернарные операторы
        // Для компилятора это прямой сигнал генерировать CMOV
        uint64_t next_u_v = (u >= v) ? u - v : v - u;
        
        // Коэффициенты x1, x2
        uint64_t x1_new = x1 - x2;
        if (x1 < x2) x1_new += m;
        
        uint64_t x2_new = x2 - x1;
        if (x2 < x1) x2_new += m;

        if (u >= v) {
            u = next_u_v;
            x1 = x1_new;
        } else {
            v = next_u_v;
            x2 = x2_new;
        }
    }

    return (v == 1) ? x2 : 0;
}

