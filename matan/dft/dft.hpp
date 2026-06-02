#pragma once
// #include <cmath>
#include <complex>
#include <cstddef>
// #include <mkl.h>
#include <numbers>
#include <utility>
#include <vector>

inline void prepare_fft(size_t& n, std::vector<size_t>& revers_idx, std::vector<std::complex<double>>& twiddles){
    [[assume((n & (n - 1)) == 0)]];
    [[assume(n >= 2)]];
    
    revers_idx.resize(n, 0);
    twiddles.resize(n / 2);

    for (size_t i = 1, j = 0; i < n; ++i){
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        revers_idx[i] = j;
    }
    
    for (size_t i = 0; i < n / 2; ++i) twiddles[i] = {std::cos(2.0 * std::numbers::pi * i / n), -std::sin(2.0 * std::numbers::pi * i / n)};
}

// Итеративный Кули-Тьюки O(nlogn)
inline void fft(std::vector<std::complex<double>>& arr, const std::vector<size_t>& revers_idx, const std::vector<std::complex<double>>& twiddles) noexcept {
    const size_t n = arr.size();
    if (n <= 1) return;
    [[assume((n & (n - 1)) == 0)]];
    [[assume(n >= 2)]];
    [[assume(revers_idx.size() == n)]];
    [[assume(twiddles.size() == n / 2)]];

    // Бит-реверсная перестановка
    for (size_t i = 0; i < n; ++i) if (i < revers_idx[i]) std::swap(arr[i], arr[revers_idx[i]]);

    // Вырожденный случай - бабочки длины 2
    for (size_t i = 0; i < n; i += 2){
        std::complex<double> u = arr[i];
        std::complex<double> v = arr[i + 1];
        arr[i] = u + v; arr[i + 1] = u - v;
    }

    // Вырожденный случай - бабочки длины 4
    for (size_t i = 0; i < n; i += 4){
        std::complex<double> u1 = arr[i];
        std::complex<double> u2 = arr[i + 1];
        std::complex<double> v1 = arr[i + 2];
        std::complex<double> v2 = {arr[i + 3].imag(), -arr[i + 3].real()};
        arr[i] = u1 + v1; arr[i + 1] = u2 + v2; 
        arr[i + 2] = u1 - v1; arr[i + 3] = u2 - v2; 
    }

    // Бабочки по каждой паре с увеличивающейся длинной 8 -> 16 -> 32...
    for (size_t len = 8; len <= n; len <<= 1){
        for (size_t i = 0; i < n; i += len){
            std::complex<double> u = arr[i], v = arr[i + len / 2];
            arr[i] = u + v;
            arr[i + len / 2] = u - v;
            #pragma clang loop vectorize(enable) interleave(enable)
            for (size_t j = 1; j < len / 2; ++j){
                std::complex<double> w = twiddles[j * (n / len)];
                std::complex<double> u = arr[i + j], v = arr[i + j + len / 2] * w;
                arr[i + j] = u + v;
                arr[i + j + len / 2] = u - v;
            }
        }
    }
}
