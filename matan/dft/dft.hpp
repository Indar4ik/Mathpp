#pragma once
#include <cmath>
#include <complex>
#include <cstddef>
// #include <mkl.h>
#include <numbers>
#include <utility>
#include <vector>

// Рекурсивный Кули-Тьюки O(nlogn)
// void fft(std::vector<std::complex<double>>& arr){
//     const size_t n = arr.size();
//     if (n <= 1) return;

//     std::vector<std::complex<double>> even(n / 2), odd(n / 2);
    
//     #pragma clang loop vectorize(enable)
//     for (size_t i = 0; i < n / 2; ++i){
//         even[i] = arr[2 * i];
//         odd[i] = arr[2 * i + 1];
//     }

//     fft(even);
//     fft(odd);

//     const double twiddle_step = -2.0 * std::numbers::pi / n; 
//     for (size_t i = 0; i < n / 2; ++i){
//         const double twiddle = twiddle_step * i;
//         const std::complex<double> w(std::cos(twiddle), std::sin(twiddle));

//         const std::complex<double> t = odd[i] * w;

//         arr[i] = even[i] + t;
//         arr[i + n / 2] = even[i] - t;
//     }
// }

// Итеративный Кули-Тьюки O(nlogn)
inline void fft_iter(std::vector<std::complex<double>>& arr){
    const size_t n = arr.size();
    if (n <= 1) return;

    // Бит-реверсная перестановка
    for (size_t i = 1, j = 0; i < n; ++i){
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(arr[i], arr[j]);
    }

    std::vector<std::complex<double>> twiddle_table(n / 2);
    for (size_t i = 0; i < n / 2; ++i){
        twiddle_table[i] = std::polar(1.0, -2.0 * std::numbers::pi * i / n);
    }

    // Бабочки по каждой паре с увеличивающейся длинной 2 -> 4 -> 8...
    for (size_t len = 2; len <= n; len <<= 1){
        for (size_t i = 0; i < n; i += len){
            std::complex<double> w(1.0, 0.0);
            for (size_t j = 0; j < len / 2; ++j){
                std::complex<double> w = twiddle_table[j * (n / len)];
                std::complex<double> u = arr[i + j], v = arr[i + j + len / 2] * w;
                arr[i + j] = u + v;
                arr[i + j + len / 2] = u - v;
            }
        }
    }
}

// DFT в лоб за O(n^2)
inline std::vector<std::complex<double>> dft(std::vector<std::complex<double>>& arr){
    const size_t n = arr.size();
    if (n <= 1) return arr;

    std::vector<std::complex<double>> ft(n);
    const double ang_len = 2 * std::numbers::pi / n;
    double ang = 0.0;
    for (size_t i = 0; i < n; ++i){
        std::complex<double> wlen(std::cos(ang), -std::sin(ang));
        std::complex<double> s(0.0, 0.0), w(1.0, 0.0);
        for (size_t j = 0; j < n; ++j){
            s += arr[j] * w; 
            w *= wlen;
        }
        ft[i] = s;
        ang += ang_len;
    }

    return ft;
}