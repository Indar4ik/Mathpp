#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <immintrin.h>
#include <span>
#include <type_traits>

// Просеивание вниз в бинарной куче.
// Дырка едет вниз, x кладётся один раз в конце: так на уровень приходится
// одна запись вместо трёх, как было бы при обмене.
template <typename T>
inline void sift_down(T *const a, size_t n, size_t i) noexcept {
    if (n < 2) return;
    const T x = a[i];
    const size_t two_children_bound = (n - 1) / 2; // узлы, у которых есть оба ребёнка
    const size_t one_child_bound = n / 2;          // узлы, у которых есть хоть один

    while (i < two_children_bound) {
        size_t child = 2 * i + 1;
        child += (a[child + 1] > a[child]);

        if (x >= a[child]) {
            a[i] = x;
            return;
        }

        a[i] = a[child];
        i = child;
    }

    // Неполный ребёнок бывает только у одного узла, и его дети -- листья,
    // поэтому здесь достаточно одного шага
    if (i < one_child_bound && a[2 * i + 1] > x) {
        a[i] = a[2 * i + 1];
        i = 2 * i + 1;
    }

    a[i] = x;
}

// Просеивание вниз в 4-арной куче.
// Уровней вчетверо меньше, чем в бинарной по основанию 2 -- а именно число
// уровней, а не число сравнений, определяет цену: каждый спуск это промах
// кэша, зависящий от предыдущего.
template <typename T>
inline void sift_down4(T *const a, size_t n, size_t i) noexcept {
    if (n < 2) return;
    const T x = a[i];
    const size_t four_children_bound = (n - 1) / 4;
    const size_t one_child_bound = (n + 2) / 4;

    while (i < four_children_bound) {
        const size_t c = 4 * i + 1;

        // Дети узлов c..c+3 лежат подряд по индексам 4c+1..4c+16, то есть в
        // 64 байтах -- одной-двух кэш-линиях. Заказываем их до того, как узнаем,
        // в какого именно ребёнка спустимся: линия успеет приехать. Это здесь
        // самая крупная оптимизация, на 16M элементов она даёт около -25%.
        if (4 * c + 16 < n) {
            __builtin_prefetch(a + 4 * c + 1, 0, 3);
            __builtin_prefetch(a + 4 * c + 9, 0, 3);
        }

        T best;
        size_t child;
        if constexpr (std::is_same_v<T, int32_t>) {
            // Четыре ребёнка -- ровно один 16-байтный вектор. Горизонтальный
            // максимум двумя парами shuffle+max, всё в регистрах. Скалярный
            // турнир вместо этого делает две загрузки по адресу, который сам
            // зависит от сравнения, и эта цепочка стоит дороже.
            const __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i *>(a + c));
            __m128i m = _mm_max_epi32(v, _mm_shuffle_epi32(v, _MM_SHUFFLE(1, 0, 3, 2)));
            m = _mm_max_epi32(m, _mm_shuffle_epi32(m, _MM_SHUFFLE(2, 3, 0, 1)));
            best = _mm_cvtsi128_si32(m);

            if (x >= best) {
                a[i] = x;
                return;
            }
            // Индекс достаём только когда действительно спускаемся
            const unsigned mask = _mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(v, m)));
            child = c + static_cast<size_t>(__builtin_ctz(mask));
        } else {
            size_t c1 = c;
            size_t c2 = c + 2;
            c1 += (a[c1 + 1] > a[c1]);
            c2 += (a[c2 + 1] > a[c2]);
            if (a[c2] > a[c1]) c1 = c2;
            best = a[c1];

            if (x >= best) {
                a[i] = x;
                return;
            }
            child = c1;
        }

        a[i] = best;
        i = child;
    }

    // Узел с неполным набором детей единственный, и его дети -- листья
    if (i < one_child_bound) {
        const size_t first = 4 * i + 1;
        const size_t last = std::min(first + 3, n - 1);
        size_t best = first;
        for (size_t c = first + 1; c <= last; ++c) {
            if (a[c] > a[best]) best = c;
        }
        if (a[best] > x) {
            a[i] = a[best];
            i = best;
        }
    }

    a[i] = x;
}

template <typename T>
inline void heapify(std::span<T> arr) noexcept {
    const size_t n = arr.size();
    if (n < 2) return;
    T *const a = arr.data();
    for (size_t i = n / 2; i-- > 0ull;) {
        sift_down(a, n, i);
    }
}

template <typename T>
inline void heapify4(std::span<T> arr) noexcept {
    const size_t n = arr.size();
    if (n < 2) return;
    T *const a = arr.data();
    for (size_t i = (n + 2) / 4; i-- > 0ull;) {
        sift_down4(a, n, i);
    }
}
