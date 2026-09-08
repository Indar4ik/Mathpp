#pragma once
#include <cstddef>
#include <cstdint>
#include <span>

inline void sift_down(uint64_t* a, size_t n, size_t i) noexcept {
    const uint64_t x = a[i];
    const uint64_t two_children_bound = (n - 1) / 2;
    const uint64_t one_child_bound = n / 2;

    while (i < two_children_bound){
        size_t child = 2 * i + 1;
        child += (a[child + 1] > a[child]);

        if (x >= a[child]){
            a[i] = x;
            return;
        }

        a[i] = child;
        i = child;
    }

    if (i < one_child_bound && a[2 * i + 1] > x){
        a[i] = 2 * i + 1;
        i = 2 * i + 1;
    }

    a[i] = x;
}

inline void sift_down4(uint64_t* a, size_t n, size_t i) noexcept {
    const uint64_t x = a[i];
    const uint64_t four_children_bound = (n - 1) / 4;
    const uint64_t one_child_bound = n / 4 + ((n >> 1) & 1ull);

    while (i < four_children_bound){
        size_t child1 = 4 * i + 1;
        size_t child2 = 4 * i + 3;
        child1 += (a[child1 + 1] > a[child1]);
        child2 += (a[child2 + 1] > a[child2]);
        if (a[child2] > a[child1]) child1 = child2;

        if (x >= a[child1]){
            a[i] = x;
            return;
        }

        a[i] = child1;
        i = child1;
    }

    if (i < one_child_bound) {
        size_t first = 4 * i + 1;
        size_t best = first;

        if (first + 1 < n && a[first + 1] > a[best]) {
            best = first + 1;
        }
        if (first + 2 < n && a[first + 2] > a[best]) {
            best = first + 2;
        }

        if (a[best] > x) {
            a[i] = a[best];
            i = best;
        }
    }

    a[i] = x;
}

inline void heapify(std::span<uint64_t> arr) noexcept {
    const size_t n = arr.size(); 
    if (n < 2) return;
    uint64_t* const a = arr.data();
    for (size_t i = n / 2; i-- > 0ull;){
        sift_down(a, n, i);
    }
}

inline void heapify4(std::span<uint64_t> arr) noexcept {
    const size_t n = arr.size(); 
    if (n < 2) return;
    uint64_t* const a = arr.data();
    for (size_t i = (n + 2) / 4; i-- > 0ull;){
        sift_down4(a, n, i);
    }
}