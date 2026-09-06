#pragma once
// #include <vector>
// #include <algorithm>
#include <span>
#include <cstdint>
#include <cstddef>

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

inline void heapify(std::span<uint64_t> arr) noexcept {
    if (arr.size() < 2) return;
    uint64_t* a = arr.data();
    const size_t n = arr.size(); 
    for (size_t i = n / 2; i-- > 0ull;){
        sift_down(a, n, i);
    }
}