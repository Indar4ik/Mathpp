#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <immintrin.h>
#include <type_traits>

namespace gf2p {

[[nodiscard]] inline constexpr uint64_t get_primitive_poly(size_t bytes) noexcept {
    switch (bytes) {
        case 1:  return 0x1d; // GF(2^8)
        case 2:  return 0x2d; // GF(2^16)
        case 4:  return 0xc5; // GF(2^32)
        case 8:  return 0x1b; // GF(2^64)
        default: return 0x00;
    }
}

struct tables8 {
    std::array<uint16_t, 256> log{510};
    std::array<uint8_t, 1021> exp{};

    constexpr tables8() {
        uint8_t v = 1;
        for (unsigned i = 0; i < 255; ++i) {
            exp[i] = exp[i + 255] = v;
            log[v] = i;
            v = (v << 1) ^ (static_cast<uint8_t>(get_primitive_poly(1)) & -(v >> 7));
        }
    }
};
// inline constexpr tables8 T8{};

template <typename T>
struct gf2p {
    static_assert(std::is_unsigned_v<T>, "T must be an unsigned integer type");
    static_assert(sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 || sizeof(T) == 8, "Unsupported field size");

    T x;

    static constexpr T POLY = static_cast<T>(get_primitive_poly(sizeof(T)));
    static constexpr size_t BITS = sizeof(T) * 8;

    constexpr gf2p() noexcept = default;
    constexpr explicit gf2p(T val) noexcept : x(val) {}

    constexpr gf2p& operator+=(gf2p y) noexcept {
        x ^= y.x;
        return *this;
    }

    constexpr gf2p& operator-=(gf2p y) noexcept {
        x ^= y.x;
        return *this;
    }

    constexpr gf2p& operator*=(gf2p y) noexcept {
        *this = *this * y;
        return *this;
    }

    [[nodiscard]] friend constexpr gf2p operator+(gf2p a, gf2p b) noexcept {
        return a += b;
    }

    [[nodiscard]] friend constexpr gf2p operator-(gf2p a, gf2p b) noexcept {
        return a -= b;
    }

    [[nodiscard]] friend constexpr gf2p operator*(gf2p a, gf2p b) noexcept {
        T op_a = a.x;
        T op_b = b.x;
        T res = 0;

        for (size_t i = 0; i < BITS; ++i) {
            T b_mask = -(op_b & 1);
            res ^= (op_a & b_mask);
            op_b >>= 1;

            T overflow_mask = -(op_a >> (BITS - 1));
            op_a = (op_a << 1) ^ (POLY & overflow_mask);
        }

        return gf2p(res);
    }

    [[nodiscard]] constexpr bool operator==(gf2p other) const noexcept { return x == other.x; }
    [[nodiscard]] constexpr bool operator!=(gf2p other) const noexcept { return x != other.x; }

};

template <typename T>
[[nodiscard]] inline constexpr gf2p<T> mulx(gf2p<T> a) noexcept {
    return gf2p<T>((a.x << 1) ^ (a.POLY & (-(a.x >> (a.BITS - 1)))));
}

using gf2p8 = gf2p<uint8_t>;
using gf2p16 = gf2p<uint16_t>;
using gf2p32 = gf2p<uint32_t>;
using gf2p64 = gf2p<uint64_t>;

// static_assert(std::is_trivially_copyable_v<gf2p<uint32_t>>, "Must be trivially copyable");

} // namespace gf2p
