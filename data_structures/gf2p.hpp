#pragma once
#include <cstddef>
#include <cstdint>
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

    [[nodiscard]] inline constexpr gf2p mulx(gf2p a) const noexcept {
        T op_a = a.x;
        T mask = -(op_a >> (BITS - 1));

        return gf2p((op_a << 1) ^ (POLY & mask));
    }
};

using gf2p8 = gf2p<uint8_t>;
using gf2p16 = gf2p<uint16_t>;
using gf2p32 = gf2p<uint32_t>;
using gf2p64 = gf2p<uint64_t>;

// static_assert(std::is_trivially_copyable_v<gf2p<uint32_t>>, "Must be trivially copyable");

} // namespace gf2p
