#ifndef UTIL_H__
#define UTIL_H__

#include <utility>
#include <functional>
#include <array>
#include <limits>
#include <type_traits>
#include <stdexcept>

#define STR_INNER_(x) #x
#define STR(x) STR_INNER_(x)

// Make a helper function so we can set breakpoints
namespace detail {
inline void thrower_(const char* msg)
{
    throw std::runtime_error(msg);
}
}

#define ensureOrThrow(expr) do { \
    if (!(expr)) { \
        detail::thrower_(__FILE__ ":" STR(__LINE__) ":'" #expr "' was not true"); \
    } \
} while(false)

namespace std {
    // We need to make our hash specialization since std::pair does not have one
    template <>
    struct hash<std::pair<int, int>> {
        size_t operator()(const std::pair<int, int>& key) const noexcept;
    };

    template <size_t N>
    struct hash<std::array<int, N>> {
        size_t operator()(const std::array<int, N>& key) const noexcept
        {
            // Code from https://codereview.stackexchange.com/q/171999
            static std::hash<int> hasher;
            size_t result = 0;
            for (int k : key) {
                result = result * 31 + hasher(k);
            }
            return result;
        }
    };
} // std


void graphErrFunc(const char *msg);

template <class Ty>
constexpr std::enable_if_t<std::is_floating_point<Ty>::value, Ty> infOrMax() {
    return std::numeric_limits<Ty>::infinity();
}

template <class Ty>
constexpr std::enable_if_t<!std::is_floating_point<Ty>::value, Ty> infOrMax() {
    return std::numeric_limits<Ty>::max();
}

template <class Ty>
constexpr std::enable_if_t<std::is_floating_point<Ty>::value, Ty> minusInfOrMin() {
    return -std::numeric_limits<Ty>::infinity();
}

template <class Ty>
constexpr std::enable_if_t<!std::is_floating_point<Ty>::value, Ty> minusInfOrMin() {
    return std::numeric_limits<Ty>::min();
}

#endif // UTIL_H__
