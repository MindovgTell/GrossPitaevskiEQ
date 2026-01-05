#pragma once

#include <cstddef>

namespace gpes {

    enum class Dimension : std::size_t {
        One = 1,
        Two = 2
    };

    template<Dimension D>
    inline constexpr std::size_t dim_v = static_cast<std::size_t>(D);  

    struct SimConfig {
        double duration;
        double dt;
        size_t num_of_steps;
    };

    struct PhysConfig {
        double a_s;
        double a_dd;
        double num_of_prt;
    };
} // namespace gpes
