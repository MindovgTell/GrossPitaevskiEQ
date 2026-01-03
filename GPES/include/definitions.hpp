#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

namespace GPES {

    enum class Dimension {
        One = 1,
        Two = 2
    };


    struct SimulationParams {
        double duration;
        double time_step;
        double a_s;
        double a_dd;
    };

    struct QuenchParams { 
    };
}

#endif