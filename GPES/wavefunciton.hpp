#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP


#include "definitions.hpp"


namespace GPES
{
    template<Dimension dim>
    class WaveFunction;

    template<>
    class WaveFunction<Dimension::One>{

    };


    template<>
    class WaveFunction<Dimension::Two>{

    };
    
} // namespace GPES


#endif
