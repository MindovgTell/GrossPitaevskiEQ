#pragma once

#include <memory>

#include "grid.hpp"
#include "wavefunction.hpp"
#include "Core/traits.hpp"

namespace gpes {

    template <Dimension D, typename Algorythm,
                typename GridT = Grid<D>,
                typename WFT = WaveFunction<D>,
                typename = std::enable_if_t<
                    trait::is_grid_v<GridT> && trait::is_wavefunction_v<WFT> && 
                    trait::is_tstep_algo_v<Algorythm, WFT>
            >>    
    class TDSESolver {
    public:
        static constexpr Dimension Dim = D;
        using GridType = GridT;
        using WFType = WFT;
        using ShrdPtrGrid = std::shared_ptr<const Grid<D>>;

        // constructor
        TDSESolver() {

        }

        void run() {
            
        }

    private:
        WFType Psi_;
        ShrdPtrGrid grid_;
        Algorythm algo_;
        SimConfig cnfg;
    };
};