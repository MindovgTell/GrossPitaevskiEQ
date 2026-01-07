#pragma once

#include <memory>
#include <stdexcept>

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

        TDSESolver(ShrdPtrGrid grid, Algorythm& algo, WFType& psi, SimConfig simcnfg)
            : Psi_(psi),
              grid_(std::move(grid)),
              algo_(algo),
              sim_config_(simcnfg) {
            if (!grid_) {
                throw std::runtime_error("TDSESolver requires a valid grid");
            }
        }

        void run() {
            for (size_t step = 0; step < sim_config_.num_of_steps; ++step) {
                algo_.step(Psi_);
            }
        }

        const WFType& wavefunction() const { return Psi_; }
        WFType& wavefunction() { return Psi_; }

    private:
        WFType Psi_
        ShrdPtrGrid grid_;
        Algorythm algo_;
        SimConfig sim_config_;
    };
};
