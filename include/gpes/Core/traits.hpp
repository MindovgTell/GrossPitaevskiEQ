/*
* In future could be separated into different files for SFINAE type traits and tag info
*/

#pragma once

#include <type_traits>
#include <utility>
#include <string_view>

namespace gpes::trait {

    template <typename...>
    using void_t = void;

    // Check Grid-like requirements with size() method
    template <typename G, typename = void>
    struct is_grid : std::false_type {};

    template <typename G>
    struct is_grid<G, void_t<
        decltype(G::Dim),
        typename G::PotVec,
        decltype(std::declval<G&>().size()) 
    >> : std::true_type {};

    template <typename G>
    inline constexpr bool is_grid_v = is_grid<G>::value;

    // Check WaveFunction-like requirements with vec(), data(), size() methods 
    template <typename WF, typename = void>
    struct is_wavefunction : std::false_type {};

    template <typename WF>
    struct is_wavefunction<WF, void_t<
        decltype(WF::Dim),
        typename WF::Scalar,
        typename WF::VectorType,
        decltype(std::declval<WF&>().vec()),
        decltype(std::declval<const WF&>().vec()),
        decltype(std::declval<WF&>().data()),
        decltype(std::declval<const WF&>().data()),
        decltype(std::declval<const WF&>().size()) 
    >> : std::true_type {};

    template <typename WF>
    inline constexpr bool is_wavefunction_v = is_wavefunction<WF>::value;

    // Stronger check: data() return type exactly matches Scalar* / const Scalar*
    template <typename WF, typename = void>
    struct has_correct_data_ptr : std::false_type {};

    template <typename WF>
    struct has_correct_data_ptr<WF, void_t<typename WF::Scalar>> {
    private:
        using S = typename WF::Scalar;
    public:
        inline constexpr bool value =
            std::is_same_v<decltype(std::declval<WF&>.data()), S*> &&
            std::is_same_v<decltype(std::declval<const WF&>.data()), const S*> &&
    };

    template <typename WF>
    inline constexpr bool has_correct_data_ptr_v = has_correct_data_ptr<WF>::value;

    // Algorithm must have static constexpr Dimension dim and step(WF&, double)
    template <typename Algo, typename WF, typename = void>
    struct is_tstep_algo : std::false_type {};

    template <typename Algo, typename WF>
    struct is_tstep_algo<Algo, WF, void_t<
        decltype(Algo::Dim),
        decltype(std::declval<Algo&>().step(std::declval<WF&>()))
    >> : std::true_type {};
    
    template <typename Algo, typename WF>
    inline constexpr bool is_tstep_algo_v = is_tstep_algo<Algo, WF>::value;

};


// Tag info for solvers 
namespace gpes::tags {

// algorithms
struct CrankNicolson {};
struct SSFM {};

// // in future could be added boundary conditions, actually support only Dirichlet
// struct Dirichlet {};
// struct Periodic {};

} // namespace gpes::tags

namespace gpes {

template <typename AlgoTag>
struct algo_traits;

template <>
struct algo_traits<tags::CrankNicolson> {
    static constexpr bool implicit = true;
    static constexpr bool requires_fft = false;
    static constexpr std::string_view name() { return "Crank-Nicolson"; }
};

template <>
struct algo_traits<tags::SSFM> {
    static constexpr bool implicit = false;
    static constexpr bool requires_fft = true;
    static constexpr std::string_view name() { return "Split-Step Fourier (SSFM)"; }
};

} // namespace gpes
