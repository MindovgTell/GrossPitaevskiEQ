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

    // Check Grid-like requirements with 1D or 2D grid APIs.
    template <typename G, typename = void>
    struct has_1d_grid_api : std::false_type {};

    template <typename G>
    struct has_1d_grid_api<G, void_t<
        decltype(std::declval<const G&>().size()),
        decltype(std::declval<const G&>().step()),
        decltype(std::declval<const G&>().start())
    >> : std::true_type {};

    template <typename G, typename = void>
    struct has_2d_grid_api : std::false_type {};

    template <typename G>
    struct has_2d_grid_api<G, void_t<
        decltype(std::declval<const G&>().size_x()),
        decltype(std::declval<const G&>().size_y()),
        decltype(std::declval<const G&>().step_x()),
        decltype(std::declval<const G&>().step_y()),
        decltype(std::declval<const G&>().start_pos_x()),
        decltype(std::declval<const G&>().start_pos_y())
    >> : std::true_type {};

    template <typename G, typename = void>
    struct is_grid : std::false_type {};

    template <typename G>
    struct is_grid<G, void_t<
        decltype(G::Dim),
        decltype(std::declval<const G&>().potential())
    >> : std::bool_constant<has_1d_grid_api<G>::value || has_2d_grid_api<G>::value> {};

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
        inline static constexpr bool value =
            std::is_same_v<decltype(std::declval<WF&>().data()), S*> &&
            std::is_same_v<decltype(std::declval<const WF&>().data()), const S*>;
    };

    template <typename WF>
    inline constexpr bool has_correct_data_ptr_v = has_correct_data_ptr<WF>::value;

    // Algorithm must have static constexpr Dimension dim, step(WF&), and last_energy()
    template <typename Algo, typename WF, typename = void>
    struct is_tstep_algo : std::false_type {};

    template <typename Algo, typename WF>
    struct is_tstep_algo<Algo, WF, void_t<
        decltype(Algo::Dim),
        decltype(std::declval<Algo&>().step(std::declval<WF&>())),
        decltype(std::declval<const Algo&>().last_energy())
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
