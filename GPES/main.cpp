#include <chrono>
#include <regex>

#include "src/grid.hpp"
#include "src/utility.hpp"
#include "src/CrankNicolson.hpp"
#include <filesystem>

void simulation1D(std::string inputfile);
void simulation2D(std::string inputfile);
void simulation2Dt(std::string inputfile);
void simulation2DD(int argc, char const *argv[]);
void t_simulation2D(int argc, char const *argv[]);


void test();
void test1D();
void test1D_c();
void test2D();
void testfftw();
void tr();

void Ares_viz_run1(int argc, char const *argv[]);
void Ares_viz_run2(int argc, char const *argv[]);
void Ares_viz_run3(int argc, char const *argv[]);
void Ares_viz_run_1D(int argc, char const *argv[]);

void Ares_viz_energy1(int argc, char const *argv[]);
void Ares_viz_energy2(int argc, char const *argv[]);
void Ares_energy_eff1(int argc, char const *argv[]);
void Ares_energy_eff2(int argc, char const *argv[]);
void Ares_energy_eff3(int argc, char const *argv[]);
void Ares_energy_eff_1D(int argc, char const *argv[]);

int main(int argc, char const *argv[]){
    auto t0 = std::chrono::steady_clock::now();
    if(argc == 14){
        simulation2DD(argc, argv);
        // simulation2D(argv[1]);
    }
    else if(argc == 10){
        simulation1D(argv[1]);
    }
    else if(argc == 1){
        // Ares_viz_run_1D(argc, argv);
        // Ares_viz_energy(argc, argv);
        // Ares_energy_eff_1D(argc,argv);
        test1D();
        // testfftw();
        // tr();
    }
    else {
        std::string executable = argv[0];
        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }
    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "Simulation have been run for " << elapsed << " s\n";
    return 0;
}

std::vector<double> load_energy_csv(const std::filesystem::path& filepath) {
    std::vector<double> data;
    std::ifstream file(filepath);

    if (!file.is_open()) {
        std::cerr << "Failed to open " << filepath << "\n";
        return data;
    }

    std::string line;
    bool first_line = true;

    while (std::getline(file, line)) {
        if (first_line) {
            first_line = false; // skip header
            continue;
        }

        std::stringstream ss(line);
        std::string index_str, value_str;

        if (!std::getline(ss, index_str, ',')) continue;
        if (!std::getline(ss, value_str, ',')) continue;

        try {
            double value = std::stod(value_str);
            data.push_back(value);
        } catch (...) {
            data.push_back(0.0); // fallback on bad data
        }
    }

    return data;
}

void Ares_viz_energy(int argc, char const *argv[]) {
    std::filesystem::path root = "/Users/egormankevich/Desktop/GrossPitaevski/AresRes/grid/Grid_Size_Test2D_2";

    for (const auto& runDir : std::filesystem::directory_iterator(root)) {
        if (!runDir.is_directory()) continue;
        if (runDir.path().filename().string().rfind("run", 0) != 0) continue;

        for (const auto& a_sDir : std::filesystem::directory_iterator(runDir)) {
            if (!a_sDir.is_directory()) continue;

            for (const auto& stateDir : std::filesystem::directory_iterator(a_sDir)) {
                if (!stateDir.is_directory()) continue;
                if (stateDir.path().filename().string().rfind("finstate_", 0) != 0) continue;

                std::cout << ">>> Folder: " << stateDir.path() << "\n";

                std::filesystem::path fin_csv = stateDir.path() / "fin.csv";
                std::filesystem::path energy_csv = stateDir.path() / "energy.csv";

                // Load energy.csv if it exists
                if (std::filesystem::exists(energy_csv)) {
                    std::vector<double> energy_data = load_energy_csv(energy_csv);
                    std::filesystem::path energy_plot_path = stateDir.path() / "energy.png";
                    std::string s =  energy_plot_path.string();
                    GPES::draw_energy_save(energy_data, s);
                    if (!energy_data.empty()) {
                        double last = energy_data.back();
                        // use `last`
                    }
                }
            }
        }
    }
}

void Ares_energy_eff1(int argc, char const *argv[]) {
    std::filesystem::path root = "/Users/egormankevich/Desktop/GrossPitaevski/AresRes/data";

    std::map<std::string, std::map<int, double>> ordered_energy_by_case;  // key: a_s{...}e_dd{...} -> map<grid_size, energy>

    for (const auto& runDir : std::filesystem::directory_iterator(root)) {
        if (!runDir.is_directory()) continue;
        if (runDir.path().filename().string().rfind("run", 0) != 0) continue;

        for (const auto& a_sDir : std::filesystem::directory_iterator(runDir)) {
            if (!a_sDir.is_directory()) continue;

            for (const auto& stateDir : std::filesystem::directory_iterator(a_sDir)) {
                if (!stateDir.is_directory()) continue;
                const std::string folder = stateDir.path().filename().string();
                if (folder.rfind("finstate_", 0) != 0) continue;

                std::cout << ">>> Folder: " << stateDir.path() << "\n";

                std::filesystem::path energy_csv = stateDir.path() / "energy.csv";
                if (std::filesystem::exists(energy_csv)) {
                    std::vector<double> energy_data = load_energy_csv(energy_csv);
                    // std::filesystem::path energy_plot_path = stateDir.path() / "energy.png";
                    // GPES::draw_energy_save(energy_data, energy_plot_path.string());

                    if (!energy_data.empty()) {
                        double last_energy = energy_data.back();

                        // extract key: a_s{...}e_dd{...} and grid_size
                        std::string key;
                        int grid_size = -1;
                        std::regex pattern("a_s([0-9.]+)e_dd([0-9.]+)_grid_size([0-9]+)");
                        std::smatch match;
                        if (std::regex_search(folder, match, pattern) && match.size() == 4) {
                            key = "a_s" + match[1].str() + "e_dd" + match[2].str();
                            grid_size = std::stoi(match[3].str());
                        } else {
                            key = "unknown";
                        }

                        if (grid_size > 0) {
                            ordered_energy_by_case[key][grid_size] = last_energy;
                        }
                    }
                }
            }
        }
    }

    // // Optional: print collected results
    // for (const auto& [key, energy_map] : ordered_energy_by_case) {
    //     std::cout << "Group: " << key << " with " << energy_map.size() << " values\n";
    //     for (const auto& [grid, energy] : energy_map) {
    //         std::cout << "   Grid size " << grid << ": " << energy << "\n";
    //     }
    // }

 // Save grouped energy vectors as PNGs
    std::filesystem::path summary_dir = root / "summary_energy_plots";
    std::filesystem::create_directories(summary_dir);

    for (const auto& [key, energy_map] : ordered_energy_by_case) {
        std::vector<double> sorted_energies;
        for (const auto& [grid, energy] : energy_map) {
            sorted_energies.push_back(energy);
        }

        if (!sorted_energies.empty()) {
            double final_val = sorted_energies.back();
            if (final_val != 0.0) {
                for (auto& e : sorted_energies) {
                    e /= final_val;
                }
            }
        }

        std::filesystem::path out_plot = summary_dir / (key + ".png");
        std::string sss = out_plot.string();
        GPES::draw_energy_save_ars(sorted_energies, sss);
        std::filesystem::path out_csv = summary_dir / (key + ".csv");
        std::string ssss = out_csv.string();
        GPES::savecsv_vec(ssss, sorted_energies);
        std::cout << "Saved grouped energy plot: " << out_plot << "\n";
    }
}

void Ares_energy_eff2(int argc, char const *argv[]) {
    std::filesystem::path root = "/Users/egormankevich/Desktop/GrossPitaevski/AresRes/grid/Grid_Size_Test2D_2";

    std::map<std::string, std::map<int, double>> ordered_energy_by_case;  // key: a_s{...}e_dd{...} -> map<grid_size, energy>

    // Outer loop: a_s0.1, a_s0.01, ...
    for (const auto& a_sDir : std::filesystem::directory_iterator(root)) {
        if (!a_sDir.is_directory()) continue;
        if (a_sDir.path().filename().string().rfind("a_s", 0) != 0) continue;

        // Inner loop: finstate_... directories
        for (const auto& stateDir : std::filesystem::directory_iterator(a_sDir)) {
            if (!stateDir.is_directory()) continue;
            const std::string folder = stateDir.path().filename().string();
            if (folder.rfind("finstate_", 0) != 0) continue;

            std::cout << ">>> Folder: " << stateDir.path() << "\n";

            std::filesystem::path energy_csv = stateDir.path() / "energy.csv";
            if (std::filesystem::exists(energy_csv)) {
                std::vector<double> energy_data = load_energy_csv(energy_csv);

                if (!energy_data.empty()) {
                    double last_energy = energy_data.back();

                    // extract key: a_s{...}e_dd{...} and grid_size
                    std::string key;
                    int grid_size = -1;
                    std::regex pattern("a_s([0-9.]+)e_dd([0-9.]+)_grid_size([0-9]+)");
                    std::smatch match;
                    if (std::regex_search(folder, match, pattern) && match.size() == 4) {
                        key = "a_s" + match[1].str() + "e_dd" + match[2].str();
                        grid_size = std::stoi(match[3].str());
                    } else {
                        key = "unknown";
                    }

                    if (grid_size > 0) {
                        ordered_energy_by_case[key][grid_size] = last_energy;
                    }
                }
            }
        }
    }

 // Save grouped energy vectors as PNGs
    std::filesystem::path summary_dir = root / "summary_energy_plots";
    std::filesystem::create_directories(summary_dir);

    for (const auto& [key, energy_map] : ordered_energy_by_case) {
        std::vector<double> sorted_energies;
        for (const auto& [grid, energy] : energy_map) {
            sorted_energies.push_back(energy);
        }

        if (!sorted_energies.empty()) {
            double final_val = sorted_energies.back();
            if (final_val != 0.0) {
                for (auto& e : sorted_energies) {
                    e /= final_val;
                }
            }
        }

        std::filesystem::path out_plot = summary_dir / (key + ".png");
        std::string sss = out_plot.string();
        GPES::draw_energy_save_ars(sorted_energies, sss);
        std::filesystem::path out_csv = summary_dir / (key + ".csv");
        std::string ssss = out_csv.string();
        GPES::savecsv_vec(ssss, sorted_energies);
        std::cout << "Saved grouped energy plot: " << out_plot << "\n";
    }
}

void Ares_energy_eff3(int argc, char const *argv[]) {
    std::filesystem::path root = "/Users/egormankevich/Desktop/GrossPitaevski/AresRes/sim";

    // key: e_dd -> map<a_s, energy>
    std::map<std::string, std::map<double, double>> energy_by_e_dd; 

    // Outer loop: a_s0.1, a_s0.01, ...
    for (const auto& aDir : std::filesystem::directory_iterator(root)) {
        if (!aDir.is_directory()) continue;
        std::string dir_name = aDir.path().filename().string();
        if (dir_name.rfind("a_s", 0) != 0) continue;

        std::string a_s_str = dir_name.substr(4); // Extract value after "a_s"
        double a_s_value;
        try {
            a_s_value = std::stod(a_s_str);
        } catch (const std::exception& e) {
            std::cerr << "Invalid a_s value in folder name: " << a_s_str << "\n";
            continue;
        }

        // Inner loop: finstate_... directories
        for (const auto& stateDir : std::filesystem::directory_iterator(aDir)) {
            if (!stateDir.is_directory()) continue;
            const std::string folder = stateDir.path().filename().string();
            if (folder.rfind("finstate_", 0) != 0) continue;

            std::cout << ">>> Folder: " << stateDir.path() << "\n";

            std::filesystem::path energy_csv = stateDir.path() / "energy.csv";
            if (std::filesystem::exists(energy_csv)) {
                std::vector<double> energy_data = load_energy_csv(energy_csv);

                if (!energy_data.empty()) {
                    double last_energy = energy_data.back();

                    // extract e_dd value
                    std::regex pattern("e_dd([0-9eE+\-\.]+)");
                    std::smatch match;
                    if (std::regex_search(folder, match, pattern) && match.size() == 2) {
                        std::string e_dd_key = "e_dd" + match[1].str();
                        energy_by_e_dd[e_dd_key][a_s_value] = last_energy;
                    } else {
                        std::cerr << "⚠️ Could not extract e_dd from folder: " << folder << "\n";
                    }
                }
            }
        }
    }

    // Save energy vectors grouped by e_dd
    std::filesystem::path summary_dir = root / "summary_energy_by_edd";
    std::filesystem::create_directories(summary_dir);

    for (const auto& [e_dd_key, a_s_map] : energy_by_e_dd) {
        std::vector<double> a_s_vals;
        std::vector<double> energies;

        for (const auto& [a_s_val, energy] : a_s_map) {
            a_s_vals.push_back(a_s_val);
            energies.push_back(energy);
        }

        // Save to CSV
        std::filesystem::path out_csv = summary_dir / (e_dd_key + "_vs_a_s.csv");
        std::ofstream file(out_csv);
        if (!file.is_open()) {
            std::cerr << "Failed to write file: " << out_csv << "\n";
            continue;
        }
        file << "a_s,energy\n";
        for (size_t i = 0; i < a_s_vals.size(); ++i) {
            file << a_s_vals[i] << "," << energies[i] << "\n";
        }
        file.close();

        std::cout << "Saved energy values for " << e_dd_key << " to: " << out_csv << "\n";
    }
}
void Ares_viz_run1(int argc, char const *argv[]){
    std::filesystem::path root = "/Users/egormankevich/Desktop/GrossPitaevski/AresRes/data";

    // Outer loop: run1, run2, ...
    for (const auto& runDir : std::filesystem::directory_iterator(root)) {
        if (!runDir.is_directory()) continue;
        if (runDir.path().filename().string().rfind("run", 0) != 0) continue;

        std::cout << ">> Entering " << runDir.path() << "\n";

        // Middle loop: a_s0.1, a_s0.01, ...
        for (const auto& a_sDir : std::filesystem::directory_iterator(runDir)) {
            if (!a_sDir.is_directory()) continue;
            std::cout << " > Scanning " << a_sDir.path() << "\n";

            // Inner loop: finstate_... directories
            for (const auto& stateDir : std::filesystem::directory_iterator(a_sDir)) {
                if (!stateDir.is_directory()) continue;
                const std::string& stateName = stateDir.path().filename().string();

                if (stateName.rfind("finstate_", 0) != 0) continue;

                std::cout << "   >> Processing folder: " << stateDir.path() << "\n";

                // Now search for fin.csv inside finstate_... folder
                for (const auto& fileEntry : std::filesystem::directory_iterator(stateDir)) {
                    if (!fileEntry.is_regular_file()) continue;

                    const std::string fn = fileEntry.path().filename().string();

                    if (fn == "fin.csv") {
                        std::cout << "     - Found fin.csv\n";

                        GPES::WaveFunction<Dimension::Two> Psi;
                        Psi.readcsv(fileEntry.path().string());
                        // std::cout << fileEntry.path().string() << std::endl;
                        std::string outPath = (stateDir.path() / "fin.png").string();
                        GPES::heatmap_save(Psi, outPath);

                        std::cout << "     - Saved visualization to " << outPath << "\n";
                    }
                }
            }
        }
    }

}

void Ares_viz_run2(int argc, char const *argv[]) {
    std::filesystem::path root = "/Users/egormankevich/Desktop/GrossPitaevski/AresRes/grid/Grid_Size_Test2D_2";

    // Outer loop: a_s0.1, a_s0.01, ...
    for (const auto& a_sDir : std::filesystem::directory_iterator(root)) {
        if (!a_sDir.is_directory()) continue;
        if (a_sDir.path().filename().string().rfind("a_s", 0) != 0) continue;

        std::cout << ">> Entering " << a_sDir.path() << "\n";

        // Inner loop: finstate_... directories
        for (const auto& stateDir : std::filesystem::directory_iterator(a_sDir)) {
            if (!stateDir.is_directory()) continue;
            const std::string& stateName = stateDir.path().filename().string();

            if (stateName.rfind("finstate_", 0) != 0) continue;

            std::cout << "   >> Processing folder: " << stateDir.path() << "\n";

            // Now search for fin.csv inside finstate_... folder
            for (const auto& fileEntry : std::filesystem::directory_iterator(stateDir)) {
                if (!fileEntry.is_regular_file()) continue;

                const std::string fn = fileEntry.path().filename().string();

                if (fn == "fin.csv") {
                    std::cout << "     - Found fin.csv\n";

                    GPES::WaveFunction<Dimension::Two> Psi;
                    Psi.readcsv(fileEntry.path().string());
                    std::string outPath = (stateDir.path() / "fin.png").string();
                    GPES::heatmap_save(Psi, outPath);

                    std::cout << "     - Saved visualization to " << outPath << "\n";
                }
            }
        }
    }
}

void Ares_viz_run3(int argc, char const *argv[]) {
    std::filesystem::path root = "/Users/egormankevich/Desktop/GrossPitaevski/AresRes/Dy";

    // Outer loop: a_s0.1, a_s0.01, ...
    for (const auto& a_sDir : std::filesystem::directory_iterator(root)) {
        if (!a_sDir.is_directory()) continue;
        if (a_sDir.path().filename().string().rfind("a_s", 0) != 0) continue;

        std::cout << ">> Entering " << a_sDir.path() << "\n";

        // Inner loop: finstate_... directories
        for (const auto& stateDir : std::filesystem::directory_iterator(a_sDir)) {
            if (!stateDir.is_directory()) continue;
            const std::string& stateName = stateDir.path().filename().string();

            if (stateName.rfind("finstate_", 0) != 0) continue;

            std::cout << "   >> Processing folder: " << stateDir.path() << "\n";

            // Now search for fin.csv inside finstate_... folder
            for (const auto& fileEntry : std::filesystem::directory_iterator(stateDir)) {
                if (!fileEntry.is_regular_file()) continue;

                const std::string fn = fileEntry.path().filename().string();

                if (fn == "fin.csv") {
                    std::cout << "     - Found fin.csv\n";

                    GPES::WaveFunction<Dimension::Two> Psi;
                    Psi.readcsv(fileEntry.path().string());
                    std::string outPath = (stateDir.path() / "fin.png").string();
                    GPES::heatmap_save(Psi, outPath);
                    Eigen::VectorXd PsiX = Psi.get_x_slice();
                    Eigen::VectorXd PsiY = Psi.get_y_slice();
                    std::string outPathX = (stateDir.path() / "SliceX.png").string();
                    GPES::draw_save(PsiX, outPathX);
                    std::string outPathY = (stateDir.path() / "SliceY.png").string();
                    GPES::draw_save(PsiY, outPathY);

                    std::cout << "     - Saved visualization to " << outPath << "\n";
                }
            }
        }
    }
}


void Ares_viz_run_1D(int argc, char const *argv[]) {
    std::filesystem::path root = "/Users/egormankevich/Desktop/GrossPitaevski/AresRes/Dy/Dy_1D";

    // Outer loop: a_s0.1, a_s0.01, ...
    for (const auto& a_sDir : std::filesystem::directory_iterator(root)) {
        if (!a_sDir.is_directory()) continue;
        if (a_sDir.path().filename().string().rfind("a_s", 0) != 0) continue;

        std::cout << ">> Entering " << a_sDir.path() << "\n";

        // Inner loop: finstate_... directories
        for (const auto& stateDir : std::filesystem::directory_iterator(a_sDir)) {
            if (!stateDir.is_directory()) continue;
            const std::string& stateName = stateDir.path().filename().string();

            if (stateName.rfind("finstate_", 0) != 0) continue;

            std::cout << "   >> Processing folder: " << stateDir.path() << "\n";

            // Now search for fin.csv inside finstate_... folder
            for (const auto& fileEntry : std::filesystem::directory_iterator(stateDir)) {
                if (!fileEntry.is_regular_file()) continue;

                const std::string fn = fileEntry.path().filename().string();

                if (fn == "fin.csv") {
                    std::cout << "     - Found fin.csv\n";

                    GPES::WaveFunction<Dimension::One> Psi;
                    Psi.readcsv(fileEntry.path().string());
                    std::string outPath = (stateDir.path() / "fin.png").string();
                    GPES::draw_save(Psi, outPath);

                    std::cout << "     - Saved visualization to " << outPath << "\n";
                }
            }
        }
    }
}

void Ares_energy_eff_1D(int argc, char const *argv[]) {
    std::filesystem::path root = "/Users/egormankevich/Desktop/GrossPitaevski/AresRes/Grid_Size_Test1D_1";

    // key: a_s -> map<grid_size, energy>
    std::map<std::string, std::map<int, double>> energy_by_as;

    // Outer loop: a_s0.1, a_s0.01, ...
    for (const auto& aDir : std::filesystem::directory_iterator(root)) {
        if (!aDir.is_directory()) continue;
        std::string dir_name = aDir.path().filename().string();
        if (dir_name.rfind("a_s", 0) != 0) continue;

        std::string a_s_str = dir_name.substr(4); // Extract value after "a_s"

        // Inner loop: finstate_... directories
        for (const auto& stateDir : std::filesystem::directory_iterator(aDir)) {
            if (!stateDir.is_directory()) continue;
            const std::string folder = stateDir.path().filename().string();
            if (folder.rfind("finstate_", 0) != 0) continue;

            std::cout << ">>> Folder: " << stateDir.path() << "\n";

            std::filesystem::path energy_csv = stateDir.path() / "energy.csv";
            if (std::filesystem::exists(energy_csv)) {
                std::vector<double> energy_data = load_energy_csv(energy_csv);

                if (!energy_data.empty()) {
                    double last_energy = energy_data.back();

                    // extract grid size from folder name
                    std::regex pattern("grid_size([0-9]+)");
                    std::smatch match;
                    if (std::regex_search(folder, match, pattern) && match.size() == 2) {
                        int grid_size = std::stoi(match[1].str());
                        energy_by_as[a_s_str][grid_size] = last_energy;
                    } else {
                        std::cerr << "⚠️ Could not extract grid size from folder: " << folder << "\n";
                    }
                }
            }
        }
    }

    // Save energy vectors grouped by a_s
    std::filesystem::path summary_dir = root / "summary_energy_by_as";
    std::filesystem::create_directories(summary_dir);

    for (const auto& [a_s_key, grid_map] : energy_by_as) {
        std::vector<std::pair<int, double>> sorted_data(grid_map.begin(), grid_map.end());
        std::sort(sorted_data.begin(), sorted_data.end());

        std::filesystem::path out_csv = summary_dir / ("a_s" + a_s_key + ".csv");
        std::ofstream file(out_csv);
        if (!file.is_open()) {
            std::cerr << "Failed to write file: " << out_csv << "\n";
            continue;
        }
        file << "grid_size,energy\n";
        for (const auto& [grid_size, energy] : sorted_data) {
            file << grid_size << "," << energy << "\n";
        }
        file.close();

        std::cout << "Saved energy values for a_s=" << a_s_key << " to: " << out_csv << "\n";
    }
}

void Ares_ener_run1(){
    std::filesystem::path root = std::filesystem::path("/Users/egormankevich/Desktop/GrossPitaevski/AresRes/2D/Fin");

    // Outer loop: each finstate_* directory
    for (auto const& dirEntry : std::filesystem::directory_iterator(root)) {
        if (!dirEntry.is_directory()) 
            continue;

        auto const& dirName = dirEntry.path().filename().string();
        if (dirName.rfind("finstate_", 0) != 0) 
            continue;  // skip non-matching dirs

        std::cout << ">> Found folder: " << dirEntry.path() << "\n";

        std::cout << "=== Iterating folder: " << dirName << " ===\n";
        // Inner loop: files within that folder
        for (auto const& fileEntry : std::filesystem::directory_iterator(dirEntry.path())) {
            if (!fileEntry.is_regular_file())
                continue;

            std::string fn = fileEntry.path().filename().string();
            if (fn == "fin.csv") {
                std::cout << "    Found fin.csv — generating wavefunction image\n";
                GPES::WaveFunction<Dimension::Two> Psi;
                // 1. Read the data
                Psi.readcsv(fileEntry.path().string());
                // 3. Draw and save visualization as fin.png in this folder
                std::string outPath = (dirEntry.path() / "energy.png").string();
                // GPES::CrankNicolson<Dimension::Two> solver();
                // double TM_energy = solver.calc_state_energy(Psi);
                // std::vector<double> E = solver.get_vec_Energy();
                // std::vector<double> TM_en(E.size(), TM_energy);

                std::cout << "    Saved visualization to " << outPath << "\n";
                continue;
            }

            // match state_step_{number}.csv
            constexpr char const* prefix = "state_step_";
            constexpr char const* suffix = ".csv";
            if (fn.rfind(prefix, 0) == 0 && 
                fn.size() > strlen(prefix) + strlen(suffix) &&
                fn.substr(fn.size() - strlen(suffix)) == suffix)
            {
                std::string num = fn.substr(
                    strlen(prefix),
                    fn.size() - strlen(prefix) - strlen(suffix)
                );
                std::cout << "    processing " << fn 
                          << " (step = " << num << ")\n";

                // … your per-file code goes here …
            }
        }
    }

}

// void simulation2DD(int argc, char const *argv[]) {
//     GPES::Parameters p = GPES::ReadSimParameters2D(argc, argv);
//     if(p.empty()){
//         std::cout << "No parameters" << std::endl;
//         return;
//     }
//     GPES::Grid<Dimension::Two> grid(p.at("grid_size_x"),p.at("grid_size_y"), p.at("start_x"), p.at("start_y"));
//     grid.set_harmonic_potential(p.at("omega_x"), p.at("omega_y"));
//     grid.set_z_freq(1000);
//     double a_dd = p.at("e_dd") * p.at("a_s");
//     GPES::WaveFunction<Dimension::Two> Psi(grid,p.at("a_s"), a_dd, p.at("Num_of_Molecules"));
//     GPES::WaveFunction<Dimension::Two> Psi2(grid, p.at("a_s"), a_dd, p.at("Num_of_Molecules")); //(Grid, )
//     Psi.set_state_Gauss(0.0,0.0,p.at("sigma_x"),p.at("sigma_y"));
//     Psi2.set_state_TF(0.0, 0.0);
//     GPES::CrankNicolson<Dimension::Two> solver(Psi2,grid,p.at("time_step"),0.1);
//     std::cout << "g_scattering: " << solver.get_g_s() << std::endl;
//     std::cout << "g_lhy: " << solver.get_g_lhy() << std::endl;
//     std::cout << "C_dd: " << solver.get_C_dd() << std::endl;
//     std::string output_dir = argv[argc-1];
//     solver.simulation(output_dir);
//     GPES::WaveFunction<Dimension::Two> Fin;
//     solver.get_final_state(Fin);
//     // Saving final state
//     std::string output_file = output_dir + "/fin.csv";
//     Fin.savecsv_state(output_file);
// }

// void simulation1D(std::string inputfile){
//     std::ifstream input_data(inputfile);
//     std::string line;
//     double Problem,grid_size, deltat, T, start, number_of_mol, a_s, e_dd, x_c, sigma_x, omega_t;
//     std::vector<double> energies_for_diff_grid;
//     energies_for_diff_grid.reserve(8);
//     std::getline(input_data,line);
//     // Loop for running simulation and getting all the fin states
//     while(std::getline(input_data,line)) {//Skip first line in file
//         try { 
//             auto t0 = std::chrono::steady_clock::now();
//             std::stringstream str_stream(line);
//             if (str_stream >> Problem >> grid_size >> deltat >> T >> start >> number_of_mol >> e_dd >> a_s >> x_c >> sigma_x >> omega_t) {
//                     int width = 10;
//                 std::cout << std::setw(width) << "Problem" << std::setw(width) << "grid_size" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "start" << std::setw(width)
//                 << "N_of_mol"  << std::setw(width) << "e_dd" << std::setw(width) << "a_s" << std::setw(width) << "x_c" << std::setw(width) << "sigma_x" << std::setw(width) << "omega_t" << std::endl;
//                 std::cout << std::setw(width) << Problem << std::setw(width) << grid_size << std::setw(width) << deltat << std::setw(width) << T << std::setw(width) << start << std::setw(width)
//                 << number_of_mol << std::setw(width) << e_dd << std::setw(width) << a_s << std::setw(width) << x_c << std::setw(width) << sigma_x << std::setw(width) << omega_t << std::endl;
//             }
//             else {
//                 std::cout << "Error reading values from string!" << std::endl;
//                 return;
//             }
//             double a_dd = e_dd * a_s; 
//             //initialize grid
//             GPES::Grid<Dimension::One> grid(grid_size, start); // grid(Number_of_grid_nodes, starting_position)
//             grid.set_harmonic_potential(1);
//             grid.set_transverse(omega_t);
//             //initialize wavefunction
//             GPES::WaveFunction<Dimension::One> Psi(grid, a_s, a_dd, number_of_mol);
//             GPES::WaveFunction<Dimension::One> Psi2(grid, a_s, a_dd, number_of_mol); //(Grid, )
//             Psi.set_state_Gauss(x_c, sigma_x);
//             Psi2.set_state_TF(x_c);
//             //initialize simulator 
//             GPES::CrankNicolson<Dimension::One> solver(grid,Psi2, deltat, T);
//             //run simulation
//             solver.print_param_of_eq();
//             solver.simulation();
//             //get and save final state to csv
//             GPES::WaveFunction<Dimension::One> Fin;
//             solver.get_final_state(Fin);
//             Fin.savecsv_state("../../res/Fin_states/sim.csv");
//             //Calculating energy the energy during simulation
//             double TM_energy = solver.calc_state_energy(Psi2);
//             std::vector<double> E = solver.get_vec_Energy();
//             std::vector<double> TM_en(E.size(), TM_energy);
//             energies_for_diff_grid.push_back(E.back());
//             std::cout << E.back() << std::endl;
//             auto t1 = std::chrono::steady_clock::now();
//             double elapsed = std::chrono::duration<double>(t1 - t0).count();
//             std::cout << "Simulation have been run for " << elapsed << " s\n";
//             GPES::draw(Psi, Fin, Psi2);
//             GPES::draw_energy(E, TM_en);
//         }
//         catch (const std::exception& ex) {           // catches runtime_error and any std::exception
//             std::cout << "Error while saving wave-function: " << ex.what() << '\n';
//             return;                    // tell the OS something went wrong
//         }
//         catch (...)                                 // catches non-standard exceptions, just in case
//         {
//             std::cout << "Unknown error while saving wave-function.\n";
//             return;
//         }
//     }
//     // GPES::draw_energy(energies_for_diff_grid);
// }

// void simulation2D(std::string inputfile){
//     std::ifstream input_data(inputfile);
//     std::string line;
//     double Problem, grid_size_x, grid_size_y, deltat, T, start_x, start_y, number_of_mol, a_s, a_dd, x_c, y_c, sigma_x, sigma_y, omega_x, omega_y;
//     std::vector<double> energies_for_diff_grid;
//     energies_for_diff_grid.reserve(8);
//     std::getline(input_data,line);
//     while(std::getline(input_data,line)) {//Skip first line in file
//         std::stringstream str_stream(line);
//         if (str_stream >> Problem >> grid_size_x >> grid_size_y >> deltat >> T >> start_x >> start_y >> number_of_mol >> a_s >> a_dd >> x_c >> y_c >> sigma_x >> sigma_y >> omega_x >> omega_y) {
//             int width = 10;
//             std::cout << std::setw(width) << "Problem" << std::setw(width) << "grid_size_x" << std::setw(width) << "grid_size_y" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width)
//             << "start_x" << std::setw(width) << "start_y" << std::setw(width) << "N_of_mol"  << std::setw(width) << "a_s" << std::setw(width) << "a_dd" << std::setw(width)
//             << "x_c"  << std::setw(width) << "y_c" << std::setw(width) << "sigma_x" << std::setw(width) << "sigma_y" << std::setw(width) << "omega_x"<< std::setw(width) << "omega_y" << std::endl;
//             std::cout << std::setw(width) << Problem << std::setw(width) << grid_size_x << std::setw(width) << grid_size_y << std::setw(width) << deltat << std::setw(width) << T << std::setw(width)
//             << start_x << std::setw(width) << start_y << std::setw(width) << number_of_mol << std::setw(width) << a_s << std::setw(width) << a_dd << std::setw(width) << x_c << std::setw(width)
//             << y_c << std::setw(width) << sigma_x  << std::setw(width) << sigma_y << std::setw(width) << omega_x << std::setw(width) << omega_y << std::endl;
//         }
//         else
//             std::cout << "Error reading values from string!" << std::endl;
//         GPES::Grid<Dimension::Two> grid(grid_size_x,grid_size_y, start_x, start_y);
//         grid.set_harmonic_potential(omega_x, omega_y);
//         GPES::WaveFunction<Dimension::Two> Psi(grid,a_s, a_dd,number_of_mol);
//         GPES::WaveFunction<Dimension::Two> Psi2(grid, a_s, a_dd, number_of_mol); //(Grid, )
//         Psi.set_state_Gauss(x_c,y_c,sigma_x,sigma_y);
//         Psi2.set_state_TF(x_c, y_c);
//         GPES::CrankNicolson<Dimension::Two> solver(Psi2,grid,deltat,T);
//         // std::cout << "g_scattering: " << solver.get_g_s() << std::endl;
//         // std::cout << "g_lhy: " << solver.get_g_lhy() << std::endl;
//         // std::cout << "C_dd: " << solver.get_C_dd() << std::endl;
//         solver.simulation();
//         GPES::WaveFunction<Dimension::Two> Fin;
//         solver.get_final_state(Fin);
//         GPES::heatmap(Psi2, "Initial state");
//         GPES::heatmap(Fin, "Final state");
//         GPES::heatmap(Psi2, "Thomas-Fermi limit");
//         //Calculating energy the energy during simulation
//         double TM_energy = solver.calc_state_energy(Psi2);
//         std::vector<double> E = solver.get_vec_Energy();
//         std::vector<double> TM_en(E.size(), TM_energy);
//         energies_for_diff_grid.push_back(E.back());
//         GPES::draw_energy(E, TM_en);
//     }
//     GPES::draw_energy(energies_for_diff_grid);
// }
 
// void simulation2Dt(std::string inputfile){
//     std::ifstream input_data(inputfile);
//     std::string line;
//     double RunNum, grid_size_x, grid_size_y, deltat, T, start_x, start_y, number_of_mol, e_dd, a_s, x_c, y_c, sigma_x, sigma_y, omega_x, omega_y;
//     std::vector<double> energies_for_diff_grid;
//     energies_for_diff_grid.reserve(8);
//     std::getline(input_data,line);
//     while(std::getline(input_data,line)) {//Skip first line in file
//         std::stringstream str_stream(line);
//         if (str_stream >> RunNum >> grid_size_x >> grid_size_y >> deltat >> T >> start_x >> start_y >> number_of_mol >> e_dd >> a_s >> x_c >> y_c >> sigma_x >> sigma_y >> omega_x >> omega_y) {
//             int width = 10;
//             std::cout << std::setw(width) << "Problem" << std::setw(width) << "grid_size_x" << std::setw(width) << "grid_size_y" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width)
//             << "start_x" << std::setw(width) << "start_y" << std::setw(width) << "N_of_mol"  << std::setw(width) << "e_dd" << std::setw(width) << "a_s" << std::setw(width)
//             << "x_c"  << std::setw(width) << "y_c" << std::setw(width) << "sigma_x" << std::setw(width) << "sigma_y" << std::setw(width) << "omega_x"<< std::setw(width) << "omega_y" << std::endl;
//             std::cout << std::setw(width) << RunNum << std::setw(width) << grid_size_x << std::setw(width) << grid_size_y << std::setw(width) << deltat << std::setw(width) << T << std::setw(width)
//             << start_x << std::setw(width) << start_y << std::setw(width) << number_of_mol << std::setw(width) << e_dd << std::setw(width) << a_s << std::setw(width) << x_c << std::setw(width)
//             << y_c << std::setw(width) << sigma_x  << std::setw(width) << sigma_y << std::setw(width) << omega_x << std::setw(width) << omega_y << std::endl;
//         }
//         else
//             std::cout << "Error reading values from string!" << std::endl;
//         double a_dd = e_dd * a_s;
//         GPES::Grid<Dimension::Two> grid(grid_size_x,grid_size_y, start_x, start_y);
//         grid.set_harmonic_potential(omega_x, omega_y);
//         GPES::WaveFunction<Dimension::Two> Psi(grid,a_s, a_dd,number_of_mol);
//         GPES::WaveFunction<Dimension::Two> Psi2(grid, a_s, a_dd, number_of_mol); //(Grid, )
//         Psi.set_state_Gauss(x_c,y_c,sigma_x,sigma_y);
//         Psi2.set_state_TF(x_c, y_c);
//         GPES::CrankNicolson<Dimension::Two> solver(Psi,grid,deltat,T);
//         // // std::cout << "g_scattering: " << solver.get_g_s() << std::endl;
//         // // std::cout << "g_lhy: " << solver.get_g_lhy() << std::endl;
//         // // std::cout << "C_dd: " << solver.get_C_dd() << std::endl;
//         solver.simulation();
//         GPES::WaveFunction<Dimension::Two> Fin;
//         solver.get_final_state(Fin);
//         Eigen::VectorXd x_Fin_slice = Fin.get_x_slice();
//         Eigen::VectorXd y_Fin_slice = Fin.get_y_slice();
//         GPES::heatmap(Psi, "Initial state");
//         GPES::heatmap(Fin, "Final state");
//         GPES::draw(x_Fin_slice, "x", "Density");
//         GPES::draw(y_Fin_slice, "y", "Density");
//         GPES::heatmap(Psi2, "Thomas-Fermi limit");
//         //Calculating energy the energy during simulation
//         double TM_energy = solver.calc_state_energy(Psi2);
//         std::vector<double> E = solver.get_vec_Energy();
//         std::vector<double> TM_en(E.size(), TM_energy);
//         energies_for_diff_grid.push_back(E.back());
//         GPES::draw_energy(E, TM_en);
//     }
//     GPES::draw_energy(energies_for_diff_grid);
// }

void test(){
    try {
        //initialize the grid
        GPES::Grid<Dimension::Two> grid(300,300, -20, -20);
        grid.set_harmonic_potential(1,5);
        grid.set_z_freq(5);
        //initialize the wavefunction
        GPES::WaveFunction<Dimension::Two> Phi(grid, 0.00373, 0.00542, 5000);
        Phi.set_state_Gauss(0,0,1.5,1.5);
        Phi.print_params();
        GPES::CrankNicolson<Dimension::Two> solver(Phi,grid,0.001,0.1);
        solver.print_param_of_eq();
        std::string outdir = "/Users/egormankevich/Desktop/GrossPitaevski/res/2D/test/test_a_s0.004_e_dd0.0047_N1000_11";
        solver.simulation(outdir);
        GPES::WaveFunction<Dimension::Two> Fin;
        solver.get_final_state(Fin);
        GPES::heatmap(Fin);
        std::string output_fin = outdir + "/fin.png";
        GPES::heatmap_save(Fin,output_fin);
        double TM_energy = solver.calc_state_energy(Phi);
        std::vector<double> E = solver.get_vec_Energy();
        std::vector<double> TM_en(E.size(), TM_energy);
        GPES::draw_energy(E,TM_en);
        Eigen::VectorXd FinX = Fin.get_x_slice();
        Eigen::VectorXd FinY = Fin.get_y_slice();
        Eigen::VectorXd Psi2X = Phi.get_x_slice();
        Eigen::VectorXd Psi2Y = Phi.get_y_slice();
        std::string output_sliceX_png = outdir + "/sliceX.png";
        GPES::draw_save(FinX, output_sliceX_png);
        std::string output_sliceY_png = outdir + "/sliceY.png";
        GPES::draw_save(FinY, output_sliceY_png);
        std::string output_energy = outdir + "/energy.png";
        GPES::draw_energy_save(E, TM_en, output_energy);

    }
    catch (const std::exception& ex) {           // catches runtime_error and any std::exception
        std::cout << "Error while saving wave-function: " << ex.what() << '\n';
        return;                    // tell the OS something went wrong
    }
    catch (...)                                 // catches non-standard exceptions, just in case
    {
        std::cout << "Unknown error while saving wave-function.\n";
        return;
    }
}


void test2D(){
    try {

        GPES::Grid<Dimension::Two> grid(800,800, -5, -5);
        grid.set_harmonic_potential(1,1);
        grid.set_z_freq(10);
        //initialize the wavefunction
        GPES::WaveFunction<Dimension::Two> Phi(grid, 0.00337, 0.00534, 15000);
        Phi.set_state_Gauss(0,0,1,1);
        Phi.print_params();
        GPES::CrankNicolson<Dimension::Two> solver(Phi,grid,0.001,0.1);
        Eigen::VectorXd v = solver.get_U_ddi();
        GPES::heatmap_pot(v,"Potential");
        GPES::savecsv_vec("/Users/egormankevich/Desktop/GrossPitaevski/res/2D/potential23.csv", v);
    }
    catch (const std::exception& ex) {           // catches runtime_error and any std::exception
        std::cout << "Error while saving wave-function: " << ex.what() << '\n';
        return;                    // tell the OS something went wrong
    }
    catch (...)                                 // catches non-standard exceptions, just in case
    {
        std::cout << "Unknown error while saving wave-function.\n";
        return;
    }
}



void test1D(){
    try {
        //initialize the grid
        GPES::Grid<Dimension::One> grid(1200, -10);
        grid.set_harmonic_potential(1);
        grid.set_transverse(2);

            // initialize the wavefunction
        GPES::WaveFunction<Dimension::One> Phi(grid, 0.0042, 0.00534, 8000);
        Phi.set_state_Gauss(0,0.5);
        Phi.print_params();

        // GPES::WaveFunction<Dimension::One> Phi;
        // Phi.readcsv("/Users/egormankevich/Desktop/GrossPitaevski/res/1D/test_a_s0.004/fin.csv");
        // Phi.set_Num(4000);

        GPES::CrankNicolson<Dimension::One> solver(grid, Phi, 0.001, 0.1);
        solver.print_param_of_eq();
        std::string outdir = "/Users/egormankevich/Desktop/GrossPitaevski/res/1D/test_a_s0.004_12";
        solver.simulation(outdir);
        GPES::WaveFunction<Dimension::One> Fin;
        solver.get_final_state(Fin);

        GPES::draw(Fin);
        std::string output_fin = outdir + "/fin.png";
        GPES::draw_save(Fin,output_fin);
        
        double TM_energy = solver.calc_state_energy(Phi);
        std::vector<double> E = solver.get_vec_Energy();
        std::vector<double> TM_en(E.size(), TM_energy);
        GPES::draw_energy(E,TM_en);
        std::string output_energy = outdir + "/energy.png";
        GPES::draw_energy_save(E, TM_en, output_energy);
    }
    catch (const std::exception& ex) {           // catches runtime_error and any std::exception
        std::cout << "Error while saving wave-function: " << ex.what() << '\n';
        return;                    // tell the OS something went wrong
    }
    catch (...)                                 // catches non-standard exceptions, just in case
    {
        std::cout << "Unknown error while saving wave-function.\n";
        return;
    }
}

void test1D_c(){
    try {
        //initialize the grid
        std::string inputdir = "/Users/egormankevich/Desktop/GrossPitaevski/res/1D/test_a_s0.004_7/fin.csv";
        GPES::Grid<Dimension::One> grid(300, -20);
        grid.set_harmonic_potential(1);
        grid.set_transverse(5);

            // initialize the wavefunction
        GPES::WaveFunction<Dimension::One> Phi(grid, 0.00371, 0.00534, 30000);
        Phi.set_state_Gauss(0,5);
        Phi.readcsv(inputdir);

        Phi.print_params();


        // GPES::WaveFunction<Dimension::One> Phi;
        // Phi.readcsv("/Users/egormankevich/Desktop/GrossPitaevski/res/1D/test_a_s0.004/fin.csv");
        // Phi.set_Num(4000);

        GPES::CrankNicolson<Dimension::One> solver(grid, Phi, 0.001, 0.1);
        solver.print_param_of_eq();
        std::string outdir = "/Users/egormankevich/Desktop/GrossPitaevski/res/1D/test_a_s0.004_7_1";
        solver.simulation(outdir);
        GPES::WaveFunction<Dimension::One> Fin;
        solver.get_final_state(Fin);

        GPES::draw(Fin);
        std::string output_fin = outdir + "/fin.png";
        GPES::draw_save(Fin,output_fin);
        
        double TM_energy = solver.calc_state_energy(Phi);
        std::vector<double> E = solver.get_vec_Energy();
        std::vector<double> TM_en(E.size(), TM_energy);
        GPES::draw_energy(E,TM_en);
        std::string output_energy = outdir + "/energy.png";
        GPES::draw_energy_save(E, TM_en, output_energy);
    }
    catch (const std::exception& ex) {           // catches runtime_error and any std::exception
        std::cout << "Error while saving wave-function: " << ex.what() << '\n';
        return;                    // tell the OS something went wrong
    }
    catch (...)                                 // catches non-standard exceptions, just in case
    {
        std::cout << "Unknown error while saving wave-function.\n";
        return;
    }
}

void testfftw(){
    GPES::WaveFunction<Dimension::One> Psi;
    // Psi.readcsv("/Users/egormankevich/Desktop/GrossPitaevski/AresRes/sim/a_s0.002/finstate_e_dd1.2_NumOfMol1000/fin.csv");
    Psi.readcsv("/Users/egormankevich/Desktop/GrossPitaevski/res/1D/test_a_s0.004_11/fin.csv");
    // Psi.set_state_Gauss(0,0,5,5);
    GPES::draw(Psi);
    Eigen::VectorXcd K = Psi.momentum_space_transform();
    Psi.set_vec(K);
    // GPES::draw(Psi);
    GPES::savecsv_vec("/Users/egormankevich/Desktop/GrossPitaevski/res/1D/test_a_s0.004_11/FFT.csv", Psi.prob());
}




std::vector<double> V_dd(GPES::Grid<Dimension::One>& grid, double l_perp){
    double step = grid.get_step_size();
    double start = grid.get_start_position();
    double N = grid.get_size_of_grid();

    std::vector<double> V(N,0.0);

    // double l_perp = 1.0; // /std::sqrt(1000);

    for(int i = 0; i < N; ++i){
        double x = start + step * i;

        double xi = std::abs(x) / l_perp;
        double xi2 = xi*xi;

        double alfa = std::abs(xi) / (std::sqrt(2.0));
        double sqalfa = 0.5 * std::pow(xi, 2);


        double erfc_val = std::erfc(alfa);
        double exp_val = std::exp(sqalfa);

        if(l_perp == 10.0)
            V[i] = 4.0 / std::pow(std::abs(x),3);
        else
            V[i] = -1.* (2.0 * std::abs(xi) - std::sqrt(2.0 * M_PI) * (1.0 + xi2) * exp_val * erfc_val); //-0.75 * 0.00537 /(std::pow(l_perp,3)) *
    }

    return V;
}

double V_dd_num(GPES::Grid<Dimension::One>& grid, double x){
    double step = grid.get_step_size();
    double start = grid.get_start_position();
    double N = grid.get_size_of_grid();

    double V;

    double l_perp = 1.0/std::sqrt(1000);

    double xi = std::abs(x) / l_perp;
    double xi2 = xi*xi;

    double alfa = std::abs(xi) / (std::sqrt(2.0));
    double sqalfa = 0.5 * std::pow(xi, 2);


    double erfc_val = std::erfc(alfa);
    double exp_val = std::exp(sqalfa);

    V = -0.75 * 0.00537 /(std::pow(l_perp,3)) *(2.0 * std::abs(xi) - std::sqrt(2.0 * M_PI) * (1.0 + xi2) * exp_val * erfc_val);

    return V;
}

void tr(){

    GPES::Grid<Dimension::One> grid(300, -5);
    std::vector<double> V1 = V_dd(grid, 1.0);
    std::vector<double> V2 = V_dd(grid, 0.4);
    std::vector<double> V3 = V_dd(grid, 10.0);
    // double Vn = V_dd_num(grid,10);
    GPES::drawe(V1, V2, V3, "$x / l_x$", "$V_{\\mathrm{dd}}(|x|/{l_{\\perp}})$");
}
// void t_simulation2D(int argc, char const *argv[]) {
//     GPES::Parameters p = GPES::ReadSimParameters(argc, argv);
//     if(p.empty()){
//         std::cout << "No parameters" << std::endl;
//         return;
//     }
//     GPES::Grid<Dimension::Two> grid(p.at("grid_size_x"),p.at("grid_size_y"), p.at("start_x"), p.at("start_y"));
//     grid.set_harmonic_potential(p.at("omega_x"), p.at("omega_y"));
//     grid.set_z_freq(1000);
//     double a_dd = p.at("e_dd") * p.at("a_s");
//     GPES::WaveFunction<Dimension::Two> Psi(grid,p.at("a_s"), a_dd, p.at("Num_of_Molecules"));
//     GPES::WaveFunction<Dimension::Two> Psi2(grid, p.at("a_s"), a_dd, p.at("Num_of_Molecules")); //(Grid, )
//     Psi.set_state_Gauss(0.0,0.0,p.at("sigma_x"),p.at("sigma_y"));
//     Psi2.set_state_TF(0.0, 0.0);
//     GPES::CrankNicolson<Dimension::Two> solver(Psi2,grid,p.at("time_step"),0.1);
//     std::cout << "g_scattering: " << solver.get_g_s() << std::endl;
//     std::cout << "g_lhy: " << solver.get_g_lhy() << std::endl;
//     std::cout << "C_dd: " << solver.get_C_dd() << std::endl;
//     solver.simulation();
//     GPES::WaveFunction<Dimension::Two> Fin;
//     solver.get_final_state(Fin);
//     // Saving final state
//     std::string output_dir = argv[argc-1];
//     std::string output_file = output_dir + "/fin.csv";
//     Psi.savecsv_state(output_file);
//     //Calculating energy the energy during simulation
//     double TM_energy = solver.calc_state_energy(Psi2);
//     std::vector<double> E = solver.get_vec_Energy();
//     std::vector<double> TM_en(E.size(), TM_energy);
//     //Draw densities of the func and energy by evolution step
//     std::string output_png = output_dir + "/fin.png";
//     GPES::heatmap_save(Fin, output_png);
//     Eigen::VectorXd FinX = Fin.get_x_slice();
//     Eigen::VectorXd FinY = Fin.get_y_slice();
//     Eigen::VectorXd Psi2X = Psi2.get_x_slice();
//     Eigen::VectorXd Psi2Y = Psi2.get_y_slice();
//     std::string output_sliceX_png = output_dir + "/sliceX.png";
//     GPES::draw_save(FinX, Psi2X, output_sliceX_png);
//     std::string output_sliceY_png = output_dir + "/sliceY.png";
//     GPES::draw_save(FinY, Psi2Y,output_sliceY_png);
//     std::string output_energy = output_dir + "/energy.png";
//     GPES::draw_energy_save(E, TM_en, output_energy);
// }