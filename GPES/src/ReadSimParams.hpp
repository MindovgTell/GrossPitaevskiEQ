#ifndef READSIMPARAMS_HPP
#define READSIMPARAMS_HPP

#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <string>
// #include <sstream>
#include <vector>
// #include <fstream>
#include "definitions.hpp"

namespace GPES
{

// using ParameterList         =       std::map<std::string, std::vector<double>>;
using Parameters            =       std::unordered_map<std::string, double>;


template<Dimension Dim>
class ReadSimulationParameters;

template<>
class ReadSimulationParameters<Dimension::One> {
public: 
    struct SimulationParameters{
        double number_of_molecules;
        double grid_size_X;
        double delta_t;
        double start_position_X;
        double e_dd, a_s;
        double x_c;
        double sigma_x;
        double omega_x, omega_transverse;
    };
private:
    SimulationParameters _params;
    void parseFile(const std::string& filename);
    void validateHeader(const std::string& headerLine);

public:
    // explicit ReadSimulationParameters(int argc, char const *argv[]) {
    //     std::vector<double> params;
    //     params.reserve(argc - 1);

    //     for (int i = 1; i < argc; ++i) {
    //         try {
    //             double d = std::stod(argv[i]);
    //             params.push_back(d);
    //         }
    //         catch (const std::invalid_argument&) {
    //             std::cerr << "Invalid number: \"" << argv[i] << "\"\n";
    //             return 1;
    //         }
    //         catch (const std::out_of_range&) {
    //             std::cerr << "Number out of range: \"" << argv[i] << "\"\n";
    //             return 1;
    //         }
    //     }

    // // now params[0]…params.back() are your doubles
    // for (double d : params) 
    //     std::cout << d << "\n";


    //     // std::ifstream input_data(filename);
    //     // double RunNum, grid_size_x, grid_size_y, deltat, T, start_x, start_y, number_of_mol, e_dd, a_s, x_c, y_c, sigma_x, sigma_y, omega_x, omega_y;

    //     // std::getline(input_data,line);
    //     // while(std::getline(input_data,line)) {//Skip first line in file

    //     //     std::stringstream str_stream(line);

    //     //     str_stream >> RunNum >> grid_size_x >> grid_size_y >> deltat >> T >> start_x >> start_y >> number_of_mol >> e_dd >> a_s >> x_c >> y_c >> sigma_x >> sigma_y >> omega_x >> omega_y;
    //     // }
    // }

    // // Returns the parsed parameter sets
    // const std::vector<SimulationParameters>& getParameters() const;
};


template<>
class ReadSimulationParameters<Dimension::Two> {
private:
    struct SimulationParameters{
        double number_of_molecules;
        double runNum;
        double gridSizeX, gridSizeY;
        double deltaT, totalTime;
        double startX, startY;
        double numberOfMolecules;
        double e_dd, a_s;
        double x_c, y_c;
        double sigma_x, sigma_y;
        double omega_x, omega_y;
    };

public: 

};

Parameters ReadSimParameters(int& argc, char const *argv[]){
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <num1> [num2] …\n";
        return {};
    }   
    
    std::vector<double> params;
    params.reserve(argc - 2);

    for (int i = 1; i < argc-1; ++i) {
        try {
            double d = std::stod(argv[i]);
            params.push_back(d);
        }
        catch (const std::invalid_argument&) {
            std::cout << "Invalid number: \"" << argv[i] << "\"\n";
            return {};
        }
        catch (const std::out_of_range&) {
            std::cout << "Number out of range: \"" << argv[i] << "\"\n";
            return {};
        }
    }

    Parameters parameters{
        {"grid_size_x",        params[0] },
        {"grid_size_y",        params[1] },
        {"time_step",          params[2] },
        {"start_x",            params[3] },
        {"start_y",            params[4] },
        {"Num_of_Molecules",   params[5] },
        {"e_dd",               params[6] },
        {"a_s",                params[7] },
        {"sigma_x",            params[8] },
        {"sigma_y",            params[9] },
        {"omega_x",            params[10]},
        {"omega_y",            params[11]}
    };

    // now params[0]…params.back() are your doubles
    int width = 15;
    std::cout << "Parameters where readed successfully" << std::endl;
    for (auto d : parameters) 
        std::cout << std::setw(width) << d.first;
    std::cout << '\n';
    for (auto d : parameters)
        std::cout << std::setw(width) << d.second;
    std::cout << '\n';

    return parameters;
}

} // namespace GPES



#endif