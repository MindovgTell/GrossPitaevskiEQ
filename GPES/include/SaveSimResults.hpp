#pragma once 


namespace GPES
{
    
class SaveSimResults {
private: 


public:

};



void savecsv_vec(std::string file_path, std::vector<double>& v){
    // 1) Ensure parent directory exists
    std::string dir = parent_dir(file_path);
    if (!dir.empty()) {
        // mkdir -p dir
        std::string cmd = "mkdir -p '" + dir + "'";
        if (std::system(cmd.c_str()) != 0) {
            throw std::runtime_error("Failed to create directory: " + dir);
        }
    }

    std::ofstream file(file_path, std::ios::out /*| std::ios::trunc is implicit*/);
    if (!file) {
        std::cerr << "Error: cannot open " << file_path << " for writing\n";
        return;
    }

    /* ---------- 3) the actual wave-function ---------- */
    file << "index,value\n";
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        file << i << ',' << v[i] << '\n';
    }

    std::cout << "Vector have been saved "<< file_path << std::endl;
}

void savecsv_vec(std::string file_path, Eigen::VectorXd v){
    // 1) Ensure parent directory exists
    std::string dir = parent_dir(file_path);
    if (!dir.empty()) {
        // mkdir -p dir
        std::string cmd = "mkdir -p '" + dir + "'";
        if (std::system(cmd.c_str()) != 0) {
            throw std::runtime_error("Failed to create directory: " + dir);
        }
    }

    std::ofstream file(file_path, std::ios::out /*| std::ios::trunc is implicit*/);
    if (!file) {
        std::cerr << "Error: cannot open " << file_path << " for writing\n";
        return;
    }

    /* ---------- 3) the actual wave-function ---------- */
    file << "index,value\n";
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        file << i << ',' << v[i] << '\n';
    }

    std::cout << "Vector have been saved "<< file_path << std::endl;
}

void savecsv_vec(std::string file_path, Eigen::VectorXcd v){
    // 1) Ensure parent directory exists
    std::string dir = parent_dir(file_path);
    if (!dir.empty()) {
        // mkdir -p dir
        std::string cmd = "mkdir -p '" + dir + "'";
        if (std::system(cmd.c_str()) != 0) {
            throw std::runtime_error("Failed to create directory: " + dir);
        }
    }

    std::ofstream file(file_path, std::ios::out /*| std::ios::trunc is implicit*/);
    if (!file) {
        std::cerr << "Error: cannot open " << file_path << " for writing\n";
        return;
    }

    /* ---------- 3) the actual wave-function ---------- */
    file << "index,value\n";
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        file << i << ',' << v[i] << '\n';
    }

    std::cout << "Vector have been saved "<< file_path << std::endl;
}

} // namespace GPES
