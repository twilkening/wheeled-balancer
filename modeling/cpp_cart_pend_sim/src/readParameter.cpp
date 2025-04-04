#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Helper function to read a parameter from a line
bool readParameter(std::ifstream& file, double& param, const std::string paramName) {
    std::string line;
    if (std::getline(file, line)) {
        std::stringstream ss(line.substr(0, line.find("//")));
        ss >> param;
        if (ss.fail()) {
            std::cerr << "Error reading " << paramName << std::endl;
            return false;
        }
        return true;
    } else {
        std::cerr << "Error reading " << paramName << std::endl;
        return false;
    }
}