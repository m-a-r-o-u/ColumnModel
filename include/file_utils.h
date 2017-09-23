#pragma once
#include <fstream>
#include <string>
#include <iostream>

inline std::ifstream loadfile(std::string fname){
    std::ifstream infile (fname);
    if (!infile.good()) {
        std::cout << "can't find: " << fname << std::endl;
    }
    return infile;
}
