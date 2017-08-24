#pragma once
#include <iostream>
#include <ctime>
#include <string>

inline std::string time_stamp(){
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    std::time (&rawtime);
    timeinfo = std::localtime(&rawtime);
    
    std::strftime(buffer, sizeof(buffer),"%Y-%m-%d_%I:%M:%S", timeinfo);
    return {buffer};
}
