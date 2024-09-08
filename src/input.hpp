#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace input {

    inline std::vector<std::string> split(std::string str) {
        std::vector<std::string> tokens;
        std::istringstream iss(str);
        std::string token;
        while(iss >> token) {
            tokens.push_back(token);
        }
        return tokens;
    }

    inline std::vector<std::string> split(std::string str, char delimiter) {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(str);
        while (std::getline(tokenStream, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }
}