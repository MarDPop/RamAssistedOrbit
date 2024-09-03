#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <array>

namespace Constants
{
    inline constexpr double EARTH_ROTATION_RATE = 7.2921150e-5;
    
};

class GForce
{
public:
    static constexpr double G = 9.806;
    static inline constexpr std::array<double, 7> MAX_G_TIME = {0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 60.0};
    static inline constexpr std::array<double, 7> MAX_G_BACK = {32.0, 27.0, 20.0, 15.0, 11.0, 9.0, 6.0};
    static inline constexpr std::array<double, 7> MAX_G_UP = {20.0, 14.0, 10.0, 8.0, 2.8, 2.6, 2.5};
    static inline constexpr std::array<double, 7> MAX_G_DOWN = {10.0, 7.0, 5.0, 4.0, 1.3, 1.2, 1.0};
};

#endif