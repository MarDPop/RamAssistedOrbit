#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <array>
#include <vector>
#include <cmath>
#include <string>

struct Air
{
    /* STATIC PROPERTIES */
    /**
     * Universal Gas Constant in SI units (J / K mol) (exact)
     */
    static constexpr double GAS_CONSTANT = 8.31446261815324; 
    /**
     * Dry Air Molecular Weight in SI units (kg / mol)
     */
    static constexpr double DRY_AIR_MW = 0.028964917; 

    /**
     * Definition of absolute zero (ie conversion from Kelvin to Celsius)
     * (Exact by International Bureau of Weights and Measures in 2007)
     */
    static constexpr double ABSOLUTE_ZERO_CELSIUS = -273.15;

    static constexpr unsigned CP_TABLE_ENTRIES = 32;
    /**
     * Kelvin start of the CP and ENTHALPY tables
     */
    static constexpr double CP_TABLE_KELVIN_START = 275;

    /**
     * inverse of the table spacing ie (1/delta) where 1/25 = 0.04
     */
    static constexpr double CP_TABLE_KELVIN_SCALE = 0.04;

    /**
     * Values of CP for dry air from CP_TABLE_KELVIN_START to CP_TABLE_KELVIN_SCALE + CP_TABLE_ENTRIES/CP_TABLE_KELVIN_SCALE
     */
    static const float CP_TABLE[CP_TABLE_ENTRIES]; // Joules
    static const float ENTHALPY_TABLE[CP_TABLE_ENTRIES]; // Joules

    /* STATIC METHODS */

    // ISENTROPIC FLOW STUFF

    static double choked_flow(double total_pressure, double total_temperature, double molecular_weight, double gamma) 
    {
        double g1 = (gamma + 1.0)*0.5;
        return total_pressure*sqrt(molecular_weight*gamma/(Air::GAS_CONSTANT*total_temperature))*pow(g1, g1/(1.0 - gamma));
    }

    /**
     * Mass flow through normalized area at throat (kg/ s m2)
     * @param total_pressure
     * @param total_temperature
     */
    static double choked_flow(double total_pressure, double total_temperature) 
    {
        constexpr double AIR_CONSTANT = 0.040414899585753; // ((1.4 + 1)*0.5)^(-(1.4 + 1)*0.5/(1.4 - 1))*sqrt(1.4/ 287.05)
        return AIR_CONSTANT*total_pressure/sqrt(total_temperature);
    }

    /**
     * returns the ratio Tt/T where Tt is total temperature and T is temperature for AIR
     * @param mach
     */
    static double isentropic_temperature_ratio(double mach)
    {
        return 1.0 + 0.2*mach*mach;
    }

    /**
     * returns the ratio Pt/P where Pt is total temperature and T is temperature for AIR
     * @param mach
     */
    static double isentropic_pressure_ratio(double mach)
    {
        double tmp = isentropic_temperature_ratio(mach);
        return tmp*tmp*tmp*sqrt(tmp);
    }

    /**
     * returns the ratio Pt/P where Pt is total temperature and T is temperature for AIR
     * @param mach
     */
    static double isentropic_density_ratio(double mach)
    {
        double tmp = isentropic_temperature_ratio(mach);
        return tmp*tmp*sqrt(tmp);
    }

    /**
     * returns the ratio Astar/A where Astar is the critical area and A is the area to evaluate
     * @param mach
     */
    static double isentropic_area_ratio(double mach);

    /**
     * returns the ratio Tt/T where Tt is total temperature and T is temperature for any gamma
     * @param mach
     * @param gamma
     */
    static double isentropic_temperature_ratio(double mach, double gamma)
    {
        return 1.0 + 0.5*(gamma - 1.0)*mach*mach;
    }

    /**
     * returns the ratio Pt/P where Pt is total pressure and P is static pressure for any gamma
     * @param mach
     * @param gamma
     */
    static double isentropic_pressure_ratio(double mach, double gamma)
    {
        return pow(isentropic_temperature_ratio(mach, gamma), gamma/(gamma-1.0));
    }

    /**
     * returns the ratio Pt/P where Pt is total pressure and T is pressure for any gamma
     * @param pressure_ratio Pt/P
     * @param gamma
     */
    static double mach_from_isentropic_pressure_ratio(double pressure_ratio, double gamma)
    {
        const double g1 = gamma - 1;
        return sqrt( 2.0/g1*(pow(pressure_ratio, g1/gamma) - 1.0));
    }

    /**
     * returns the ratio Astar/A where Astar is the critical area and A is the area to evaluate
     * @param mach
     * @param gamma ratio of specific heats
     */
    static double isentropic_area_ratio(double mach, double gamma);

    /**
     * Helper function for supersonic_mach_area_ratio and subsonic_mach_area_ratio
     * Performs Newton iteration to compute the mach associated for a given area ratio given an initial guess
     * @param area_ratio area ratio you wish to compute the mach number for (Astar/A)
     * @param gamma ratio of specific heats
     * @param mach_guess initial guess (must be relatively close, chance that a supersonic guess may become subsonic if not close)
     * @return mach number for area ratio
     */
    static double mach_area_ratio(double area_ratio, double gamma, double mach_guess);

    /**
     * Get the mach number for a given area ratio (only use for supersonic)
     * @param area_ratio area ratio you wish to compute the mach number for (Astar/A)
     * @param gamma ratio of specific heats
     * @return mach number for area ratio
     */
    static double supersonic_mach_area_ratio(double area_ratio, double gamma);

    /**
     * Get the mach number for a given area ratio (only use for subsonic)
     * @param area_ratio area ratio you wish to compute the mach number for  (Astar/A)
     * @param gamma ratio of specific heats
     * @return mach number for area ratio
     */
    static double subsonic_mach_area_ratio(double area_ratio, double gamma);

    // SHOCK STUFF
    struct ShockQuantities
    {
        double pressure_ratio; // pressure after shock / pressure before shock
        double density_ratio; // density after shock / density before shock
        double mach_sq; // mach squared after shock
    };

    static ShockQuantities normal_shock(double mach, double gamma);

    static double normal_shock_temperature_ratio(double pressure_ratio, double density_ratio)
    {
        return pressure_ratio/density_ratio;
    }

    static double normal_shock_total_pressure_ratio(double pressure_ratio, double density_ratio, double gamma)
    {
        const double ex = 1.0/(gamma - 1.0);
        return pow(density_ratio, ex*gamma)*pow(pressure_ratio, -ex);
    }

    struct ObliqueShockQuantities
    {
        ShockQuantities shock;
        double deflection_angle;
    };

    static ObliqueShockQuantities oblique_shock(double mach, double gamma, double shockAngle);

    /**
     * Calculates the isothermal pressure ratio from base geopotential height
     * @param kelvin
     * @return dynamic pressure in Pa s
     */
    static double dynamic_viscosity_sutherland(double kelvin)
    {
        return (5.058e-8 - 1.295e-11 * kelvin) * kelvin + 4.371e-6;
    }

    /**
     * Gets the partial pressure of water at temperature
     * @param kelvin
     * @return pressure in Pa
     */
    static double partial_pressure_H20(double kelvin)
    {
        double celsius = kelvin + ABSOLUTE_ZERO_CELSIUS;
        return 610.76*exp(17.27*celsius/(celsius+237.5));
    }

    /**
     * Retrieves the linearly interpolated value of constant pressure specific heat capacity of dry air
     * in Joules
     * @param temperature in kelvin
     */
    static float cp_air(double temperature);

    /**
     * Retrieves the linearly interpolated value of the enthalpy *increase* of dry air
     * in Joules from CP_TABLE_KELVIN_START (275 K)
     * @param temperature in kelvin
     */
    static float enthalpy_air(double temperature);

    /**
     * Number of doubles in the Air struct, (to treat as array)
     */
    static constexpr unsigned NVALUES = 7;

    /* OBJECT PROPERTIES */
    double pressure;
    double density;
    double inv_sound_speed;
    double temperature;
    double dynamic_viscosity;
    double gamma;
    double specific_gas_constant;

    Air(){}

    Air(double* data) : pressure(data[0]), density(data[1]), inv_sound_speed(data[2]), temperature(data[3]), 
        dynamic_viscosity(data[4]), gamma(data[5]), specific_gas_constant(data[6]) {}

    Air(double pressure_, double density_, double inv_sound_speed_, double temperature_, double dynamic_viscosity_,
        double gamma_, double specific_gas_constant_) : pressure(pressure_), density(density_), inv_sound_speed(inv_sound_speed_), 
        temperature(temperature_), dynamic_viscosity(dynamic_viscosity_), gamma(gamma_), 
        specific_gas_constant(specific_gas_constant_) {}

    Air(double pressure_, double temperature_) : pressure(pressure_), specific_gas_constant(Air::GAS_CONSTANT/Air::DRY_AIR_MW),
        density(pressure/(specific_gas_constant*temperature)), 
        inv_sound_speed(1.0/sqrt(gamma*specific_gas_constant*temperature)),
        temperature(temperature_), dynamic_viscosity(dynamic_viscosity_sutherland(temperature)), gamma(1.4) {}

    Air(double pressure_, double temperature_, double gamma_ = 1.4, double mw = 0.028137, 
        double dynamic_viscosity_ = 1.61e-5) : pressure(pressure_), specific_gas_constant(Air::GAS_CONSTANT/mw), 
        density(pressure_/(specific_gas_constant*temperature_)), inv_sound_speed(1.0/sqrt(gamma*specific_gas_constant*temperature)),
        temperature(temperature_), gamma(gamma_) {}

    /* OBJECT METHODS */
    double* data()
    {
        return &pressure;
    }

    const double* data() const
    {
        return &pressure;
    }

};

struct Atmosphere 
{
    static const Air VACUUM;

    static const Air SEA_LEVEL;

    /**
     * Calculates the isothermal pressure ratio from base geopotential height
     * @param z
     * @param R0
     * @return
     */
    static double geometric2geopotential(double z, const double R0 = 6371)
    {
        return R0 * z / (R0 + z);
    }

    /**
     * Calculates the isothermal pressure ratio from base geopotential height
     * @param H
     * @param R0
     * @return
     */
    static double geopotential2geometric(double H, const double R0 = 6371)
    {
        return R0 * H / (R0 - H);
    }

    /**
     * calculates the ratio of static pressures using the barometric formula (assuming constant temperature)
     * @param h height to evaluate (meters)
     * @param href reference height for which g0, and tref are valid (meters)
     * @param tref reference temperature (Kelvin)
     * @param R Gas constant
     * @param g0 gravity constant (meters/s2)
     */
    static double pressure_ratio_isothermal(double h, double href, double tref, double R, double g0 = 9.806)
    {
        return exp(g0*(href - h)/(R*tref));
    }

    /**
     * calculates the ratio of static pressures using the barometric formula (assuming linear lapse Rate)
     * @param h height to evaluate (meters)
     * @param href (reference height for which g0, and tref are valid)
     * @param tref reference temperature (Kelvin)
     * @param lapseRate change in temperature vs altitude (Kelvin/meter)
     * @param R Gas constant (meters)
     * @param g0 gravity constant (meters/s2)
     */
    static double pressure_ratio_linearthermal(double h, double href, double tref, double lapseRate, double R, double g0 = 9.806)
    {
        return pow(1.0 - lapseRate*(h - href)/tref,g0/(lapseRate*R));
    }

    virtual void set_air(double height, Air& air) const = 0;

    virtual double get_max_height() const = 0;
};

class AtmosphereLinearTable : public Atmosphere
{
    const std::vector<Air> _table;

    const std::vector<Air> _linear_factors;

    double _min_height;

    double _max_height;

    double _height_increment;

    double _inv_height_increment;

    static std::vector<Air> compute_linear_factors(const std::vector<Air>& table,
        double height_increment);

public:

    static const double US_1976_HEIGHTS[10];
    static const double US_1976_LAPSE_RATE[10];
    static const double NRLMSISE_HEIGHTS[10];
    static const double NRLMSISE_LAPSE_RATE[10];

    enum class STD_ATMOSPHERES
    {
        US_1976,
        ICAO,
        NRLMSISE
    };

    AtmosphereLinearTable(std::vector<Air> table, double min_height, double height_increment);

    /**
     * @brief Create atmopshere table from arbitrary scale height, reference pressure, reference temperature
     * 
     * @param scale_height 
     * @param reference_pressure 
     * @param reference_temperature 
     * @param max_height 
     * @param height_increment 
     * @return AtmosphereLinearTable 
     */
    static AtmosphereLinearTable create(double scale_height, double reference_pressure, double reference_temperature,
        double max_height, double height_increment, double min_height = 0.0);

    /**
     * @brief Creates an table from general US 1976 Standard Atmosphere
     * 
     * @param sl_tmp 
     * @param p0 
     * @param g0 
     * @return AtmosphereLinearTable 
     */
    static AtmosphereLinearTable create(double sl_tmp = 288.15, double p0 = 101325, double g0 = 9.80665);

    static AtmosphereLinearTable create(STD_ATMOSPHERES atm, double height_increment, 
        double ground_temperature = 288.15, double ground_pressure = 101325.0, double g0 = 9.806);

    static AtmosphereLinearTable create(const std::string& file);

    double get_max_height() const override
    {
        return _max_height;
    }

    double get_height_increment() const
    {
        return _height_increment;
    }

    void set_air(double height, Air& air) const override;
};

#endif