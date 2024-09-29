#include "atmosphere.hpp"

#include "earth.hpp"
#include "input.hpp"

#include <fstream>

const float Air::CP_TABLE[32] = 
{
    1005,
    1006,
    1007,
    1009,
    1012,
    1015,
    1018,
    1022,
    1024,
    1031,
    1036,
    1041,
    1046,
    1051,
    1057,
    1063,
    1069,
    1075,
    1080,
    1086,
    1092,
    1098,
    1104,
    1109,
    1115,
    1120,
    1126,
    1131,
    1136,
    1141,
    1146,
    1150
};

const float Air::ENTHALPY_TABLE[32] = 
{
0,
25137.5,
50300,
75500,
100762.5,
126100,
151512.5,
177012.5,
202587.5,
228275,
254112.5,
280075,
306162.5,
332375,
358725,
385225,
411875,
438675,
465612.5,
492687.5,
519912.5,
547287.5,
574812.5,
602475,
630275,
658212.5,
686287.5,
714500,
742837.5,
771300,
799887.5,
828587.5
};

double Air::isentropic_area_ratio(double mach)
{
    auto tmp = 1.2/(1.0 + 0.2*mach*mach);
    return mach*tmp*tmp*tmp;
}

double Air::isentropic_area_ratio(double mach, double gamma)
{
    auto g1 = (gamma + 1.0)*0.5;
    auto g2 = (gamma - 1.0);
    auto tmp = 1.0 + 0.5*g2*mach*mach;
    return mach*pow(g1/tmp,g1/g2);
}

double Air::mach_area_ratio(double area_ratio, double gamma, double mach_guess)
{
    auto g1 = (gamma + 1.0)*0.5;
    auto g2 = (gamma - 1.0);
    auto exp = g1/g2;
    auto g3 = g2*0.5;
    auto g4 = 1.0/g1;
    for(int iter = 0; iter < 20; iter++)
    {
        auto M2 = mach_guess*mach_guess;
        auto tmp = 1.0 + g3*M2;
        auto dM = tmp*(area_ratio*pow(tmp*g4, exp) - mach_guess)/(tmp - g1*M2);
        mach_guess += dM;
        if(fabs(dM) < 1e-6)
        {
            break;
        }
    }
    return mach_guess;
}

double Air::supersonic_mach_area_ratio(double area_ratio, double gamma)
{
    if(area_ratio > 0.25000001)
    {
        static constexpr double MACH[16] = {
            1.0,
            1.26592402979173,
            1.39299837292862,
            1.50044885442846,
            1.59970844014855,
            1.69545040476657,
            1.79035559520279,
            1.88634322285565,
            1.98504910077651,
            2.08807941000668,
            2.19719812165219,
            2.31451972255102,
            2.44276484552305,
            2.58565970347403,
            2.74863500620618,
            2.94017916931260
        };

        auto frac = (1.0 - area_ratio)*20.0; // 20 is 1/0.05 (0.05 is increment of area_ratio in table)
        auto idx = static_cast<unsigned>(frac);
        frac -= static_cast<double>(idx);
        auto mach = MACH[idx] + frac*(MACH[idx + 1] - MACH[idx]);
        return mach_area_ratio(area_ratio, gamma, mach);
    }
    else
    {
        double mach = std::min(0.05/area_ratio + 3.85 - 4.5*area_ratio, 4.0 + 0.25/sqrt(area_ratio));
        return mach_area_ratio(area_ratio, gamma, mach);
    }  
}

double Air::subsonic_mach_area_ratio(double area_ratio, double gamma)
{
    auto dx = 1.0 - area_ratio; // decent guess ( approximately inverse parabola)
    return mach_area_ratio(area_ratio, gamma, 1.0 - sqrt(dx));
}

Air::ShockQuantities Air::normal_shock(double mach, double gamma)
{
    Air::ShockQuantities output;
    auto M2 = mach*mach;
    auto g1 = gamma - 1.0;
    auto g2 = gamma + 1.0;

    auto f1 = 2*gamma*M2 - g1;
    auto f2 = g1*M2 + 2;

    output.pressure_ratio = f1/g2;
    output.density_ratio = g2*M2/f2;
    output.mach_sq = f2/f1;

    return output;
}

Air::ObliqueShockQuantities Air::oblique_shock(double mach, double gamma, double shockAngle)
{
    Air::ObliqueShockQuantities output;
    double cTheta = cos(shockAngle);
    double sTheta = sin(shockAngle);
    double m1 = cTheta*mach;
    output.shock = normal_shock(m1, gamma);
    double mTangent = sTheta*mach;
    output.deflection_angle = atan(mTangent/sqrt(output.shock.mach_sq));
    return output;
}

float Air::cp_air(double temperature)
{
    // assume 1 bar and temperature > 220 K temperature < 1050 K
    if(temperature < CP_TABLE_KELVIN_START)
    {
        return CP_TABLE[0];
    }
    if(temperature >= 1050)
    {
        return CP_TABLE[31] + (temperature - 1050.0)*0.0001375;
    }
    float f = (temperature - CP_TABLE_KELVIN_START)*CP_TABLE_KELVIN_SCALE;
    auto idx = static_cast<unsigned>(f);
    return CP_TABLE[idx] + (f - static_cast<float>(idx))*(CP_TABLE[idx+1] - CP_TABLE[idx]);
}

float Air::enthalpy_air(double temperature)
{
    // assume 1 bar and temperature > 220 K temperature < 1050 K
    if(temperature < CP_TABLE_KELVIN_START)
    {
        return -CP_TABLE[0]*temperature;
    }
    if(temperature >= 1050)
    {
        double dt = (temperature - 1050.0);
        return ENTHALPY_TABLE[31] + dt*(CP_TABLE[31] + dt*6.875E-05);
    }
    float f = (temperature - CP_TABLE_KELVIN_START)*CP_TABLE_KELVIN_SCALE;
    auto idx = static_cast<unsigned>(f);
    return ENTHALPY_TABLE[idx] + (f - static_cast<float>(idx))*(ENTHALPY_TABLE[idx+1] - ENTHALPY_TABLE[idx]);
}


const Air Atmosphere::VACUUM = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0};

const Air Atmosphere::SEA_LEVEL = {101325.0, 1.225, 1.0/340.294, 288.15, 1.608e-5, 1.4, 287.0528};

const double AtmosphereLinearTable::US_1976_HEIGHTS[10] = 
{
    0,
    11000,
    20000,
    32000,
    47000,
    51000,
    71000,
    84852,
    90000,
    105000
};

const double AtmosphereLinearTable::US_1976_LAPSE_RATE[10] = 
{
    -0.0065,
    0,
    0.001,
    0.0028,
    0,
    -0.0028,
    -0.002,
    0,
    0.004,
    0.005
};

const double AtmosphereLinearTable::NRLMSISE_HEIGHTS[10] = 
{
    0
};

const double AtmosphereLinearTable::NRLMSISE_LAPSE_RATE[10] = 
{
    0
};


std::vector<Air> AtmosphereLinearTable::compute_linear_factors(const std::vector<Air>& table, 
    double height_increment)
{
    std::vector<Air> linear_factors;
    linear_factors.resize(table.size());
    double inv_height_increment = 1.0/height_increment;

    for(auto i = 1u; i < table.size(); i++)
    {
        auto* dv = linear_factors[i-1].data();
        const auto* v0 = table[i-1].data();
        const auto* v1 = table[i].data();
        for(auto j = 0u; j < Air::NVALUES; j++)
        {
            dv[j] = (v1[j] - v0[j])*inv_height_increment;
        }
    }
    return linear_factors;
}

AtmosphereLinearTable::AtmosphereLinearTable(std::vector<Air> table, double min_height, 
    double height_increment) : _table(table), _linear_factors(compute_linear_factors(table, height_increment)), 
    _min_height(min_height), _max_height(min_height + height_increment*table.size()), _height_increment(height_increment), 
    _inv_height_increment(1.0/_height_increment) {}

AtmosphereLinearTable AtmosphereLinearTable::create(double scale_height, double reference_pressure, double reference_temperature,
    double max_height, double height_increment, double min_height)
{
    constexpr double MW_DRY = 0.0289652;

    Air air;
    air.temperature = reference_temperature;
    air.dynamic_viscosity = Air::dynamic_viscosity_sutherland(reference_temperature);
    air.gamma = 1.4;
    air.specific_gas_constant = Air::GAS_CONSTANT / MW_DRY;
    air.inv_sound_speed = 1.0/sqrt(air.gamma*air.temperature*air.specific_gas_constant);

    double R = 6371000 + scale_height;
    double g0 = Earth::MU/(R*R);
    const double scale_factor = -g0/(air.specific_gas_constant*reference_temperature);

    std::vector<Air> table;
    table.reserve(static_cast<size_t>((max_height - min_height)/height_increment) + 1);
    double h = min_height;
    while(h < max_height)
    {
        Air input = air;
        input.pressure = reference_pressure*exp((scale_height - h)*scale_factor);
        input.density = input.pressure/(air.specific_gas_constant*reference_temperature);
        table.push_back(input);
        h += height_increment;
    }

    return AtmosphereLinearTable(table, min_height, height_increment);
}

AtmosphereLinearTable AtmosphereLinearTable::create(STD_ATMOSPHERES, double height_increment, double ground_temperature, double ground_pressure, double g0)
{
    // TODO: use STD_ATMOSPHERES
    constexpr double SPACE_LINE = 100000;

    std::vector<Air> table;
    table.reserve(static_cast<size_t>(SPACE_LINE/height_increment) + 1);

    double h = 0;
    double h_geo_potential = 0;
    double tref = ground_temperature;
    double pref = ground_pressure;
    int atmIdx = 0;
    bool isothermal = false;
    double avg_gas_const = 0;
    double avg_gamma = 0;
    int layer_count = 0;
    while(h < SPACE_LINE)
    {
        Air air;

        h_geo_potential = geometric2geopotential(h);

        if(h_geo_potential < US_1976_HEIGHTS[4])
        {
            air.specific_gas_constant = Air::GAS_CONSTANT / Air::DRY_AIR_MW;
            air.gamma = 1.4;
        }
        else
        {
            double factor = (h_geo_potential - US_1976_HEIGHTS[4])/SPACE_LINE;
            factor *= factor;
            air.specific_gas_constant = Air::GAS_CONSTANT / (Air::DRY_AIR_MW*(1 - factor) + factor*0.021);
            air.gamma = 1.4*(1 - factor) + factor*1.63;
        }

        if(atmIdx < 10 && h_geo_potential > US_1976_HEIGHTS[atmIdx + 1])
        {
            avg_gas_const *= 1.0/layer_count;
            avg_gamma *= 1.0/layer_count;

            pref *= (isothermal ? Atmosphere::pressure_ratio_isothermal(US_1976_HEIGHTS[atmIdx + 1], US_1976_HEIGHTS[atmIdx], tref, avg_gas_const, g0) : 
                Atmosphere::pressure_ratio_linearthermal(US_1976_HEIGHTS[atmIdx + 1], US_1976_HEIGHTS[atmIdx], tref, US_1976_LAPSE_RATE[atmIdx], avg_gas_const, g0));

            tref += US_1976_LAPSE_RATE[atmIdx]*(US_1976_HEIGHTS[atmIdx + 1] - US_1976_HEIGHTS[atmIdx]);

            atmIdx++;
            isothermal = abs(US_1976_LAPSE_RATE[atmIdx]) < 1e-6;
            avg_gas_const = 0;
            avg_gamma = 0;
            layer_count = 0;
        }

        avg_gas_const += air.specific_gas_constant;
        avg_gamma += air.gamma;
        layer_count++;

        if(isothermal)
        {
            air.temperature = tref;
            air.pressure = pref*Atmosphere::pressure_ratio_isothermal(h_geo_potential, US_1976_HEIGHTS[atmIdx], tref, avg_gas_const/layer_count, g0);
        }
        else
        {
            air.temperature = tref + (h_geo_potential - US_1976_HEIGHTS[atmIdx])*US_1976_LAPSE_RATE[atmIdx];
            air.pressure = pref*Atmosphere::pressure_ratio_linearthermal(h_geo_potential, US_1976_HEIGHTS[atmIdx], tref, US_1976_LAPSE_RATE[atmIdx], avg_gas_const/layer_count, g0);
        }
        air.inv_sound_speed = 1.0/sqrt(air.gamma*air.specific_gas_constant*air.temperature);
        air.density = air.pressure/(air.specific_gas_constant*air.temperature);

        table.push_back(air);

        h += height_increment;
        
    }

    return AtmosphereLinearTable(table, 0.0, height_increment);
}

AtmosphereLinearTable AtmosphereLinearTable::create(const std::string& filename)
{
    std::ifstream file(filename);
    if(!file.is_open())
    {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    std::string line;
    if(!std::getline(file, line)) 
    {
        throw std::runtime_error("Invalid file format: " + filename);
    }
    // skipped header
    std::vector<Air> table;
    while(std::getline(file, line))
    {
        std::vector<std::string> values = input::split(line);
        if(values.size()!= Air::NVALUES) 
        {
            break;
        }
        Air air;
        double h = std::stod(values[0]);
        double temperature = std::stod(values[1]);
    }
    // TODO

    return AtmosphereLinearTable(table, 1.0, 1.0);
}

void AtmosphereLinearTable::saveAsTable(std::string filename, double height_increment) const
{
    std::ofstream file(filename);
    if(!file.is_open())
    {
        std::cerr << "Error: Unable to open file: " << filename << std::endl;
        return;
    }
    double h = _min_height;
    Air air;
    while(h < _max_height)
    {
        this->set_air(h, air);
        file << h << " " << air.temperature << " " << air.pressure << " " 
            << air.density << " " << 1.0/air.inv_sound_speed << " " 
            << air.dynamic_viscosity << " " << AtmosphereLinearTable::geometric2geopotential(h) << "\n";
        h += height_increment;
    }
}

void AtmosphereLinearTable::set_air(double height, Air& air) const
{
    if(height > _max_height)
    {
        air = VACUUM;
        return;
    }
    height = std::max(height, 0.0);

    auto delta = height*_inv_height_increment;
    const auto idx = static_cast<unsigned>(delta);
    delta -= idx;

    auto* data = air.data();
    const auto* v = _table[idx].data();
    const auto* dv = _linear_factors[idx].data();
    for(auto i = 0u; i < Air::NVALUES; i++)
    {
        data[i] = v[i] + delta*dv[i];
    }
}