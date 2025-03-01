#include "constants.hpp"
#include "functions.hpp"
#include "propulsion.hpp"

#ifdef DEBUG
    #include <iostream>
#endif

ThrusterTabulated::ThrusterTabulated(const std::vector<double>& pressures,
    const std::vector<double>& temperatures, const std::vector<double>& machs,
    const std::vector<double>& throttles, const std::vector<double>& beta_angle,
    const std::vector<double>& alpha_angle, const std::vector<double>& thrust_values,
    double min_mass_rate, double max_mass_rate) : 
        _pressures(pressures), _dpressures(functions::diff(pressures)), 
        _temperatures(temperatures), _dtemperatures(functions::diff(temperatures)), 
        _machs(machs), _dmachs(functions::diff(machs)),
        _throttles(throttles), _dthrottles(functions::diff(throttles)),
        _beta_angle(beta_angle), _dbeta_angle(functions::diff(beta_angle)),
        _alpha_angle(alpha_angle), _dalpha_angle(functions::diff(alpha_angle)),
        _thrust_values(thrust_values), _min_mass_rate(min_mass_rate), _max_mass_rate(max_mass_rate) 
{
    layer_offset[5] = 1;
    layer_offset[4] = beta_angle.size();
    layer_offset[3] = layer_offset[3]*alpha_angle.size();
    layer_offset[2] = layer_offset[2]*throttles.size();
    layer_offset[1] = layer_offset[1]*machs.size();
    layer_offset[0] = layer_offset[0]*temperatures.size();

    unsigned total_nodes = layer_offset[0]*pressures.size();
    if(thrust_values.size() != total_nodes)
    {
        throw std::runtime_error("Thrust values size does not match the number of nodes");
    }
}

void ThrusterTabulated::update_thrust(const Air& air, const AeroQuantities& aero)  
{
    const unsigned LAYERS = 6;
    double factors[LAYERS];
    unsigned index[LAYERS];
    index[0] = std::lower_bound(_pressures.begin(), _pressures.end(), air.pressure) - _pressures.begin();
    index[1] = std::lower_bound(_temperatures.begin(), _temperatures.end(), air.temperature) - _temperatures.begin();
    index[2] = std::lower_bound(_machs.begin(), _machs.end(), aero.mach) - _machs.begin();
    index[3] = std::lower_bound(_throttles.begin(), _throttles.end(), _throttle) - _throttles.begin();
    index[4] = std::lower_bound(_beta_angle.begin(), _beta_angle.end(), fabs(aero.beta_angle)) - _beta_angle.begin();
    index[5] = std::lower_bound(_alpha_angle.begin(), _alpha_angle.end(), aero.alpha_angle) - _alpha_angle.begin();

    factors[0] = (air.pressure - _pressures[index[0]])*_dpressures[index[0]];   
    factors[1] = (air.temperature - _temperatures[index[1]])*_dtemperatures[index[1]];
    factors[2] = (aero.mach - _machs[index[2]])*_dmachs[index[2]];
    factors[3] = (_throttle - _throttles[index[3]])*_dthrottles[index[3]];
    factors[4] = (fabs(aero.beta_angle) - _beta_angle[index[4]])*_dbeta_angle[index[4]];
    factors[5] = (aero.alpha_angle - _alpha_angle[index[5]])*_dalpha_angle[index[5]];

    constexpr unsigned BRANCHES = 2 << LAYERS;
    unsigned thrust_idx[BRANCHES];
    const unsigned top_idx = index[0]*layer_offset[0] + index[1]*layer_offset[1] + index[2]*layer_offset[2] 
        + index[3]*layer_offset[3] + index[4]*layer_offset[4] + index[5]*layer_offset[5];

    std::fill_n(thrust_idx, BRANCHES, top_idx);

    for(unsigned i = 32; i < BRANCHES; i++)
    {
        thrust_idx[i] += layer_offset[0];
    }

    for(unsigned i = 16; i < BRANCHES; i+=32) 
    {
        for(unsigned j = 0; j < 16; j++)
        {
            thrust_idx[i + j] += layer_offset[1];
        }
    }

    for(unsigned i = 8; i < BRANCHES; i+=16) 
    {
        for(unsigned j = 0; j < 8; j++)
        {
            thrust_idx[i + j] += layer_offset[2];
        }
    }

    for(unsigned i = 4; i < BRANCHES; i+=8) 
    {
        for(unsigned j = 0; j < 4; j++)
        {
            thrust_idx[i + j] += layer_offset[3];
        }
    }

    for(unsigned i = 2; i < BRANCHES; i+=4) 
    {
        for(unsigned j = 0; j < 2; j++)
        {
            thrust_idx[i + j] += layer_offset[4];
        }
    }

    for(unsigned i = 0; i < BRANCHES; i+=2)
    {
        thrust_idx[i]++;
    }

    double thrust_layer5[32];
    for(unsigned i = 0; i < 32; i++)
    {
        unsigned idx = i*2;
        thrust_layer5[i] = _thrust_values[thrust_idx[idx]] + factors[5]*(_thrust_values[thrust_idx[idx + 1]] - _thrust_values[thrust_idx[idx]]);
    }

    double thrust_layer4[16];
    for(unsigned i = 0; i < 16; i++)
    {
        unsigned idx = i*2;
        thrust_layer4[i] = thrust_layer5[idx] + factors[4]*(thrust_layer5[idx + 1] - thrust_layer5[idx]);
    }

    double thrust_layer3[8];
    for(unsigned i = 0; i < 8; i++)
    {
        unsigned idx = i*2;
        thrust_layer3[i] = thrust_layer4[idx] + factors[3]*(thrust_layer4[idx + 1] - thrust_layer4[idx]);
    }   
    
    double thrust_layer2[4];
    for(unsigned i = 0; i < 4; i++)
    {
        unsigned idx = i*2;
        thrust_layer2[i] = thrust_layer3[idx] + factors[2]*(thrust_layer3[idx + 1] - thrust_layer3[idx]);
    }

    double thrust_layer1[2];
    for(unsigned i = 0; i < 2; i++)
    {
        unsigned idx = i*2;
        thrust_layer1[i] = thrust_layer2[idx] + factors[1]*(thrust_layer2[idx + 1] - thrust_layer2[idx]);
    }

    _thrust = thrust_layer1[0] + factors[0]*(thrust_layer1[1] - thrust_layer1[0]);  
}

double shock_exit_pressure_ratio(double exit_mach)
{
    double p_ratio_exit = Air::isentropic_pressure_ratio(exit_mach);
    double p_ratio_shock = 1.1666666666666666666666666666*exit_mach*exit_mach - 0.16666666666666666666666;
    return p_ratio_shock/p_ratio_exit;
}

RamjetSimpleModel::RamjetSimpleModel(double throat_area, 
        double min_exit_area, 
        double max_exit_area,
        double diffusor_efficiency, 
        double burner_pressure_ratio, 
        double burner_efficiency,
        double nozzle_efficiency, 
        double dry_mass, 
        double max_mass_rate,
        double combustion_temperature,
        double fuel_heating_value, 
        double max_fuel_air_ratio, 
        double gamma_combustion, 
        double mw_combustion) :
    _throat_area(throat_area), 
    _min_exit_area(min_exit_area),
    _max_exit_area(max_exit_area),
    _diffusor_efficiency(diffusor_efficiency),
    _burner_pressure_ratio(burner_pressure_ratio),
    _burner_efficiency(burner_efficiency),
    _nozzle_efficiency(nozzle_efficiency), 
    _max_mass_rate(max_mass_rate),
    _combustion_temperature(combustion_temperature), 
    _fuel_heating_value(fuel_heating_value), 
    _max_fuel_air_ratio(max_fuel_air_ratio), 
    _gamma_combustion(gamma_combustion),
    _cp_4(Air::GAS_CONSTANT*gamma_combustion/(gamma_combustion - 1.0)/mw_combustion),
    _specific_enthalpy_4(combustion_temperature*_cp_4),
    _max_exit_mach(Air::supersonic_mach_area_ratio(_max_exit_area, _gamma_combustion)),
    _min_exit_mach(Air::supersonic_mach_area_ratio(_min_exit_area, _gamma_combustion)),
    _max_exit_area_temperature_ratio(Air::isentropic_temperature_ratio(_max_exit_mach, _gamma_combustion)),
    _min_exit_area_temperature_ratio(Air::isentropic_temperature_ratio(_min_exit_mach, _gamma_combustion)),
    _max_exit_area_pressure_ratio(Air::isentropic_pressure_ratio(_max_exit_mach, _gamma_combustion)),
    _min_exit_area_pressure_ratio(Air::isentropic_pressure_ratio(_min_exit_mach, _gamma_combustion)),
    _exit_const(Air::GAS_CONSTANT*gamma_combustion/mw_combustion),
    _nH(burner_efficiency*fuel_heating_value),
    _cp_ratio4(1.0/(_cp_4*(1.0 + fuel_heating_value))),
    _h04_max(_cp_4*combustion_temperature),
    _g1_4(2/(_gamma_combustion - 1.0)),
    _g2_4((gamma_combustion - 1.0)/_gamma_combustion)
{
    _dry_mass = dry_mass;
}

void RamjetSimpleModel::update_thrust(const Air& air, const AeroQuantities& aero) 
{
    const auto mach_const = 0.2*aero.mach*aero.mach;

    // inlet
    const auto ram_recovery_factor = RamjetSimpleModel::ram_recovery_factor(aero.mach);
    const auto beta_fake = 1.0 + _diffusor_efficiency*mach_const;
    const auto p_total_2 = air.pressure*ram_recovery_factor*(beta_fake*beta_fake*beta_fake*sqrt(beta_fake));
    constexpr double t_total_inlet_loss = 0.998;
    const auto t_total_2 = air.temperature*(1.0 + mach_const)*t_total_inlet_loss;
    
    // Get Fuel Rate
    const auto air_mass_rate = _throat_area*Air::choked_flow(p_total_2, t_total_2);

    constexpr double cp_air = 1005.0;
    double total_enthalpy_2 = cp_air*t_total_2;
    double t_total_4_stoichiometric = (_max_fuel_air_ratio*_nH + total_enthalpy_2)*_cp_ratio4;

    double t_total_4 = t_total_4_stoichiometric;
    double f = _max_fuel_air_ratio;
    if (t_total_4_stoichiometric > _combustion_temperature) 
    {
        t_total_4 = _combustion_temperature;
        f = (_h04_max - total_enthalpy_2)/(_nH - _h04_max);
    }

    if(f < 0.0) 
    {
        _mass_rate = 0.0;
        _thrust = 0.0;
        return;
    }
  
    _mass_rate = air_mass_rate*f;

    double total_pressure_4 = p_total_2*_burner_pressure_ratio;

    double p_ratio_ideal = total_pressure_4/air.pressure;

    double M_sq_exit_ideal = (pow(p_ratio_ideal,_g2_4) - 1.0)*_g1_4;
    double M_exit_ideal = sqrt(M_sq_exit_ideal);
    double A_exit_ideal = Air::isentropic_area_ratio(M_exit_ideal, _gamma_combustion);

    if (A_exit_ideal > _max_exit_area)
    {
        _A_exit = _max_exit_area;
        _M_exit = _max_exit_mach;
    }
    else if (A_exit_ideal < _min_exit_area)
    {
        _A_exit = _min_exit_area;
        _M_exit = _min_exit_mach;
    }
    else
    {
        _A_exit = A_exit_ideal;
        _M_exit = M_exit_ideal;
    }

    double beta = 1.0/Air::isentropic_temperature_ratio(_M_exit, _gamma_combustion);

    _T_exit = t_total_4*beta;
    _v_exit = sqrt(_exit_const*_T_exit)*_M_exit*_nozzle_efficiency;
    _p_exit = total_pressure_4*pow(beta, 1.0/_g2_4);
    
    double ST_ideal = (1 + f)*_v_exit - aero.airspeed;
    _thrust = air_mass_rate*ST_ideal + (_p_exit - air.pressure)*_A_exit;
}

RamjetReal::RamjetReal( double max_mass_rate, double throat_area, double nominal_exit_area, 
        double min_exit_area, double max_exit_area, double max_intake_area,
        double heating_value_fuel, double fuel_air_ratio, 
        double combustion_efficiency, double combustor_pressure_ratio, double nozzle_efficiency, 
        double nozzle_pressure_ratio, double adiabatic_efficiency) : 
    _heating_value_fuel(heating_value_fuel),
    _fuel_air_ratio(fuel_air_ratio),
    _max_mass_rate(max_mass_rate),
    _max_air_ingest(1.25*max_mass_rate/fuel_air_ratio),
    _combustion_efficiency(combustion_efficiency),
    _combustor_pressure_ratio(combustor_pressure_ratio),
    _nozzle_efficiency(nozzle_efficiency),
    _nozzle_pressure_ratio(nozzle_pressure_ratio),
    _adiabatic_efficiency(adiabatic_efficiency),
    _throat_area(throat_area),
    _max_intake_area(max_intake_area),
    _critical_mach(Air::supersonic_mach_area_ratio(throat_area/_max_intake_area, 1.4)),
    _nominal_exit_area(nominal_exit_area),
    _max_exit_area(max_exit_area),
    _max_exit_mach(Air::supersonic_mach_area_ratio(throat_area/max_exit_area, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _max_mach_exit_pressure_ratio(nozzle_pressure_ratio/Air::isentropic_pressure_ratio(_max_exit_mach, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _max_mach_exit_temperature_ratio(adiabatic_efficiency/Air::isentropic_temperature_ratio(_max_exit_mach, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _min_exit_area(min_exit_area),
    _min_exit_mach(Air::supersonic_mach_area_ratio(throat_area/min_exit_area, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _min_mach_exit_pressure_ratio(nozzle_pressure_ratio/Air::isentropic_pressure_ratio(_min_exit_mach, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _min_mach_exit_temperature_ratio(adiabatic_efficiency/Air::isentropic_temperature_ratio(_min_exit_mach, GAMMA_COMBUSTION_PRODUCTS_KEROSENE))
    {}

RamjetReal RamjetReal::create(double thrust2weight, double altitude, double exit_mach, double cruise_mach,
    double lift2drag, double mass, const double thrust_margin, const double mass_rate_margin)
{
    const double weight = mass*GForce::G;
    const double drag = weight/lift2drag;
    const double desired_thrust = drag*(1.0 + thrust_margin);

    AtmosphereLinearTable atm = AtmosphereLinearTable::create(AtmosphereLinearTable::STD_ATMOSPHERES::US_1976, 100);
    Air air;
    atm.set_air(altitude, air);
    AeroQuantities aero;
    const double cruise_speed = cruise_mach/air.inv_sound_speed;
    Eigen::Vector3d velocity(cruise_speed, 0, 0);
    Eigen::Matrix3d CS = Eigen::Matrix3d::Identity();
    aero.update(air, velocity, CS);

    constexpr double AIR_GAMMA = 1.4;
    constexpr double APPROX_EXIT_VELOCITY = 2100;
    constexpr double APPROX_FUEL_AIR = 1.1;

    const double exit_area_ratio = 1.0/Air::isentropic_area_ratio(exit_mach, GAMMA_COMBUSTION_PRODUCTS_KEROSENE);
    const double intake_area_ratio = Air::isentropic_area_ratio(cruise_mach, AIR_GAMMA);
    
    double throat_area = desired_thrust*intake_area_ratio/(air.density*cruise_speed*APPROX_FUEL_AIR*(APPROX_EXIT_VELOCITY - cruise_speed));
    
    constexpr int MAXITERATIONS = 10;
    constexpr double AREA_FRACTION = 0.05;
    constexpr double MIN_THROAT_AREA_FRACTION = 0.2;
    double exit_area = 0.0;
    double mass_rate = 0.0;
    for(int iter = 0; iter < MAXITERATIONS; iter++)
    {
        exit_area = exit_area_ratio*throat_area;
        RamjetReal ramjet(1e10, throat_area, exit_area,exit_area,exit_area, throat_area*5e3);
        ramjet.update_thrust(air, aero);
        double current_thrust = ramjet.get_thrust();

        double dA = throat_area*AREA_FRACTION;
        double throat_area1 = throat_area + dA;
        exit_area = exit_area_ratio*throat_area1;
        RamjetReal ramjet2(1e10, throat_area1, exit_area,exit_area,exit_area, throat_area*5e3);
        ramjet2.update_thrust(air, aero);
        double more_thrust = ramjet2.get_thrust();

        double dArea = dA*(current_thrust - desired_thrust)/(more_thrust - current_thrust);

        throat_area = std::max(MIN_THROAT_AREA_FRACTION*throat_area, throat_area - dArea);

        if(fabs(dArea*2) < dA)
        {
            mass_rate = ramjet.get_mass_rate();
            break;
        }
    }
    double intakeAreaIdeal = throat_area/Air::isentropic_area_ratio(cruise_mach, AIR_GAMMA);
    RamjetReal ramjet(mass_rate*mass_rate_margin, throat_area, exit_area, exit_area*0.5, exit_area*2, intakeAreaIdeal*1.2);
    ramjet._dry_mass = desired_thrust/(GForce::G*thrust2weight);
    ramjet.update_thrust(air, aero);
    return ramjet;
}

bool RamjetReal::can_turn_on(const Air& air, const AeroQuantities& aero) const 
{
    double intake_area = _throat_area/Air::isentropic_area_ratio(aero.mach);
    const auto mdot_air = aero.airspeed*air.density*std::min(intake_area,_max_intake_area)*(-aero.air_body_vector.x());
    return mdot_air < _max_air_ingest;
}

void RamjetReal::update_thrust(const Air& air, const AeroQuantities& aero) 
{    
    const auto g1 = (air.gamma + 1.0)*0.5;
    const auto g2 = (air.gamma - 1.0);
    const auto beta = 1.0 + 0.5*g2*aero.mach*aero.mach;

    #ifdef DEBUG
        std::cout << "Mach: " << aero.mach << " ";
    #endif

    // Expandable intake area to match throat
    const auto isentropic_area_ratio = pow(beta/g1,g1/g2)/aero.mach;

    // Check Intake Area
    double intake_area = std::min(_throat_area*isentropic_area_ratio, _max_intake_area); // ideal intake area
    #ifdef DEBUG
    if(_throat_area*isentropic_area_ratio > _max_intake_area)
    {
        // if required intake area too big penalize
        std::cout << "at max intake area\n";
        
    }
    #endif
    // const double shock_losses = 1.0 - (mach > _critical_mach)*(mach - _critical_mach)*0.02 - mach*0.01;

    const double diffuser_losses = std::max(0.1, 0.98 - 0.06*(aero.mach - 1.0)*sqrt(aero.mach - 1.0));

    // Get Air flowing through intake
    const auto mdot_air = aero.airspeed*air.density*intake_area*(-aero.air_body_vector.x());

    // Compute required fuel rate to match
    const double ideal_mass_rate = mdot_air*_fuel_air_ratio;
    _mass_rate = std::min(ideal_mass_rate, _max_mass_rate);
    double stoichiometry_penalty = _mass_rate/ideal_mass_rate;
    #ifdef DEBUG
    if(ideal_mass_rate > _max_mass_rate)
    {
        std::cout << "At max mass rate!\n";
    }
    std::cout << "shock losses: " << diffuser_losses << " ";
    std::cout << "intake area: " << intake_area << " ";
    std::cout << "Air mass rate: " << mdot_air << " ";
    #endif

    const auto mdot_out = _mass_rate + mdot_air;

    // Get combustion conditions
    // Get isentropic compression ( total temperature )
    const auto total_temperature = air.temperature*beta;
    double pressure_ratio = beta*beta*beta*sqrt(beta);
    const auto ptotal_combustor = air.pressure*pressure_ratio*diffuser_losses;

    constexpr double mach_combustor = 0.2;
    const auto beta_combustor = Air::isentropic_temperature_ratio(mach_combustor);

    const auto temperature_into_combustor = total_temperature / beta_combustor;
    // Add energy from burning fuel
    // https://www.engineeringtoolbox.com/air-specific-heat-capacity-d_705.html
    const auto enthalpy_rate_in_combustion_chamber = Air::enthalpy_air(temperature_into_combustor)*mdot_air;
    const auto enthalpy_rate_burn = _mass_rate*_heating_value_fuel*stoichiometry_penalty;
    const auto enthalpy_rate_out_combustion_chamber = enthalpy_rate_in_combustion_chamber + enthalpy_rate_burn;

    //https://ntrs.nasa.gov/api/citations/19740019632/downloads/19740019632.pdf

    constexpr double RG = R_GAS_COMBUSTION_PRODUCTS_KEROSENE*GAMMA_COMBUSTION_PRODUCTS_KEROSENE;
    constexpr double CP_COMBUSTION_PRODUCTS = RG/(GAMMA_COMBUSTION_PRODUCTS_KEROSENE - 1.0);
    
    const auto t_total_out_combustion_chamber = enthalpy_rate_out_combustion_chamber*_combustion_efficiency
        *beta_combustor/(mdot_out*CP_COMBUSTION_PRODUCTS);

    // Get Nozzle Conditions
    const auto max_pressure_exit = ptotal_combustor*_min_mach_exit_pressure_ratio;
    const auto min_pressure_exit = ptotal_combustor*_max_mach_exit_pressure_ratio;

    const auto pressure_exit = std::clamp(air.pressure, min_pressure_exit, max_pressure_exit);

    pressure_ratio = ptotal_combustor/pressure_exit;
    constexpr double BETA_COEF_KEROSENE = (GAMMA_COMBUSTION_PRODUCTS_KEROSENE - 1.0)*0.5;
    constexpr double INV_BETA_COEF_KEROSENE = 1.0/BETA_COEF_KEROSENE;
    const auto exit_mach_sq = INV_BETA_COEF_KEROSENE*(pow(pressure_ratio, INV_EXP_COMBUSTION_PRODUCTS_KEROSENE) - 1.0);
    
    const auto tmp = 1.0/(1.0 + BETA_COEF_KEROSENE*exit_mach_sq);
    const auto t_out = t_total_out_combustion_chamber*tmp*_adiabatic_efficiency;
    const auto v_out = sqrt(RG*t_out*exit_mach_sq);
    constexpr double mix_g1 = (GAMMA_COMBUSTION_PRODUCTS_KEROSENE + 1.0)*0.5;
    constexpr double mix_g2 = mix_g1/(GAMMA_COMBUSTION_PRODUCTS_KEROSENE - 1.0);
    const auto exit_area = _throat_area*sqrt(exit_mach_sq)*pow(mix_g1*tmp, mix_g2);

    _thrust = mdot_out*v_out*_nozzle_efficiency - mdot_air*aero.airspeed + exit_area*(pressure_exit - air.pressure); // p_intake ~= pambient

    _thrust = std::max(_thrust, 0.0);
}


Rocket::Rocket(const std::vector<double>& pressures, const std::vector<double>& thrusts, 
        double max_mass_rate, double min_mass_rate, int max_number_restarts) :
        _max_mass_rate(max_mass_rate), _min_mass_rate(min_mass_rate), 
        _dmass_rate_throttle(max_mass_rate - min_mass_rate), _inv_max_mass_rate(1.0/max_mass_rate),
        _max_number_restarts(max_number_restarts)
{
    _pressures = pressures;
    _thrusts = thrusts;
    _dthrusts.resize(pressures.size());
    for(auto i = 1u; i < pressures.size(); i++)
    {
        _dthrusts[i-1] = (thrusts[i] - thrusts[i-1]) / (pressures[i] - pressures[i-1]);
    }
}

void Rocket::update_thrust(const Air& air, const AeroQuantities&)
{
    if(air.pressure < _pressures[0])
    {
        _thrust = _thrusts[0];
    }
    else if(air.pressure >= _pressures.back())
    {
        _thrust = _thrusts.back();
    }
    else
    {
        for(auto idx = 0u; idx < _pressures.size(); idx++)
        {
            if(_pressures[idx] < air.pressure)
            {
                double dp = air.pressure - _pressures[idx];
                _thrust = (_thrusts[idx] + dp*_dthrusts[idx]);
                break;
            }
        }
    }
    double throttle_fraction = this->_mass_rate*this->_inv_max_mass_rate;
    _thrust *= throttle_fraction;
}


void DoubleGimbalControl::move(double , double dt)
{
    Eigen::Vector3d crossP = _thrust_vector_body.cross(_commanded_thrust_vector_body);
    const double currentAngle = crossP.norm();
    double angle = std::min(currentAngle, max_slew_rate*dt);

    Eigen::AngleAxisd rotation(angle, crossP*(1.0/currentAngle));

    _thrust_vector_body = rotation*_thrust_vector_body;
}

void DoubleGimbalControl::set_commanded_thrust_vector_body(const Eigen::Vector3d& vector)
{
    if(vector.x() < min_x_component)
    {
        double factor = max_t_component/sqrt(vector[1]*vector[1] + vector[2]*vector[2]);
        _commanded_thrust_vector_body[0] = min_x_component;
        _commanded_thrust_vector_body[1] = vector[1]*factor;
        _commanded_thrust_vector_body[2] = vector[2]*factor;
    }
    else
    {
        _commanded_thrust_vector_body = vector;
    }
}
