#include "constants.hpp"
#include "propulsion.hpp"

#ifdef DEBUG
    #include <iostream>
#endif

RamjetFixedInlet::RamjetFixedInlet(double throat_area, double exit_area,
        double nozzle_efficiency, double nozzle_ram_ratio, double dry_mass,
        double combustion_temperature, double fuel_heating_value, double max_fuel_air_ratio, 
        double gamma_combustion, double mw_combustion) :
    _throat_area(throat_area), _exit_area(exit_area),
    _nozzle_efficiency(nozzle_efficiency), _nozzle_ram_ratio(nozzle_ram_ratio),
    _combustion_temperature(combustion_temperature), _fuel_heating_value(fuel_heating_value), 
    _max_fuel_air_ratio(max_fuel_air_ratio), _gamma_combustion(gamma_combustion), 
    _cp_4(Air::GAS_CONSTANT*gamma_combustion/(gamma_combustion - 1.0)/mw_combustion),
    _specific_enthalpy_4(combustion_temperature*_cp_4),
    _exit_mach(Air::supersonic_mach_area_ratio(throat_area/exit_area, gamma_combustion)),
    _t_exit_ratio(1.0/Air::isentropic_temperature_ratio(_exit_mach, gamma_combustion)),
    _p_exit_ratio(pow(_t_exit_ratio, gamma_combustion/(gamma_combustion - 1.0))),
    _exit_const(Air::GAS_CONSTANT*gamma_combustion/mw_combustion)
{
    _dry_mass = dry_mass;
}

void RamjetFixedInlet::update_thrust(const Air& air, const AeroQuantities& aero) 
{
    constexpr double MIN_MACH = 1.1;
    if(aero.mach < MIN_MACH) 
    {
        _thrust = 0;
        _mass_rate = 0;
        return;
    }

    const auto g2 = (air.gamma - 1.0);
    const auto g4 = air.gamma/g2;
    const auto beta = 1.0 + 0.5*g2*aero.mach*aero.mach;

    const auto t_total_0 = air.temperature*beta;
    const auto p_total_0 = air.pressure*pow(beta, g4);
    const auto p_total_loss_diffusor = ramRatio(aero.mach);
    const auto p_total_2 = p_total_0*p_total_loss_diffusor;

    const auto air_mass_rate = _throat_area*Air::choked_flow(p_total_2, t_total_0);

    const auto specific_enthalpy_2 = air.specific_gas_constant*g4*t_total_0;
    
    const auto fuel_mass_rate_desired = air_mass_rate*(_specific_enthalpy_4 - specific_enthalpy_2)
        /(_fuel_heating_value - _specific_enthalpy_4);

    _mass_rate = std::min(_max_fuel_air_ratio*air_mass_rate, fuel_mass_rate_desired);

    const auto mixed_mass_rate = _mass_rate + air_mass_rate;

    const auto t_total_4 = (_mass_rate*_fuel_heating_value + air_mass_rate*specific_enthalpy_2)
        /(_cp_4*mixed_mass_rate);

    constexpr double t_total_nozzle_loss = 0.98;
    const auto t_total_6 = t_total_nozzle_loss*t_total_4;
    const auto p_total_6 = _nozzle_ram_ratio*p_total_2;
    const auto t_6 = t_total_6*_t_exit_ratio;
    const auto p_6 = p_total_6*_p_exit_ratio;
    
    const auto v_exit = _exit_mach*sqrt(_exit_const*t_6);

    _thrust = std::max(mixed_mass_rate*v_exit*_nozzle_efficiency + _exit_area*(p_6 - air.pressure),0.0);
}

RamjetFixedInlet RamjetFixedInlet::create(const double mach, const double altitude, const double thrust, 
        const double thrust2weight, const double nozzle_ram_ratio)
{
    AtmosphereLinearTable atm = AtmosphereLinearTable::create(AtmosphereLinearTable::STD_ATMOSPHERES::US_1976, 100);
    Air air;
    atm.set_air(altitude, air);
    AeroQuantities aero;
    const double cruise_speed = mach/air.inv_sound_speed;
    Eigen::Vector3d velocity(cruise_speed, 0, 0);
    Eigen::Matrix3d CS = Eigen::Matrix3d::Identity();
    aero.update(air, velocity, CS);

    const double p0 = air.pressure*Air::isentropic_pressure_ratio(aero.mach);

    const double t0 = air.temperature*Air::isentropic_temperature_ratio(aero.mach);

    const double pratio = RamjetFixedInlet::ramRatio(aero.mach);

    const double p02 = p0*pratio;

    const double mass_flux = p02/Air::choked_flow(p02, t0);

    const double p04 = p02*nozzle_ram_ratio;

    constexpr double gamma_products = 1.3;
    constexpr double mw_products = 0.042;

    const double mach_exit = Air::mach_from_isentropic_pressure_ratio(p04/air.pressure, gamma_products);

    const double beta_exit = Air::isentropic_temperature_ratio(mach_exit, gamma_products);
    const double t04 = 2800;
    const double t_exit = t04/beta_exit;
    const double v_exit = mach_exit*sqrt(gamma_products*Air::GAS_CONSTANT/mw_products*t_exit);
    const double nozzle_efficiency = 0.95;

    const double mass_flow = thrust/(nozzle_efficiency*v_exit);

    const double throat_area = mass_flow/mass_flux;

    const double exit_area = throat_area/Air::isentropic_area_ratio(mach_exit, gamma_products);

    return RamjetFixedInlet(throat_area, exit_area, nozzle_efficiency, 0.98, thrust/thrust2weight/GForce::G);
}

RamjetVariableInlet::RamjetVariableInlet( double max_mass_rate, double throat_area, double exit_area, double max_intake_area,
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
    _max_exit_area(exit_area),
    _max_exit_mach(Air::supersonic_mach_area_ratio(throat_area/exit_area, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _max_mach_exit_pressure_ratio(nozzle_pressure_ratio/Air::isentropic_pressure_ratio(_max_exit_mach, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _max_mach_exit_temperature_ratio(adiabatic_efficiency/Air::isentropic_temperature_ratio(_max_exit_mach, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _min_exit_area(exit_area*0.5),
    _min_exit_mach(Air::supersonic_mach_area_ratio(throat_area/_min_exit_area, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _min_mach_exit_pressure_ratio(nozzle_pressure_ratio/Air::isentropic_pressure_ratio(_min_exit_mach, GAMMA_COMBUSTION_PRODUCTS_KEROSENE)),
    _min_mach_exit_temperature_ratio(adiabatic_efficiency/Air::isentropic_temperature_ratio(_min_exit_mach, GAMMA_COMBUSTION_PRODUCTS_KEROSENE))
    {}

RamjetVariableInlet RamjetVariableInlet::create(double thrust2weight, double altitude, double exit_mach, double cruise_mach,
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
        RamjetVariableInlet ramjet(1e10, throat_area, exit_area, throat_area*5e3);
        ramjet.update_thrust(air, aero);
        double current_thrust = ramjet.get_thrust();

        double dA = throat_area*AREA_FRACTION;
        double throat_area1 = throat_area + dA;
        exit_area = exit_area_ratio*throat_area1;
        RamjetVariableInlet ramjet2(1e10, throat_area1, exit_area, throat_area*5e3);
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
    RamjetVariableInlet ramjet(mass_rate*mass_rate_margin, throat_area, exit_area, intakeAreaIdeal*1.2);
    ramjet._dry_mass = desired_thrust/(GForce::G*thrust2weight);
    ramjet.update_thrust(air, aero);
    return ramjet;
}

bool RamjetVariableInlet::can_turn_on(const Air& air, const AeroQuantities& aero) const 
{
    double intake_area = _throat_area/Air::isentropic_area_ratio(aero.mach);
    const auto mdot_air = aero.airspeed*air.density*std::min(intake_area,_max_intake_area)*(-aero.air_body_vector.x());
    return mdot_air < _max_air_ingest;
}

void RamjetVariableInlet::update_thrust(const Air& air, const AeroQuantities& aero) 
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
