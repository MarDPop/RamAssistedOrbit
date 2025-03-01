#ifndef PROPULSION_H
#define PROPULSION_H

#include "aerodynamics.hpp"
#include "atmosphere.hpp"
#include "constants.hpp"
#include "functions.hpp"

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Dense>

#include <cmath>
#include <exception>
#include <memory>

/**
 * Base class for all thrusters, provides required data fields to update
 * If not inherited assumes that the thruster provides fixed values and can be
 * restarted freely
 */
class Thruster
{
protected:

    /**
     * Mass of the Thruster (kg)
     */
    double _dry_mass = 0.0;

    /**
     * Current mass rate of fuel to thruster (kg/s)
     */
    double _mass_rate = 0.0;

    /**
     * Current thrust (N)
     */
    double _thrust = 0.0;

    /**
     * Flag if the thruster is running (if false, )
     */
    bool _running = true;

public:

    /**
     * default constructor doesn't initialize anything
     */
    Thruster(){}

    /**
     * For really simple thrusters
     */
    Thruster(double thrust, double isp) : _mass_rate(thrust/(isp*GForce::G)), _thrust(thrust)  {}

    virtual ~Thruster(){}

    /**
     * Gets the mass of the thruster
     */
    [[nodiscard]] double get_dry_mass() const
    {
        return _dry_mass;
    }

    /**
     * Gets the current mass rate (kg/s)
     */
    [[nodiscard]] double get_mass_rate() const
    {
        return _running*_mass_rate;
    }

    /**
     * Gets the current Thrust (N)
     */
    [[nodiscard]] double get_thrust() const
    {
        return _running*_thrust;
    }

    /**
     * Checks if the thruster is currently running
     */
    [[nodiscard]] bool is_running() const
    {
        return _running;
    }

    /**
     * Starts the thruster
     */
    virtual void start()
    {
        this->_running = true;
    }

    /**
     * Sets the throttle position for the thruster.
     *
     * @param throttle The new throttle position, typically a value between 0.0 (off) and 1.0 (full throttle).
     *                  This value should be clamped to the valid range for the specific thruster.
     *
     * @note This function is intended to be overridden for real thrusters to provide functionality.
     *       For virtual thrusters, it does nothing.
     */
    virtual void set_throttle(double throttle) 
    {
        // do nothing, implement for real thrusters
    }

    /**
     * Stops the thruster
     * Note this might need to be overriden for solid rockets
     */
    virtual void stop()
    {
        this->_running = false;
    }

    /**
     * 
     */
    virtual void update_thrust([[maybe_unused]] const Air& air, [[maybe_unused]] const AeroQuantities& aero) {}
};

class ThrusterTabulated : public virtual Thruster 
{
    const std::vector<double> _pressures;

    const std::vector<double> _dpressures;

    const std::vector<double> _temperatures;

    const std::vector<double> _dtemperatures;   

    const std::vector<double> _machs;

    const std::vector<double> _dmachs;

    const std::vector<double> _throttles;

    const std::vector<double> _dthrottles;

    const std::vector<double> _beta_angle;

    const std::vector<double> _dbeta_angle;

    const std::vector<double> _alpha_angle;

    const std::vector<double> _dalpha_angle;

    const std::vector<double> _thrust_values;

    unsigned layer_offset[6];

    const double _min_mass_rate;

    const double _max_mass_rate;

    double _throttle = 1.0;

public: 

    ThrusterTabulated(const std::vector<double>& pressures,
        const std::vector<double>& temperatures, const std::vector<double>& machs,
        const std::vector<double>& throttles, const std::vector<double>& beta_angle,
        const std::vector<double>& alpha_angle, const std::vector<double>& thrust_values,
        double min_mass_rate, double max_mass_rate);

    void set_throttle(double throttle) override
    {
        this->_throttle = std::clamp(throttle, 0.0, 1.0);
        this->_mass_rate = _min_mass_rate + (_max_mass_rate - _min_mass_rate)*_throttle;
    }
    
    void update_thrust(const Air& air, const AeroQuantities& aero) override;
};

class Ramjet : public virtual Thruster
{

public:

    virtual ~Ramjet() {}

    virtual bool can_operate(double air_pressure, double mach) const = 0;
};

class RamjetSimpleModel final : public virtual Ramjet
{
    const double _throat_area;

    const double _min_exit_area;

    const double _max_exit_area;

    const double _diffusor_efficiency;

    const double _burner_pressure_ratio;

    const double _burner_efficiency;

    const double _nozzle_efficiency; // thrust losses due to friction and off axis velocity

    const double _combustion_temperature;

    const double _fuel_heating_value;

    const double _max_fuel_air_ratio;

    const double _max_mass_rate;

    const double _gamma_combustion;

    const double _cp_4;

    const double _g1_4;

    const double _g2_4;

    const double _specific_enthalpy_4;

    const double _min_exit_mach;

    const double _max_exit_mach;

    const double _max_exit_area_temperature_ratio;

    const double _min_exit_area_temperature_ratio;    

    const double _max_exit_area_pressure_ratio;

    const double _min_exit_area_pressure_ratio;

    const double _exit_const;

    const double _nH;

    const double _cp_ratio4;

    const double _h04_max;

    double _A_exit = 1.0;

    double _v_exit;

    double _M_exit;

    double _p_exit;

    double _T_exit;

public:

    RamjetSimpleModel(double throat_area, 
        double min_exit_area, 
        double max_exit_area,
        double diffusor_efficiency, // typical range 0.90 - 0.95
        double burner_pressure_ratio, // typical range 0.95 - 0.99
        double burner_efficiency, // typical range 0.91 - 0.98
        double nozzle_efficiency, // typical range 0.95 - 0.99
        double dry_mass, 
        double max_mass_rate,
        double combustion_temperature = 2800, 
        double fuel_heating_value = 45e6, 
        double max_fuel_air_ratio = 0.066, 
        double gamma_combustion = 1.3, 
        double mw_combustion = 0.040);

    RamjetSimpleModel(const RamjetSimpleModel& copy) = default;
    RamjetSimpleModel(RamjetSimpleModel&& move) noexcept = default;
    RamjetSimpleModel& operator=(const RamjetSimpleModel& copy) = default;
    RamjetSimpleModel& operator=(RamjetSimpleModel&& move) noexcept = default;

    static double ram_recovery_factor(double mach)
    {
        return std::max(1.0 - 0.01*mach*mach, 0.5);
    }

    bool can_operate(double air_pressure, double mach) const override
    {
        constexpr double MIN_MACH = 1.0;
        return mach > MIN_MACH;
    }

    void update_thrust(const Air& air, const AeroQuantities& aero) override;
};

class RamjetReal final : public virtual Ramjet
{
    /**
     * Heating value of the fuel in J/kg
     */
    const double _heating_value_fuel;

    /**
     * Optimum fuel to air ratio
     */
    const double _fuel_air_ratio;

    /**
     * Maximum fuel rate the engine can provide (kg/s)
     */
    const double _max_mass_rate;

    /**
     * Maximum air rate before engine can turn on (kg/s)
     */
    const double _max_air_ingest;

    /**
     * Effiency of the combustor [accounts for heat and pressure losses] (0.9 - 0.97 typical)
     */
    const double _combustion_efficiency;

    /**
     * 
     */
    const double _combustor_pressure_ratio;

    /**
     * Efficiency of the nozzle [axial velocity losses] (0.92 - 0.99 typical)
     */
    const double _nozzle_efficiency;

    /**
     * 
     */
    const double _nozzle_pressure_ratio;

    /**
     * 
     */
    const double _adiabatic_efficiency;

    /**
     * Area of the throat (m^2)
     */
    const double _throat_area;

    /**
     * Maximum size of the intake (m^2)
     */
    const double _max_intake_area;

    /**
     * mach number at which the inlet no longer can isentropically diffuse
     */
    const double _critical_mach;

     /**
     * Area of the nozzle exit (m^2)
     */
    const double _max_exit_area;

    /**
     * Optimal mach number of exit 
     */
    const double _max_exit_mach;

    /**
     * Pressure of the nozzle exit from total pressure at combustor
     */
    const double _max_mach_exit_pressure_ratio;

    /**
     * 
     */
    const double _max_mach_exit_temperature_ratio;

    /**
     * @brief 
     * 
     */
    const double _nominal_exit_area;

     /**
     * Area of the nozzle exit (m^2)
     */
    const double _min_exit_area;

    /**
     * Optimal mach number of exit 
     */
    const double _min_exit_mach;

    /**
     * Pressure of the nozzle exit from total pressure at combustor
     */
    const double _min_mach_exit_pressure_ratio;

    /**
     * 
     */
    const double _min_mach_exit_temperature_ratio;

public:

    enum STAGE : int {
        AMBIENT = 0,
        DIFFUSOR = 1,
        COMBUSTOR_ENTRANCE = 2,
        COMBUSTOR = 3,
        COMBUSTOR_EXIT = 4,
        NOZZLE = 5,
        NOZZLE_EXIT = 6,
        NUM_STAGES 
    };

    static constexpr double KEROSENE_HEATING_VALUE = 47e6; // J / kg
    static constexpr double KEROSENE_FUEL_AIR_RATIO = 1.0/15.0;

    static constexpr double R_GAS_COMBUSTION_PRODUCTS_KEROSENE = 269.1;
    static constexpr double GAMMA_COMBUSTION_PRODUCTS_KEROSENE = 1.3;
    static constexpr double INV_EXP_COMBUSTION_PRODUCTS_KEROSENE = 0.230769230769231;

    /**
     * @param max_mass_rate kg/s
     * @param throat_area m^2
     * @param exit_are m^2
     * @param heating_value_fuel J/kg
     */
    RamjetReal( double max_mass_rate, double throat_area, double nominal_exit_area, 
        double min_exit_area, double max_exit_area, double max_intake_area,
        double heating_value_fuel = KEROSENE_HEATING_VALUE, double fuel_air_ratio = KEROSENE_FUEL_AIR_RATIO, 
        double combustion_efficiency = 1.0, double combustor_pressure_ratio = 1.0, double nozzle_efficiency = 1.0,
        double nozzle_pressure_ratio = 1.0, double adiabatic_efficiency = 1.0);

    RamjetReal(const RamjetReal& copy) = delete; // try not to do this
    RamjetReal(RamjetReal&& copy) = default;
    RamjetReal& operator=(const RamjetReal& copy) = default; // this is ok
    RamjetReal& operator=(RamjetReal&& copy) = default;

    /**
     * Creates a ramjet
     */
    static RamjetReal create(double thrust2weight, double altitude, double exit_mach, double cruiseSpeed,
        double lift2drag, double mass, double thrust_margin = 1.1, double mass_rate_margin = 2.0);

    /**
     * 
     */
    [[nodiscard]] bool can_turn_on(const Air& air, const AeroQuantities& aero) const;

    /**
     * 
     * 
     */
    void update_thrust(const Air& air, const AeroQuantities& aero) override;

    bool can_operate(double air_pressure, double mach) const override
    {
        return mach > _min_exit_mach;
    }
};

/**
 * Class rocket is a simple thruster that has a pressure table defined performance
 * 
 */
class Rocket final : public virtual Thruster
{
    /**
     * pressures for table (Pa)
     */
    std::vector<double> _pressures;

    /**
     * Thrusts for table (N)
     */
    std::vector<double> _thrusts;

    /**
     * Linear interpolate factors for table (N/Pa)
     */
    std::vector<double> _dthrusts;

    /**
     * mass rate for full throttle (kg/s)
     */
    const double _max_mass_rate;

    /**
     * mass rate for min throttle (kg/s)
     */
    const double _min_mass_rate;

    /**
     * change in mass for throttle range (kg/s)
     */
    const double _dmass_rate_throttle;

    /**
     * To help with throttle factor
     */
    const double _inv_max_mass_rate;

    /**
     * Current throttle setting (remember 0 does not equal off)
     */
    double _throttle = 1.0;

    /**
     * Max Number of times this thruster can be restarted (typically 1)
     */
    const int _max_number_restarts;

    /**
     * Number of times this thruster has been restarted
     */
    int _restarted = 0;

public:

    Rocket(const std::vector<double>& pressures, const std::vector<double>& thrusts, 
        double max_mass_rate, double min_mass_rate, int max_number_restarts);

    /**
     * Sets throttle for the rocket (if throttlable)
     * @param throttle (range: 0 - 1.0)
     */
    void set_throttle(double throttle)
    {
        this->_throttle = std::clamp(throttle, 0.0, 1.0);
        this->_mass_rate = _min_mass_rate + _dmass_rate_throttle*this->_throttle;
    }

    /**
     * Resets the restarts back to 0
     */
    void reset_restarts()
    {
        this->_restarted = 0;
    }

    /**
     * Starts the engine if possible (won't start if max number of restarts reached)
     */
    void start() final 
    {
        this->_restarted++;
        this->_running = this->_restarted <= this->_max_number_restarts;
    }

    /**
     * Updates thrust and mass rate due to altitude
     */
    void update_thrust(const Air& air, const AeroQuantities& aero) override;
};

class ThrusterControl
{
protected:

    Eigen::Vector3d _thrust_vector_body{};

    Eigen::Vector3d _commanded_thrust_vector_body{};

    double _old_time{};

public:

    virtual ~ThrusterControl(){}

    const Eigen::Vector3d& get_thrust_vector_body() const
    {
        return _thrust_vector_body;
    }

    const Eigen::Vector3d& get_commanded_thrust_vector_body() const
    {
        return _commanded_thrust_vector_body;
    }

    virtual void set_commanded_thrust_vector_body(const Eigen::Vector3d& vector)
    {
        _commanded_thrust_vector_body = vector;
    }

    virtual void move([[maybe_unused]] double time, [[maybe_unused]] double dt) 
    {
        _thrust_vector_body = _commanded_thrust_vector_body;
    }

    void update(double time) 
    {
        move(time, time - _old_time);
        _old_time = time;
    }

};

/**
 * A controllable thrust vectored rocket
 */
class DoubleGimbalControl final : public virtual ThrusterControl
{

public: 

    const double min_x_component;

    const double max_t_component;

    const double max_slew_rate;

    DoubleGimbalControl(double max_angle, double max_slew_rate_) :
        min_x_component(cos(max_angle)), max_t_component(sin(max_angle)),
        max_slew_rate(max_slew_rate_)  {}

    DoubleGimbalControl(const DoubleGimbalControl& other) : min_x_component(other.min_x_component),
        max_t_component(other.max_t_component), max_slew_rate(other.max_slew_rate) {}

    void move(double time, double dt);

    void set_commanded_thrust_vector_body(const Eigen::Vector3d& vector) override;

};

template<class T, class C>
class PropulsionSystem_t
{
    Eigen::Vector3d _body_force;

    Eigen::Vector3d _body_moment;

    T _thruster;

    C _control;

public:

    const Eigen::Vector3d position;

    PropulsionSystem_t(T thruster, C control, const Eigen::Vector3d& position_) : _thruster(std::move(thruster)),
        _control(std::move(control)), position(position_) {}

    const Eigen::Vector3d& get_body_force() const
    {
        return _body_force;
    }

    const Eigen::Vector3d& get_body_moment() const
    {
        return _body_moment;
    }

    void update(const Air& air, const AeroQuantities& aero, const Eigen::Vector3d& COM, double time) 
    {
        _thruster.update_thrust(air, aero);
        _control.update(time);
        _body_force = _thruster.get_thrust()*_control.get_thrust_vector_body();
        _body_moment = (position - COM).cross(_body_force);
    }
};

class PropulsionSystem final
{
    Eigen::Vector3d _body_force = Eigen::Vector3d(0,0,0);

    Eigen::Vector3d _body_moment = Eigen::Vector3d(0,0,0);

    Eigen::Vector3d _position = Eigen::Vector3d(0,0,0);

public:

    std::unique_ptr<Thruster> const thruster;

    std::unique_ptr<ThrusterControl> const control;

    PropulsionSystem(std::unique_ptr<Thruster> thruster_, std::unique_ptr<ThrusterControl> control_) :
        thruster(std::move(thruster_)), control(std::move(control_)) 
    {
        if(!thruster || !control)
        {
            throw "";
        }
    }

    const Eigen::Vector3d& get_body_force() const
    {
        return _body_force;
    }

    const Eigen::Vector3d& get_body_moment() const
    {
        return _body_moment;
    }

    void set_position(const Eigen::Vector3d& position)
    {
        _position = position;
    }

    void update_forces(const Air& air, const AeroQuantities& aero, const Eigen::Vector3d& COM, double time) 
    {
        thruster->update_thrust(air, aero);
        control->update(time);
        _body_force = thruster->get_thrust()*control->get_thrust_vector_body();
        _body_moment = (_position - COM).cross(_body_force);
    }
};

#endif