#ifndef ODE_H
#define ODE_H

#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <type_traits>
#include <vector>

template<unsigned NSTATES>
class Dynamics
{
public:
    static constexpr unsigned N_STATES = NSTATES;

    virtual void operator()(const std::array<double, NSTATES>& x, const double t, 
        std::array<double, NSTATES>& dx) = 0;

    [[nodiscard]] virtual operator bool() const 
    {
        return true;
    }

    [[nodiscard]] virtual unsigned num_data() const 
    {
        return 0u;
    }

    [[nodiscard]] virtual std::unique_ptr<double[]> get_data() const
    {
        return nullptr;
    }
};

template<unsigned NSTATES>
struct Recording
{
    std::vector<std::array<double, NSTATES>> states;
    std::vector<std::unique_ptr<double[]>> output;
    std::vector<double> times;
};

template<unsigned NSTATES>
struct TimestepOptions
{
    double max_stepsize = 1.0;
    double min_stepsize = 1e-3;
    std::array<double, NSTATES> inv_absolute_error = {1e-6};
};

template<unsigned NSTATES>
struct RunOptions
{
    std::array<double, NSTATES> initial_state = {0};
    double initial_time = 0.0;
    double final_time = 1.0;
    double recording_interval = 1.0;
    TimestepOptions<NSTATES> timestep;
    bool (*stop)(const std::array<double, NSTATES>&) = [](const std::array<double, NSTATES>&){ return false; };
    unsigned max_steps = 1000000;
};

template<typename T>
class ODE
{
protected:
    //static_assert(std::is_convertible<T*, Dynamics<T::N_STATES>*>::value, "Derived must inherit dynamics as public");

    T& _f;

    std::array<double, T::N_STATES> _x;

    double _t;

    double _dt;

    Recording<T::N_STATES> _recording;

    void euler_step(const double dt)
    {
        std::array<double, T::N_STATES> dx;
        _f(_x, _t, dx);
        for(auto i = 0u; i < T::N_STATES; i++)
        {
            _x[i] += dx[i]*dt;
        }
        _t += _dt;
    }

    virtual void _setup([[maybe_unused]] const RunOptions<T::N_STATES>& options) {}

    virtual void _step([[maybe_unused]] const TimestepOptions<T::N_STATES>& options)
    {
        euler_step(_dt);
    }

public:

    ODE(T& f) :
        _f(f) {}

    void run(const RunOptions<T::N_STATES>& options)
    {
        assert(options.stop != nullptr);

        _recording = Recording<T::N_STATES>();
        double time_to_record = options.initial_time;
        _x = options.initial_state;
        _t = options.initial_time;
        _dt = options.timestep.min_stepsize;

        _setup(options);

        unsigned step = 0;
        while(step++ < options.max_steps && _t < options.final_time && _f && !options.stop(_x))
        {
            // record if necessary
            if(_t >= time_to_record)
            {
                _recording.times.push_back(_t);
                _recording.output.push_back(std::move(_f.get_data()));
                _recording.states.push_back(_x);
                time_to_record += options.recording_interval;
            }

            // perform step
            _step(options.timestep);
        }

        // back track if needed
        if(_t > options.final_time)
        {
            euler_step(options.final_time - _t);
            _t = options.final_time;
        }

        // get final state and time no matter what
        _recording.times.push_back(_t);
        _recording.output.push_back(std::move(_f.get_data()));
        _recording.states.push_back(_x);
    }

    [[nodiscard]] const Recording<T::N_STATES>& recording() const
    {
        return _recording;
    }

};

template<typename T>
class ODE_RK4 final : public virtual ODE<T>
{
    std::array<double, T::N_STATES> _k1;
    std::array<double, T::N_STATES> _k2;
    std::array<double, T::N_STATES> _k3;
    std::array<double, T::N_STATES> _k4;

    std::array<double, T::N_STATES> _x0;

    double _dt_half;
    double _dt_sixth;

    void _setup(const RunOptions<T::N_STATES>& ) override
    {
        _dt_half = this->_dt*0.5;
        _dt_sixth = this->_dt*0.166666666666666666;
    }

    void _step(const TimestepOptions<T::N_STATES>&) override
    {
        _x0 = this->_x;

        this->_f(this->_x, this->_t, _k1);

        this->_t += _dt_half;

        for(auto i = 0u; i < T::N_STATES; i++)
        {
            this->_x[i] = _x0[i] + _k1[i]*_dt_half;
        }
        this->_f(this->_x, this->_t, _k2);

        for(auto i = 0u; i < T::N_STATES; i++)
        {
            this->_x[i] = _x0[i] + _k2[i]*_dt_half;
        }
        this->_f(this->_x, this->_t, _k3);

        this->_t += _dt_half;

        for(auto i = 0u; i < T::N_STATES; i++)
        {
            this->_x[i] = _x0[i] + _k3[i]*this->_dt;
        }
        this->_f(this->_x, this->_t, _k4);

        for(auto i = 0u; i < T::N_STATES; i++)
        {
            _k2[i] += _k3[i];
        }

        for(auto i = 0u; i < T::N_STATES; i++)
        {
            _k3[i] = _k1[i] + _k2[i] + _k2[i] + _k4[i];
        }

        for(auto i = 0u; i < T::N_STATES; i++)
        {
            this->_x[i] = _x0[i] + _k3[i]*_dt_sixth;
        }
    }

public:

    ODE_RK4(T& f) : 
        ODE<T>(f) {}
};

template<typename T>
class ODE_HUEN_EULER final : public virtual ODE<T>
{
    static constexpr unsigned NUM_ITER = 2u;

    std::array<double, T::N_STATES> _k0;
    std::array<double, T::N_STATES> _k1;
    std::array<double, T::N_STATES> _x0;
    std::array<double, T::N_STATES> _x1;
    double _t0;

    void _step(const TimestepOptions<T::N_STATES>& options) override
    {
        _x0 = this->_x;
        _t0 = this->_t;
        this->_f(this->_x, this->_t, _k0);
        for(auto iter = 0u; iter < NUM_ITER; iter++)
        {
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                _x1[i] = _x0[i] + _k0[i]*this->_dt;
            }

            this->_f(_x1, _t0 + this->_dt, _k1);

            const double dt_half = this->_dt*0.5;
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                this->_x[i] = _x0[i] + (_k0[i] + _k1[i])*dt_half;
            }

            // compute diff
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                _x1[i] = fabs((_x1[i] - this->_x[i])*options.inv_absolute_error[i]);
            }

            // get max diff for step size change
            double max_diff = 0;
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                max_diff = std::max(max_diff, _x1[i]);
            }

            max_diff = std::max(max_diff, 0.5); // don't grow stepsize more than 2x per step

            // Don't grow or shrink step size too much
            this->_dt = std::clamp(this->_dt/max_diff, options.min_stepsize, options.max_stepsize);
        }
        this->_t += this->_dt;
    }

public:

    ODE_HUEN_EULER(T& f) : 
        ODE<T>(f) {}
};

template<typename T>
class ODE_RK45 final : public virtual ODE<T>
{
    static constexpr unsigned MAX_ITER = 20u;

    std::array<double, T::N_STATES> _k1;
    std::array<double, T::N_STATES> _k2;
    std::array<double, T::N_STATES> _k3;
    std::array<double, T::N_STATES> _k4;
    std::array<double, T::N_STATES> _k5;

    std::array<double, T::N_STATES> _x0;
    double _t0;

    void _step(const TimestepOptions<T::N_STATES>& options) override
    {
        _x0 = this->_x;
        _t0 = this->_t;
        for(auto iter = 0u; iter < MAX_ITER; iter++)
        {

        }
    }
};


#endif