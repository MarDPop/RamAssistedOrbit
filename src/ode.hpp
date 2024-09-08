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

    std::array<double, T::N_STATES> _dx;

    double _t;

    double _dt;

    Recording<T::N_STATES> _recording;

    bool _run_options_terminated;

    bool _dynamics_terminated;

    void euler_step(const double dt)
    {
        _f(_x, _t, _dx);
        for(auto i = 0u; i < T::N_STATES; i++)
        {
            _x[i] += _dx[i]*dt;
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

        _run_options_terminated = false;
        _dynamics_terminated = false;

        _f(options.initial_state, options.initial_time, _dx);

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

        _run_options_terminated = options.stop(_x);
        _dynamics_terminated = !_f;

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

    bool dynamics_stopped() const 
    {
        return _dynamics_terminated;
    }

    bool options_stopped() const 
    {
        return _run_options_terminated;
    }

    [[nodiscard]] const Recording<T::N_STATES>& recording() const
    {
        return _recording;
    }

};

template<typename T>
class ODE_RK4 final : public virtual ODE<T>
{
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

        this->_f(this->_x, this->_t, this->_dx);

        this->_t += _dt_half;

        for(auto i = 0u; i < T::N_STATES; i++)
        {
            this->_x[i] = _x0[i] + this->_dx[i]*_dt_half;
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
            _k3[i] = this->_dx[i] + _k2[i] + _k2[i] + _k4[i];
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

    std::array<double, T::N_STATES> _k1;
    std::array<double, T::N_STATES> _x0;
    std::array<double, T::N_STATES> _x1;
    double _t0;

    void _step(const TimestepOptions<T::N_STATES>& options) override
    {
        _x0 = this->_x;
        _t0 = this->_t;
        this->_f(this->_x, this->_t, this->_dx);
        for(auto iter = 0u; iter < NUM_ITER; iter++)
        {
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                _x1[i] = _x0[i] + this->_dx[i]*this->_dt;
            }

            this->_f(_x1, _t0 + this->_dt, _k1);

            const double dt_half = this->_dt*0.5;
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                this->_x[i] = _x0[i] + (this->_dx[i] + _k1[i])*dt_half;
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
    static constexpr unsigned MAX_ITER = 3u;

    std::array<double, T::N_STATES> _k1;
    std::array<double, T::N_STATES> _k2;
    std::array<double, T::N_STATES> _k3;
    std::array<double, T::N_STATES> _k4;
    std::array<double, T::N_STATES> _k5;
    std::array<double, T::N_STATES> _k6;

    std::array<double, 4> h;

    std::array<double, T::N_STATES> _x0;
    double _t0;

    // RKF45 coefficients
    static constexpr double A[] = {1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0};
    static constexpr double B[] = {1.0/4.0, 3.0/32.0, 9.0/32.0, 1932.0/2197.0, -7200.0/2197.0, 
        7296.0/2197.0, 439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, -8.0/27.0, 2.0, 
        -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0};
    static constexpr double C[] = {16.0/135.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0};
    static constexpr double D[] = {-1.0/360.0, 128.0/4275.0, 2197.0/75240.0, -1.0/50.0, -2.0/55.0};

    void _step(const TimestepOptions<T::N_STATES>& options) override
    {
        _x0 = this->_x;
        _t0 = this->_t;
        this->_f(this->_x, this->_t, this->_dx);
        // note that multiplying by dt for the dt terms more efficient only with low number of states
        // also vectorization improves efficiency of multiplying constants by dt instead of for each state
        for(auto iter = 0u; iter < MAX_ITER; iter++)
        {
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                _k1[i] = this->_dx[i]*this->_dt;
                this->_x[i] = _x0[i] + _k1[i]*B[0];
            }
            this->_f(this->_x, _t0 + A[0]*this->_dt, _k2);

            h[1] = this->_dt*B[2];
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                this->_x[i] = _x0[i] + _k1[i]*B[1] + _k2[i]*h[1];
            }
            this->_f(this->_x, _t0 + A[1]*this->_dt, _k3);

            h[0] = this->_dt*B[4];
            h[1] = this->_dt*B[5];
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                this->_x[i] = _x0[i] + _k1[i]*B[3] + _k2[i]*h[0] + _k3[i]*h[1];
            }
            this->_f(this->_x, _t0 + A[2]*this->_dt, _k4);

            h[0] = this->_dt*B[7];
            h[1] = this->_dt*B[8];
            h[2] = this->_dt*B[9];
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                this->_x[i] = _x0[i] + _k1[i]*B[6] + _k2[i]*h[0] + _k3[i]*h[1] + _k4[i]*h[2];
            }
            this->_f(this->_x, _t0 + A[3]*this->_dt, _k5);
 
            h[0] = this->_dt*B[11];
            h[1] = this->_dt*B[12];
            h[2] = this->_dt*B[13];
            h[3] = this->_dt*B[14];
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                this->_x[i] = _x0[i] + _k1[i]*B[10] + _k2[i]*h[0] + _k3[i]*h[1] + _k4[i]*h[2]  + _k5[i]*h[3];
            }
            this->_f(this->_x, _t0 + A[4]*this->_dt, _k6);

            h[0] = this->_dt*C[1];
            h[1] = this->_dt*C[2];
            h[2] = this->_dt*C[3];
            h[3] = this->_dt*C[4];
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                this->_x[i] = _x0[i] + C[0]*_k1[i] + _k3[i]*h[0] + _k4[i]*h[1]  + _k5[i]*h[2] + _k6[i]*h[3];
            }

            h[0] = this->_dt*D[1];
            h[1] = this->_dt*D[2];
            h[2] = this->_dt*D[3];
            h[3] = this->_dt*D[4];
            double max_err = 0.0;
            for(auto i = 0u; i < T::N_STATES; i++)
            {
                double te = D[0]*_k1[i] + _k3[i]*h[0] + _k4[i]*h[1]  + _k5[i]*h[2] + _k6[i]*h[3];
                double rel_err = fabs(te*options.inv_absolute_error[i]);
                max_err = std::max(max_err, rel_err);
            }

            const double factor = sqrt(sqrt(max_err));// std::pow(max_err, -0.2);
            this->_dt = std::clamp(this->_dt*0.9/factor, options.min_stepsize, options.max_stepsize);
            constexpr double error_threshold = 0.95;
            if(factor < error_threshold) 
            {
                break;
            }
        }
    }
};


#endif