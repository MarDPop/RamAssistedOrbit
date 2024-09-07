#ifndef AERODYNAMICS_H
#define AERODYNAMICS_H

#include "atmosphere.hpp"

#include "input.hpp"

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Dense>

#include <array>

struct AeroQuantities
{
    double airspeed = 0.0;
    Eigen::Vector3d air_body_vector = Eigen::Vector3d(1,0,0);
    Eigen::Vector3d body_lift_vector = Eigen::Vector3d(0,0,1);
    double mach = 0.0;
    double dynamic_pressure = 0.0;
    double alpha_angle = 0.0; 
    double beta_angle = 0.0;
    double reynolds_number_per_meter = 0.0;

    void update(const Air& air, const Eigen::Vector3d& velocity_ecef, const Eigen::Matrix3d& orientation);
};

class Aerodynamics
{
protected:

    Eigen::Vector3d _body_force;

    Eigen::Vector3d _body_moment;

public:

    static constexpr unsigned N_COEF = 6;

    Aerodynamics() : _body_force(0,0,0), _body_moment(0,0,0) {}
    virtual ~Aerodynamics() {}

    virtual void update_forces(const AeroQuantities&) {};

    const Eigen::Vector3d& get_body_force() const
    {
        return _body_force;
    }

    const Eigen::Vector3d& get_body_moment() const
    {
        return _body_moment;
    }
};

class BallisticAerodynamics final : public virtual Aerodynamics
{
    const double _CDA;

public:

    BallisticAerodynamics(double CD, double A) : _CDA(-CD*A) {}

    void update_forces(const AeroQuantities& aero) override;

};

class Lift2DragAerodynamics final : public virtual Aerodynamics
{
    double _lift;

    double _drag;

public:

    const double _drag2lift;

    Lift2DragAerodynamics(double L2D) : _drag2lift(1.0/L2D) {}

    void set_lift(double lift)
    {
        _lift = lift;
        _drag = _drag2lift*lift;
    }

    double get_drag() const 
    {
        return _drag;
    }

    double get_lift() const
    {
        return _lift;
    }

    void update_forces(const AeroQuantities& aero) override;

};

class AerodynamicBasicCoefficients final : public virtual Aerodynamics
{

    double _elevator = 0;

public:

    struct Coef
    {
        double CD0;
        double K;
        double alpha0;
        double CL_alpha;
        double CM_alpha;
        double dCD_dEl;
        double dCM_dEl;
        double stall_angle;
        double max_elevator_deflection;
        double A_ref;
        double L_ref;
    };  

    const Coef coef_scaled;

    AerodynamicBasicCoefficients(const Coef& coef);

    double get_elevator() const
    {
        return _elevator;
    }

    void set_elevator(double elevator)
    {
        _elevator = std::clamp(elevator, -coef_scaled.max_elevator_deflection, 
            coef_scaled.max_elevator_deflection);
    }

    void update_forces(const AeroQuantities& aero) override;

};

class AerodynamicAeroGriddedTable final : public virtual Aerodynamics
{

public:

    static constexpr unsigned N_COEF = 6;

    static constexpr unsigned N_PARAMS = 3;

    static constexpr unsigned MACH = 0;

    static constexpr unsigned ALPHA = 1;

    static constexpr unsigned BETA = 2;

private:

    std::vector<std::array<double, N_PARAMS + N_COEF>> _table;

    std::array<double, N_COEF> _coef;

    std::vector<double> _mach;

    std::vector<double> _alpha;

    std::vector<double> _beta;

    unsigned _n_mach_divisions;

    unsigned _n_alpha_divisions;


    void update_coef(const AeroQuantities& aero);

public:

    AerodynamicAeroGriddedTable(std::vector<std::array<double, N_PARAMS + N_COEF>> table);

    static AerodynamicAeroGriddedTable load(std::string filename);

    void update_forces(const AeroQuantities& aero) override;

};


template<unsigned N_MACH, unsigned N_ALPHA, unsigned N_BETA>
class AerodynamicAeroTablePlane final : public virtual Aerodynamics
{

    static_assert(N_MACH > 1 && N_ALPHA > 1 && N_BETA > 1);

public:

    static constexpr unsigned N_PARAMS = 3;

    static constexpr unsigned MACH_PARAM_INDEX = 0;

    static constexpr unsigned ALPHA_PARAM_INDEX = 1;

    static constexpr unsigned BETA_PARAM_INDEX = 2;

    static constexpr unsigned idx_offsets[7] = {1, 
        N_BETA, 
        N_BETA + 1, 
        N_ALPHA, 
        N_ALPHA + 1,
        N_ALPHA + N_BETA, 
        N_ALPHA + N_BETA + 1};

    template<unsigned N>
    static constexpr std::array<double, N> get_delta(std::array<double, N> arr) 
    {
        static_assert(N > 1, "Need at least two elements");
        std::array<double, N> result;
        for(auto i = 1u; i < N; i++) 
        {
            result[i-1] = 1.0/(arr[i] - arr[i-1]);
        }
        result[N-1] = result[N-2];
        return result;
    }

    const std::vector<std::array<double, N_COEF>> _table;

    const std::array<double, N_MACH> MACH;

    const std::array<double, N_ALPHA> ALPHA;

    /**
     * With beta angles, we say it's symmetric and 0 at 0
     */
    const std::array<double, N_BETA> BETA;

    const std::array<double, N_MACH> D_MACH;

    const std::array<double, N_ALPHA> D_ALPHA;

    const std::array<double, N_BETA> D_BETA;

private:

    std::array<double, N_COEF> _coef;

    unsigned _mach_idx = 0u;

    unsigned _alpha_idx = 0u;
    
    unsigned _beta_idx = 0u;

    AerodynamicAeroTablePlane(const std::array<double, N_MACH>& mach, 
        const std::array<double, N_ALPHA>& alpha, const std::array<double, N_BETA>& ,
        const std::vector<std::array<double, N_COEF>>& table)  
        : MACH(mach), ALPHA(alpha), BETA(beta), _table(table),
        D_MACH(get_delta(mach)), D_ALPHA(get_delta(alpha)), D_BETA(get_delta(beta)) {}

    void update_coef(const AeroQuantities& aero)
    {
        // Things shouldn't change too fast between updates
        _mach_idx += _mach_idx < N_MACH - 1 & aero.mach > MACH[_mach_idx + 1];
        _mach_idx -= _mach_idx > 0 & aero.mach < MACH[_mach_idx];
        _alpha_idx += _alpha_idx < N_ALPHA - 1 & aero.alpha_angle > ALPHA[_alpha_idx + 1];
        _alpha_idx -= _alpha_idx > 0 & aero.alpha_angle < ALPHA[_alpha_idx];
        _beta_idx += _beta_idx < N_BETA - 1 & aero.beta_angle > BETA[_beta_idx + 1];
        _beta_idx -= _beta_idx > 0 & aero.beta_angle < BETA[_beta_idx];

        const unsigned idx = _mach_idx*N_BETA*N_ALPHA + _alpha_idx*N_BETA + _beta_idx;
        double* coefs[8];
        coefs[0] = _table[idx].data();
        coefs[1] = _table[idx + _idx_offsets[0]].data();
        coefs[2] = _table[idx + _idx_offsets[1]].data();
        coefs[3] = _table[idx + _idx_offsets[2]].data();
        coefs[4] = _table[idx + _idx_offsets[3]].data();
        coefs[5] = _table[idx + _idx_offsets[4]].data();
        coefs[6] = _table[idx + _idx_offsets[5]].data();
        coefs[7] = _table[idx + _idx_offsets[6]].data();

        const double f_mach = (aero.mach - MACH[_mach_idx])*D_MACH[_mach_idx];
        const double f_alpha = (aero.alpha_angle - ALPHA[_alpha_idx])*D_ALPHA[_alpha_idx];
        const double f_beta = (aero.beta_angle - BETA[_beta_idx])*D_BETA[_beta_idx];

        const double f1_mach = 1.0 - f_mach;
        const double f1_alpha = 1.0 - f_alpha;
        const double f1_beta = 1.0 - f_beta;

        //common factors
        const double fm1a1 = f1_mach*f1_alpha;
        const double fm1a = f1_mach*f_alpha;
        const double fma1 = f_mach*f1_alpha;
        const double fma = f_mach*f_alpha;

        const double factors[8] = 
        {
            fm1a1*f1_beta, fm1a1*f_beta, fm1a*f1_beta, fm1a*f_beta,
            fma1*f1_beta, fma1*f_beta, fma*f1_beta, fma*f_beta
        };  

        for(auto i = 0u; i < N_COEF; i++) 
        {
            _coef[i] = factors[0]*coefs[0][i] + factors[1]*coefs[1][i] + factors[2]*coefs[2][i] 
                + factors[3]*coefs[3][i] + factors[4]*coefs[4][i] + factors[5]*coefs[5][i] + factors[6]*coefs[6][i]
                + factors[7]*coefs[7][i];
        }
    }

public:     

    static AerodynamicAeroTablePlane load(std::string filename)
    {
        std::vector<std::array<double, N_COEF>> table;
        std::vector<double> mach, alpha, beta;

        std::ifstream input(filename);
        if (!input.is_open())
        {
            std::cerr << "Error opening file: " << filename << std::endl;
            return AerodynamicAeroTablePlane({}, {}, {}, table);
        }
        std::string line;
        while (std::getline(input, line))
        {
            auto tokens = input::split(line);

            if(tokens.size() != N_PARAMS + N_COEF) 
            {
                std::cerr << "Invalid line format: " << line << std::endl;
                throw std::runtime_error("Invalid line format");
            }
            std::array<double, N_COEF> row;
            mach.push_back(std::stod(tokens[0]));
            alpha.push_back(std::stod(tokens[1]));
            beta.push_back(std::stod(tokens[2]));
            for (size_t i = 0; i < N_COEF; ++i)
            {
                row[i] = std::stod(tokens[i + 3]);
            }
            table.push_back(row);
        }
        input.close();

        if(table.size() != N_MACH*N_ALPHA*N_BETA)
        {
            std::cerr << "Invalid table size: " << table.size() << " expected " << N_MACH*N_ALPHA*N_BETA << std::endl;
            throw std::runtime_error("Invalid table size");
        }

        // verify unique mach, alpha, beta values
        std::array<double, N_MACH> unique_mach;
        std::array<double, N_ALPHA> unique_alpha;
        std::array<double, N_BETA> unique_beta;

        // verify unique values and sort
        auto nUniqueMach = 0u;
        auto nUniqueAlpha = 0u;
        auto nUniqueBeta = 0u;
        for(auto i = 0u; i < table.size(); i++)
        {
            auto idxFound = 0u;
            for(idxFound = 0u; idxFound < nUniqueMach; idxFound++) 
            {
                if(mach[i] == unique_mach[idxFound]) break;
            }
            if(idxFound == nUniqueMach) 
            {
                unique_mach[nUniqueMach] = mach[i];
                nUniqueMach++;
            }
            for(idxFound = 0u; idxFound < nUniqueAlpha; idxFound++)
            {
                if(alpha[i] == unique_alpha[idxFound]) break;
            }
            if(idxFound == nUniqueAlpha)
            {
                unique_alpha[nUniqueAlpha] = alpha[i];
                nUniqueAlpha++;
            }
            for(idxFound = 0u; idxFound < nUniqueBeta; idxFound++)
            {
                if(beta[i] == unique_beta[idxFound]) break;
            }
            if(idxFound == nUniqueBeta)
            {
                unique_beta[nUniqueBeta] = beta[i];
                nUniqueBeta++;
            }
            if(nUniqueAlpha > N_ALPHA || nUniqueBeta > N_BETA || nUniqueMach > N_MACH)
            {
                std::cerr << "Duplicate or out of order values found" << std::endl;
                throw std::runtime_error("Duplicate or out of order values found");
            }
        }

        if(nUniqueAlpha < N_ALPHA || nUniqueBeta < N_BETA || nUniqueMach < N_MACH)
        {
            std::cerr << "Duplicate or out of order values found" << std::endl;
            throw std::runtime_error("Duplicate or out of order values found");
        }

        for(auto i = 0u; i < N_MACH; i++)
        {
            for(auto j = 0u; j < N_ALPHA; j++)
            {
                for(auto k = 0u; k < N_BETA; k++)
                {
                    auto idx = i*N_ALPHA*N_BETA + j*N_BETA + k;
                    if(beta[idx] != unique_beta[k] || alpha[idx] != unique_alpha[j] || mach[idx] != unique_mach[i])
                    {
                        std::cerr << "Out of order found: (" << mach[i] << ", " << alpha[j] << ", " << beta[k] << ")" << std::endl;
                        throw std::runtime_error("Out of order found");
                    }
                }
            }
        }

        AerodynamicAeroTablePlane(unique_mach, unique_alpha, unique_beta, table);
    }

    void update_forces(const AeroQuantities& aero) override
    {
        update_coef(aero);

        _body_force = Eigen::Vector3d(aero.dynamic_pressure*_coef[0], 
            aero.dynamic_pressure*_coef[1], aero.dynamic_pressure*_coef[2]);

        _body_moment = Eigen::Vector3d(aero.dynamic_pressure*_coef[3], 
            aero.dynamic_pressure*_coef[4], aero.dynamic_pressure*_coef[5]);
    }

};

#endif