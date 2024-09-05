#ifndef GNC_HPP
#define GNC_HPP

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Dense>

#include <memory>

class VehicleBase;

class Navigation
{

public:

    virtual void update([[maybe_unused]] const VehicleBase& vehicle, [[maybe_unused]] double time) {}

};

class Guidance
{

public:

    const Navigation& navigation;

    Guidance(const Navigation& navigation_) : navigation(navigation_) {} 

    virtual void update([[maybe_unused]] double time) {}

};


class Control
{

public:

    const Guidance& guidance;

    Control(const Guidance& guidance_) : guidance(guidance_) {} 

    virtual void update([[maybe_unused]] double time) {}
};

class VehicleGNC
{

public:

    virtual void update([[maybe_unused]] const VehicleBase& vehicle, [[maybe_unused]] double time) {}

};

class GNC final : public virtual VehicleGNC
{
public:

    const std::unique_ptr<Guidance> guidance;

    const std::unique_ptr<Navigation> navigation;

    const std::unique_ptr<Control> control;

    void update(const VehicleBase& vehicle, double time) override
    {
        navigation->update(vehicle, time);
        guidance->update(time);
        control->update(time);
    }

};

#endif