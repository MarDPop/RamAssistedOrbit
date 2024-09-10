#include "doctest.h"

#include "functions.hpp"

#include <Eigen/Dense>

TEST_CASE("test the qdot rotations") {

    Eigen::Quaterniond q1(1,0,0,0);
    Eigen::Vector3d w1(0,1,0);
    Eigen::Quaterniond qdot1;
    Eigen::Quaterniond qdot2;

    functions::quaternion_orientation_rate(q1.coeffs().data(), w1.data(), qdot1.coeffs().data());
    functions::quaternion_orientation_rate(q1, w1, qdot2);

    CHECK(fabs(qdot1.w() - qdot2.w()) < 1e-9);
    CHECK(fabs(qdot1.x() - qdot2.x()) < 1e-9);
    CHECK(fabs(qdot1.y() - qdot2.y()) < 1e-9);
    CHECK(fabs(qdot1.z() - qdot2.z()) < 1e-9);

    Eigen::Quaterniond q2(Eigen::AngleAxisd(0.2, Eigen::Vector3d::UnitX()));
    Eigen::Vector3d w2(-0.1,0.1,0.2);

    functions::quaternion_orientation_rate(q1.coeffs().data(), w1.data(), qdot1.coeffs().data());
    functions::quaternion_orientation_rate(q1, w1, qdot2);

    CHECK(fabs(qdot1.w() - qdot2.w()) < 1e-9);
    CHECK(fabs(qdot1.x() - qdot2.x()) < 1e-9);
    CHECK(fabs(qdot1.y() - qdot2.y()) < 1e-9);
    CHECK(fabs(qdot1.z() - qdot2.z()) < 1e-9);
}

TEST_CASE("test ZYRotation")
{
    double zAngle = 0.2;
    double yAngle = 0.1;
    Eigen::Matrix3d Zrotation = Eigen::AngleAxisd(zAngle, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    Eigen::Matrix3d Yrotation = Eigen::AngleAxisd(yAngle, Eigen::Vector3d::UnitY()).toRotationMatrix();

    Eigen::Matrix3d ZYrotation = Zrotation*Yrotation;

    Eigen::Matrix3d myZY;
    functions::ZY_rotation(zAngle, yAngle, myZY.data());

    Eigen::Matrix3d diff = ZYrotation - myZY;

    for(auto i = 0u; i < 9; i++)
    {
        CHECK(fabs(diff.data()[i]) < 1e-6);
    }
}