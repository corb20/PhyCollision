#include "PhyObject.h"

class PhySphere : public PhyObject
{
public:
    Eigen::Vector3d center;
    double radius;

    PhySphere() : center(Eigen::Vector3d(0, 0, 0)), radius(1) {}
    PhySphere(Eigen::Vector3d center, double radius) : center(center), radius(radius) {}

    Eigen::Vector3d SupportPosition(const Eigen::Vector3d& direction) override;
    IntesectData CollisionWith(PhyObject* other) override;
};