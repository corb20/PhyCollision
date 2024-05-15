#include "PhySphere.h"
#include "PhyCapsule.h"
#include "PhyBox.h"
#include "PhyConvex.h"

Eigen::Vector3d PhySphere::SupportPosition(const Eigen::Vector3d& direction)
{
    return center + radius * direction.normalized();
}

IntesectData PhySphere::CollisionWith(PhyObject* other)
{
    IntesectData data;
    data.isIntersecting = false;
    data.normal = Eigen::Vector3d(0, 0, 0);
    data.penetration = 0;

    // 首先检测other是否是球
    PhySphere* sphere = dynamic_cast<PhySphere*>(other);
    if (sphere)
    {
        Eigen::Vector3d direction = sphere->center - center;
        double distance = direction.norm();
        double penetration = radius + sphere->radius - distance;
        if (penetration > 0)
        {
            data.isIntersecting = true;
            data.normal = direction.normalized();
            data.penetration = penetration;
        }
        return data;
    }

    // 检测是否是Capsule
}