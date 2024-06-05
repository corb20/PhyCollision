#include "PhySphere.h"
#include "PhyCapsule.h"
#include "PhyBox.h"
#include "PhyConvex.h"
#include "PhyAlgrithom.h"
#include "Logger.h"

//速度过快，可能会穿过，这是一个需要思考和注意的问题

Eigen::Vector3d PhySphere::SupportPosition(const Eigen::Vector3d& direction)
{
    return center + radius * direction.normalized();
}

//IntesectData normal的方向是other向外的方向，表示为other最快离开重合的方向
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
            data.normal = direction;
            //需要处理normal为0的情况，球心重合
            if(!data.normal.isApprox(Eigen::Vector3d(0, 0, 0)))
                data.normal.normalize();
            data.penetration = penetration;
        }
        return data;
    }

    // 检测是否是Capsule
    PhyCapsule* capsule = dynamic_cast<PhyCapsule*>(other);
    if (capsule)
    {
        Eigen::Vector3d closestPoint;
        double dist = capsule->segment.GetClosestDist(center,closestPoint);
        double penetration = radius + capsule->radius - dist;
        if (penetration > 0)
        {
            data.isIntersecting = true;
            data.normal = (closestPoint - center);
            //需要处理normal为0的情况，球心和胶囊轴线重合,这时的normal方向很难评，需要根据物理碰撞时的速度和方向来判断，而且重合过深，属于比较劣质的情况
            if(!data.normal.isApprox(Eigen::Vector3d(0, 0, 0))){
                data.normal.normalize();
            }
            data.penetration = penetration;
        }
        return data;
    }

    // 检测是否是Box
    PhyBox* box = dynamic_cast<PhyBox*>(other);
    if (box)
    {
        Eigen::Vector3d closestPoint;
        double dist = box->GetClosestDist(center,closestPoint);
        double penetration = radius - dist;
        if (penetration > 0)
        {
            data.isIntersecting = true;
            //如果球心也进去，要做对应的取反操作
            int isNag = dist>0?1:-1;
            data.normal = (closestPoint - center)*isNag;
            if(data.normal.isApprox(Eigen::Vector3d(0, 0, 0))){
                //沿着box该边的法向量？
                data.normal = box->GetRandomNormal(closestPoint);
            }
            data.normal.normalize();
            data.penetration = penetration;
        }
        return data;
    }

    // 检测是否是Convex
    PhyConvex* convex = dynamic_cast<PhyConvex*>(other);
    if (convex)
    {
        //使用闵可夫斯基差来计算
        data = PhyA::MKFSJCollision(this,convex);
        return data;
    }

    //other不合法，返回默认值
    Logger::logError("PhySphere CollisionWith: other is not a valid PhyObject");
    return data;
}