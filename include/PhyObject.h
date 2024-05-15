#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>

using namespace Eigen;

struct IntesectData
{
    // 是否相交
    bool isIntersecting;
    // 穿透法向
    Vector3d normal;
    // 穿透深度
    float penetration;
};

class PhyObject
{
public:
    //获取支撑点坐标
    virtual Vector3d SupportPosition(const Vector3d& direction) = 0;
    // virtual void update() = 0;
    // virtual void draw() = 0;
    // virtual void applyForce(const Vector2& force) = 0;
    // virtual void applyImpulse(const Vector2& impulse) = 0;
    // virtual void setVelocity(const Vector2& velocity) = 0;
    // virtual void setMass(float mass)

};