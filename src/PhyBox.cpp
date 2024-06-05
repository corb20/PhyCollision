#include "PhyBox.h"
#include "PhyBaseStruct.h"
#include "Logger.h"

std::vector<Eigen::Vector3d> PhyBox::GetLocalVertices()
{
    std::vector<Eigen::Vector3d> vertices;
    vertices.push_back(Eigen::Vector3d(minP[0], minP[1], minP[2]));
    vertices.push_back(Eigen::Vector3d(maxP[0], minP[1], minP[2]));
    vertices.push_back(Eigen::Vector3d(maxP[0], maxP[1], minP[2]));
    vertices.push_back(Eigen::Vector3d(minP[0], maxP[1], minP[2]));
    vertices.push_back(Eigen::Vector3d(minP[0], minP[1], maxP[2]));
    vertices.push_back(Eigen::Vector3d(maxP[0], minP[1], maxP[2]));
    vertices.push_back(Eigen::Vector3d(maxP[0], maxP[1], maxP[2]));
    vertices.push_back(Eigen::Vector3d(minP[0], maxP[1], maxP[2]));
    return vertices;
}
std::vector<Eigen::Vector3d> PhyBox::GetWorldVertices()
{
    std::vector<Eigen::Vector3d> vertices = GetLocalVertices();
    for (auto& vertex : vertices)
    {
        vertex = coord.ConversionToWorld(vertex);
    }
    return vertices;
}

double PhyBox::GetClosestDist(const Eigen::Vector3d& point, Eigen::Vector3d& closestPoint){
    //将点转换到局部坐标系
    Eigen::Vector3d localPoint = coord.ConversionFromWorld(point);
    //将所有的坐标中心变成center
    Eigen::Vector3d boxCenter = (maxP + minP)/2;
    Eigen::Vector3d a = maxP - boxCenter;
    Eigen::Vector3d p = localPoint - boxCenter;
    Eigen::Vector3d mask = Vector3d(1,1,1);
    //将p全部变为正值，并且计算变换mask
    for(int i = 0;i<3;i++){
        if(p[i] < 0){
            p[i] = -p[i];
            mask[i] = -1;
        }
    }
    Eigen::Vector3d ap = p - a;
    Eigen::Vector3d apv = Vector3d(0,0,0);
    if(ap[0]>0 || ap[1]>0 || ap[2]>0)
    {
        apv = Vector3d(fmax(ap[0],0),fmax(ap[1],0),fmax(ap[2],0));
    }
    else
    {
        //全小于0，说明点在长方体内部
        int minIndex = 0;
        double minDist = ap[0];
        for(int i = 1;i<3;i++){
            //越接近0，说明越靠近长方体表面
            if(ap[i] > minDist){
                minDist = ap[i];
                minIndex = i;
            }
        }
        apv[minIndex] = minDist;
    }
    closestPoint = boxCenter + (p - apv).cwiseProduct(mask);
    //转到世界坐标系
    closestPoint = coord.ConversionToWorld(closestPoint);
    return apv.norm();
}

Eigen::Vector3d PhyBox::GetRandomNormal(Eigen::Vector3d p){
    //如果p在某个面上，并且不在边上，则返回该面的法向量
    Eigen::Vector3d localP = coord.ConversionFromWorld(p) - minP;
    Eigen::Vector3d wlh = maxP - minP;
    //超范围了或者在边界，说明是不合格的点
    if(localP.x()<0 || localP.x()>wlh.x() || localP.y()<0 || localP.y()>wlh.y() || localP.z()<0 || localP.z()>wlh.z()){
        Logger::logError("PhyBox::GetRandomNormal: The point is out of range or on the edge");
        return Vector3d(0,0,0);
    }
    //在某个面上
    if(IsEqual(localP.x(),0) && !IsEqual(localP.y(),0) && !IsEqual(localP.z(),0) && !IsEqual(localP.y(),wlh.y()) && !IsEqual(localP.z(),wlh.z())){
        return -coord.i;
    }
    if(IsEqual(localP.x(),wlh.x()) && !IsEqual(localP.y(),0) && !IsEqual(localP.z(),0) && !IsEqual(localP.y(),wlh.y()) && !IsEqual(localP.z(),wlh.z())){
        return coord.i;
    }
    if(IsEqual(localP.y(),0) && !IsEqual(localP.x(),0) && !IsEqual(localP.z(),0) && !IsEqual(localP.x(),wlh.x()) && !IsEqual(localP.z(),wlh.z())){
        return -coord.j;
    }
    if(IsEqual(localP.y(),wlh.y()) && !IsEqual(localP.x(),0) && !IsEqual(localP.z(),0) && !IsEqual(localP.x(),wlh.x()) && !IsEqual(localP.z(),wlh.z())){
        return coord.j;
    }
    if(IsEqual(localP.z(),0) && !IsEqual(localP.x(),0) && !IsEqual(localP.y(),0) && !IsEqual(localP.x(),wlh.x()) && !IsEqual(localP.y(),wlh.y())){
        return -coord.k;
    }
    if(IsEqual(localP.z(),wlh.z()) && !IsEqual(localP.x(),0) && !IsEqual(localP.y(),0) && !IsEqual(localP.x(),wlh.x()) && !IsEqual(localP.y(),wlh.y())){
        return coord.k;
    }
    //其余的一律返回(0,0,0)
    Logger::logWarning("PhyBox::GetRandomNormal: The point is inside of the box, or on the edge");
    return Vector3d(0,0,0);
}

Eigen::Vector3d PhyBox::SupportPosition(const Eigen::Vector3d& direction)
{
    //将minP，maxP变为8个顶点坐标
    std::vector<Eigen::Vector3d> vertices = GetWorldVertices();
    double max = -1e9;
    Eigen::Vector3d result;
    for (auto& v : vertices)
    {
        double dot = v.dot(direction);
        if (dot > max)
        {
            max = dot;
            result = v;
        }
    }
    return result;
}


IntesectData PhyBox::CollisionWith(PhyObject* other)
{
    return IntesectData();
}