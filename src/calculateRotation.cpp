#ifndef calculateRotation_CPP
#define calculateRotation_CPP

#include "calculateRotation.h"

#include "Eigen/Geometry"

#pragma warning(push)
#pragma warning(disable : 4018 4101 4129 4244 4267 4305 4566 4819 4996)
#include "igl/list_to_matrix.h"
#include "igl/polygon_mesh_to_triangle_mesh.h"
#include "igl/remove_unreferenced.h"
#include "igl/iterative_closest_point.h"
#pragma warning(pop)

template <typename Scalar, typename Index>
void calculateRotation(
    const Mesh<Scalar, Index> &meshIn,
    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &eulerXZY,
    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &translateXYZ)
{

    // convert to matrix
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> VX, VY;
    Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> FX, FY;
    igl::list_to_matrix(meshIn.V, VX);
    VY = VX;
    igl::polygon_mesh_to_triangle_mesh(meshIn.F, FY);

    std::vector<int> maskedF;
    for (int f = 0; f < FY.rows(); ++f)
    {
        bool masked = true;
        for (int fv = 0; fv < 3; ++fv)
        {
            masked &= (meshIn.M.at(FY(f, fv)) < 0.5);
        }
        if (masked)
        {
            maskedF.push_back(f);
        }
    }
    FX.resize(maskedF.size(), 3);
    for (int f = 0; f < FX.rows(); ++f)
    {
        FX.row(f) = FY.row(maskedF.at(f));
    }
    // mirror with YZ plane
    VY.col(0) *= -1;
    for (int f = 0; f < FY.rows(); ++f)
    {
        std::swap(FY(f, 1), FY(f, 2));
    }

    Eigen::Matrix<int, Eigen::Dynamic, 1> I;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> VX_ = VX;
    Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> FX_ = FX;
    igl::remove_unreferenced(VX_, FX_, VX, FX, I);
    igl::writeOBJ("masked.obj", VX, FX);

    Eigen::Matrix<double, 3, 3> R;
    Eigen::Matrix<double, 1, 3> t;
    igl::iterative_closest_point(VX, FX, VY, FY, 1000, 100, R, t);

    eulerXZY = R.eulerAngles(0, 2, 1);
    translateXYZ = t;
}

#endif
