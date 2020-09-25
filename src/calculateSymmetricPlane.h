#pragma once

#include <string>
#include <vector>
#include "Eigen/Core"

template <typename Scalar, typename Index>
struct Mesh
{
    std::string meshName;
    std::vector<std::vector<Scalar>> V;
    std::vector<std::vector<Index>> F;
    std::vector<std::vector<std::pair<Scalar, Scalar>>> UV;
    std::vector<std::vector<Scalar>> VC;
    std::vector<Scalar> M;
    std::vector<Index> G;
};

template <typename Scalar, typename Index>
void calculateSymmetricPlane(
	const Mesh<Scalar, Index> &meshIn,
    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &eulerXZY,
    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &translateXYZ);

#include "calculateSymmetricPlane.cpp"
