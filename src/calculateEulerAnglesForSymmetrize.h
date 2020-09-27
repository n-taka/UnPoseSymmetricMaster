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

bool calculateEulerAnglesForSymmetrize(
	const Mesh<double, int> &meshIn,
    Eigen::Matrix<double, 3, 3> &rotMatrix,
    double &reflectionX);

#include "calculateEulerAnglesForSymmetrize.cpp"
