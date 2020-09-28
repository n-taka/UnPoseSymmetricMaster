#ifndef calculateEulerAnglesForSymmetrize_CPP
#define calculateEulerAnglesForSymmetrize_CPP

#include "calculateEulerAnglesForSymmetrize.h"

#include "Meanshift.hpp"
#include "Eigen/Geometry"

#pragma warning(push)
#pragma warning(disable : 4018 4101 4129 4244 4267 4305 4566 4819 4996)
#include "igl/list_to_matrix.h"
#include "igl/polygon_mesh_to_triangle_mesh.h"
#include "igl/remove_unreferenced.h"

#include "igl/parallel_for.h"
#include "igl/per_vertex_normals.h"
#include "igl/principal_curvature.h"

#include "igl/octree.h"
#include "igl/knn.h"

#include "igl/AABB.h"
#include "igl/per_face_normals.h"
#include "igl/iterative_closest_point.h"

#include "igl/writeOBJ.h"
#include "igl/writePLY.h"
#include "igl/writeSTL.h"
#pragma warning(pop)

#include <algorithm>
#include <random>

namespace
{
    typedef struct
    {
        int vIndex;
        Eigen::Matrix<double, 1, 3> xyz;
        Eigen::Matrix<double, 1, 3> c1, c2, n;
        double k1, k2;
    } Signature;

    typedef struct
    {
        double s;
        Eigen::Matrix<double, 1, 3> R;
        Eigen::Matrix<double, 1, 3> t;
        int vi, vj;
    } Transformation;

    void calculateSignatureFromMesh(
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &V,
        const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &F,
        const std::vector<int> &sampleIdxIntoV,
        std::vector<Signature> &signature)
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> N;
        igl::per_vertex_normals(V, F, N);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> PDmax, PDmin;
        Eigen::Matrix<double, Eigen::Dynamic, 1> PVmax, PVmin;
        igl::principal_curvature(V, F, PDmax, PDmin, PVmax, PVmin);

        signature.reserve(sampleIdxIntoV.size());
        for (int sample = 0; sample < sampleIdxIntoV.size(); ++sample)
        {
            Signature s;
            s.vIndex = sampleIdxIntoV.at(sample);
            s.xyz = V.row(sampleIdxIntoV.at(sample)).template head<3>();
            s.k1 = PVmin(sampleIdxIntoV.at(sample), 0);
            s.k2 = PVmax(sampleIdxIntoV.at(sample), 0);
            s.c1 = PDmin.row(sampleIdxIntoV.at(sample)).template head<3>();
            s.c2 = PDmax.row(sampleIdxIntoV.at(sample)).template head<3>();
            s.n = s.c1.cross(s.c2);
            s.n.normalize();
            if (s.n.dot(N.row(sampleIdxIntoV.at(sample)).template head<3>()) < 0.0)
            {
                s.n *= -1;
                s.c1 *= -1;
            }
            signature.push_back(s);
        }
    }

    void prunePoints(
        std::vector<Signature> &signature,
        const double &gamma)
    {
        std::vector<Signature> tmp = signature;
        signature.clear();
        signature.reserve(tmp.size());
        for (const auto &s : tmp)
        {
            if (std::fabs(s.k1 / s.k2) < gamma)
            {
                signature.push_back(s);
            }
        }
    }

    void pairPoints(
        const std::vector<Signature> &signature,
        const int &maxSampleCount,
        std::vector<Transformation> &transformation)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 3> P;
        P.resize(signature.size(), 3);
        for (int s = 0; s < signature.size(); ++s)
        {
            const auto &sign = signature.at(s);
            // For our purpose, we use Omega_6 for pairing
            P.row(s) << sign.k1, sign.k2, 0;
        }
        std::vector<std::vector<int>> point_indices;
        Eigen::Matrix<int, Eigen::Dynamic, 8> CH;
        Eigen::Matrix<double, Eigen::Dynamic, 3> CN;
        Eigen::Matrix<double, Eigen::Dynamic, 1> W;
        igl::octree(P, point_indices, CH, CN, W);

        const int numberOfNeighbors = std::min<int>(4, P.rows());
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> I;
        igl::knn(P, numberOfNeighbors, point_indices, CH, CN, W, I);

        ////
        // Sampling.
        ////
        std::vector<int> sampledSigIdx;
        const int sampleCount = std::min<int>(maxSampleCount, signature.size());
        {
            sampledSigIdx.resize(sampleCount);
            std::iota(sampledSigIdx.begin(), sampledSigIdx.end(), 0);
            std::random_device rd;
            std::mt19937 g(rd());

            std::shuffle(sampledSigIdx.begin(), sampledSigIdx.end(), g);

            sampledSigIdx.erase(sampledSigIdx.begin() + sampleCount, sampledSigIdx.end());
        }

        for (int s = 0; s < sampledSigIdx.size(); ++s)
        {
            const int &sampleIdx = sampledSigIdx.at(s);
            const auto &signQ = signature.at(sampleIdx);
            for (int n = 0; n < numberOfNeighbors; ++n)
            {
                const auto &signP = signature.at(I(sampleIdx, n));
                if (signQ.vIndex == signP.vIndex)
                {
                    // avoid pair with itself
                    continue;
                }
                if (std::fabs(((signQ.n - signP.n).normalized()).dot((signQ.xyz - signP.xyz).normalized())) < 0.5)
                {
                    // avoid pair with largely different normal
                    continue;
                }

#if 0
                // calculate transformation
                {
                    Transformation t;
                    t.s = (signQ.k1 / signP.k1 + signQ.k2 / signP.k2) * 0.5;
                    // t.s = 1.0;
                    const Eigen::Quaternion<double> q0 = Eigen::Quaternion<double>::FromTwoVectors(signQ.n, signP.n);
                    Eigen::Matrix<double, 1, 3> rotatedQC1;
                    Eigen::Matrix<double, 1, 3> rotatedQC2;
                    // Eigen::Matrix<double, 1, 3> rotatedQN;
                    rotatedQC1 = q0 * signQ.c1;
                    rotatedQC2 = q0 * signQ.c2;
                    // rotatedQN = q0 * signQ.n;

                    const Eigen::Quaternion<double> q1 = Eigen::Quaternion<double>::FromTwoVectors(rotatedQC1, signP.c1);
                    const Eigen::Quaternion<double> q2 = Eigen::Quaternion<double>::FromTwoVectors(rotatedQC1, -signP.c1);
                    // q.w = cos(theta/2)
                    const Eigen::Matrix<double, 3, 3> R = ((((q1.w() > q2.w())) ? q1 : q2) * q0).normalized().toRotationMatrix();
                    t.R = R.eulerAngles(0, 1, 2).transpose();
                    t.t = signP.xyz - t.s * (R * signQ.xyz.transpose()).transpose();

                    t.vi = signQ.vIndex;
                    t.vj = signP.vIndex;

                    transformation.push_back(t);
                    // sigPairs.push_back(std::make_pair(signQ.index, signP.index));
                }
#endif
                // calculate transformation (with reflection)
                {
                    Signature signQReflect = signQ;
                    Signature signPReflect = signP;
                    // reflect with YZ plane
                    signQReflect.xyz(0) *= -1;
                    signQReflect.c1(0) *= -1;
                    signQReflect.c2(0) *= -1;
                    signQReflect.n(0) *= -1;
                    // make sure (c1, c2, n) form right-hand coordinate
                    signQReflect.c1 *= -1;
                    // reflect with YZ plane
                    signPReflect.xyz(0) *= -1;
                    signPReflect.c1(0) *= -1;
                    signPReflect.c2(0) *= -1;
                    signPReflect.n(0) *= -1;
                    // make sure (c1, c2, n) form right-hand coordinate
                    signPReflect.c1 *= -1;

                    Transformation t;
                    // t.s = 1.0;
                    t.s = (signQReflect.k1 / signPReflect.k1 + signQReflect.k2 / signPReflect.k2) * 0.5;
                    const Eigen::Quaternion<double> q0 = Eigen::Quaternion<double>::FromTwoVectors(signQReflect.n, signPReflect.n);
                    Eigen::Matrix<double, 1, 3> rotatedQC1;
                    Eigen::Matrix<double, 1, 3> rotatedQC2;
                    // Eigen::Matrix<double, 1, 3> rotatedQN;
                    rotatedQC1 = q0 * signQReflect.c1;
                    rotatedQC2 = q0 * signQReflect.c2;
                    // rotatedQN = q0 * signQ.n;

                    const Eigen::Quaternion<double> q1 = Eigen::Quaternion<double>::FromTwoVectors(rotatedQC1, signPReflect.c1);
                    const Eigen::Quaternion<double> q2 = Eigen::Quaternion<double>::FromTwoVectors(rotatedQC1, -signPReflect.c1);
                    // q.w = cos(theta/2)
                    const Eigen::Matrix<double, 3, 3> R = ((((q1.w() > q2.w())) ? q1 : q2) * q0).normalized().toRotationMatrix();
                    t.R = R.eulerAngles(0, 1, 2).transpose();
                    t.t = signPReflect.xyz - t.s * (R * signQReflect.xyz.transpose()).transpose();

                    t.vi = signQ.vIndex;
                    t.vj = signP.vIndex;

                    transformation.push_back(t);
                    // sigPairs.push_back(std::make_pair(signQ.index, signP.index));
                }
            }
        }
    }

    void clusteringTransformation(
        const std::vector<Transformation> &transformation,
        const double &BBDiagonalLength,
        const int &minModeCount,
        const int &maxModeCount,
        std::vector<Transformation> &modes,
        std::vector<std::vector<int>> &indices)
    {
        const double beta1 = 0.01;
        const double beta2 = 1.0 / (M_PI * M_PI);
        const double beta3 = 4.0 / (BBDiagonalLength * BBDiagonalLength);

        std::vector<std::vector<double>> X;
        X.reserve(transformation.size());
        for (const auto &t : transformation)
        {
            X.push_back({t.s, t.R(0, 0), t.R(0, 1), t.R(0, 2), t.t(0, 0), t.t(0, 1), t.t(0, 2)});
        }
        std::vector<std::vector<double>> modesV;
        std::vector<int> indexMap;
        double kernel_bandwidth = 2.0;
        while (modesV.size() > maxModeCount || modesV.size() < minModeCount)
        {
            clustering::Meanshift meanShift(0, kernel_bandwidth, 7, 0.5, beta1, beta2, beta3);
            meanShift.FindModes(X, modesV, indexMap);
            if (modesV.size() > maxModeCount)
            {
                kernel_bandwidth *= 1.5;
            }
            if (modesV.size() < minModeCount)
            {
                kernel_bandwidth *= 0.8;
            }
        }

        std::vector<std::pair<int, int>> count;
        for (int i = 0; i < modesV.size(); ++i)
        {
            count.push_back(std::make_pair(0, i));
        }
        std::vector<std::vector<int>> indicesV;
        indicesV.resize(modesV.size());
        for (int i = 0; i < indexMap.size(); ++i)
        {
            const int &index = indexMap.at(i);
            count.at(index).first -= 1;
            indicesV.at(index).push_back(i);
        }
        std::sort(count.begin(), count.end());

        // in decending order
        modes.clear();
        modes.reserve(modesV.size());
        indices.clear();
        indices.reserve(modesV.size());
        for (int m = 0; m < count.size(); ++m)
        {
            const std::vector<double> &mode = modesV.at(count.at(m).second);
            Transformation t;
            t.s = mode.at(0);
            t.R << mode.at(1), mode.at(2), mode.at(3);
            t.t << mode.at(4), mode.at(5), mode.at(6);
            modes.push_back(t);
            indices.push_back(indicesV.at(count.at(m).second));
        }
    }

    void estimateBestPlane(
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &V,
        const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &F,
        const std::vector<Transformation> &transformation,
        const std::vector<int> &indices,
        Eigen::Matrix<double, 1, Eigen::Dynamic> &point,
        Eigen::Matrix<double, 1, Eigen::Dynamic> &normal)
    {
        const double threshold = 1.0e-1;
        double minError = std::numeric_limits<double>::max();

        igl::AABB<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, 3> Ytree;
        Ytree.init(V, F);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> NY;
        igl::per_face_normals(V, F, NY);

        Eigen::Matrix<double, 3, 3> bestR;
        Eigen::Matrix<double, 1, 3> bestt;
        Eigen::Matrix<double, 1, Eigen::Dynamic> bestPoint;
        Eigen::Matrix<double, 1, Eigen::Dynamic> bestNormal;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mirrorICPV;

        for (int pairIdx = 0; pairIdx < indices.size() && minError > threshold; ++pairIdx)
        {
            const Transformation &oneOfBestT = transformation.at(indices.at(pairIdx));
            const Eigen::Matrix<double, 1, Eigen::Dynamic> vi = V.row(oneOfBestT.vi);
            const Eigen::Matrix<double, 1, Eigen::Dynamic> vj = V.row(oneOfBestT.vj);
            point = (vi + vj) * 0.5;
            normal = (vi - vj).normalized();

            // alignment with ICP
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> VX;
            Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> FX;

            // mirror with plane
            VX = V;
            for (int v = 0; v < VX.rows(); ++v)
            {
                const double height = (VX.row(v) - point).dot(normal);
                VX.row(v) -= height * 2 * normal;
            }
            FX = F;
            for (int f = 0; f < FX.rows(); ++f)
            {
                std::swap(FX(f, 1), FX(f, 2));
            }
            // igl::writeOBJ("mirror.obj", VX, FX);

            const int num_samples = 2000;
            const int max_iters = 20;

            Eigen::Matrix<double, 3, 3> R;
            Eigen::Matrix<double, 1, 3> t;
            R.setIdentity(3, 3);
            t.setConstant(1, 3, 0);
            Eigen::Matrix<double, Eigen::Dynamic, 1> sqrD;
            for (int iter = 0; iter < max_iters; iter++)
            {
                // Sample random points on X
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X;
                {
                    Eigen::VectorXi XI;
                    Eigen::MatrixXd B;
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> VXRT = (VX * R).rowwise() + t;

                    igl::random_points_on_mesh(num_samples, VXRT, FX, B, XI, X);
                }
                // Compute closest point
                Eigen::VectorXi I;
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P;
                {
                    Ytree.squared_distance(V, F, X, sqrD, I, P);
                }

                // Use better normals?
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> N;
                igl::slice(NY, I, 1, N);
                //MatrixXS N = (X - P).rowwise().normalized();
                // fit rotation,translation
                Eigen::Matrix<double, 3, 3> Rup;
                Eigen::Matrix<double, 1, 3> tup;
                // Note: Should try out Szymon Rusinkiewicz's new symmetric icp
                igl::rigid_alignment(X, P, N, Rup, tup);
                // update running rigid transformation
                R = (R * Rup).eval();
                t = (t * Rup + tup).eval();
                // Better stopping condition?
            }

            if ((sqrD.sum() / num_samples) < minError)
            {
                minError = sqrD.sum() / num_samples;
                bestR = R;
                bestt = t;
                bestPoint = point;
                bestNormal = normal;

                VX *= R;
                VX.rowwise() += t;
                mirrorICPV = VX;
                // igl::writeOBJ("mirror_ICP_best.obj", VX, FX);
            }
        }

        point = (V.row(0) + mirrorICPV.row(0)) * 0.5;
        normal = (V.row(0) - mirrorICPV.row(0)).normalized();
        // {
        //     Eigen::Matrix<double, 1, Eigen::Dynamic> base1 = normal.template head<3>().cross(V.row(0).template head<3>()).normalized();
        //     Eigen::Matrix<double, 1, Eigen::Dynamic> base0 = base1.template head<3>().cross(normal.template head<3>()).normalized();
        //     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmpV;
        //     Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> tmpF;
        //     tmpV.resize(4, 3);
        //     tmpV.row(0) = point + base0 * 10;
        //     tmpV.row(1) = point + base1 * 10;
        //     tmpV.row(2) = point - base0 * 10;
        //     tmpV.row(3) = point - base1 * 10;
        //     tmpF.resize(2, 3);
        //     tmpF.row(0) << 0, 1, 2;
        //     tmpF.row(1) << 2, 3, 0;
        //     igl::writeOBJ("plane_ICP_best.obj", tmpV, tmpF);
        // }
    }
} // namespace

bool calculateEulerAnglesForSymmetrize(
    const Mesh<double, int> &meshIn,
    Eigen::Matrix<double, 3, 3> &rotMatrix,
    Eigen::Matrix<double, 1, 3> &translateXYZ)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> V, VRaw;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> F, F_;
    igl::list_to_matrix(meshIn.V, VRaw);
    // avoid too small input
    double ratio = 100.0 / (VRaw.colwise().maxCoeff() - VRaw.colwise().minCoeff()).norm();
    VRaw *= ratio;

    igl::polygon_mesh_to_triangle_mesh(meshIn.F, F_);

    std::vector<int> maskedF;
    for (int f = 0; f < F_.rows(); ++f)
    {
        bool masked = true;
        for (int fv = 0; fv < 3; ++fv)
        {
            masked &= (meshIn.M.at(F_(f, fv)) < 0.5);
        }
        if (masked)
        {
            maskedF.push_back(f);
        }
    }
    if (maskedF.size() <= 0)
    {
        return false;
    }

    F.resize(maskedF.size(), 3);
    for (int f = 0; f < F.rows(); ++f)
    {
        F.row(f) = F_.row(maskedF.at(f));
    }
    Eigen::Matrix<int, Eigen::Dynamic, 1> I;
    F_ = F;
    igl::remove_unreferenced(VRaw, F_, V, F, I);
    // igl::writeOBJ("masked.obj", V, F);

    ////
    // find symmetry plane
    ////

    ////
    // Sampling.
    ////
    std::vector<int> sampleIdxIntoV;
    {
        // currently I use all the vertices as samples
        sampleIdxIntoV.reserve(V.rows());
        for (int i = 0; i < V.rows(); ++i)
        {
            sampleIdxIntoV.push_back(i);
        }
    }

    ////
    // (Sec. 2) Signatures and Transformations
    ////
    std::vector<Signature> signature;
    calculateSignatureFromMesh(V, F, sampleIdxIntoV, signature);
    // std::cout << "Use " << signature.size() << " samples" << std::endl;
    if (signature.size() <= 0)
    {
        return false;
    }

    ////
    // (Sec. 2.1) Point Pruning
    ////
    prunePoints(signature, 0.75);
    // std::cout << "Use " << signature.size() << " samples (after pruning)" << std::endl;
    if (signature.size() <= 0)
    {
        return false;
    }

    ////
    // (Sec. 2.2) Pairing
    ////
    // std::vector<std::pair<int, int>> sigPairs;
    std::vector<Transformation> transformation;
    pairPoints(signature, 200, transformation);
    // std::cout << "Use " << transformation.size() << " pairs (incl. reflection)" << std::endl;
    if (transformation.size() <= 0)
    {
        return false;
    }
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmpV;
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> tmpF, tmpE;
        tmpV.resize(transformation.size() * 2, 3);
        tmpE.resize(transformation.size(), 2);
        for (int v = 0; v < transformation.size(); ++v)
        {
            tmpV.row(v * 2 + 0) = V.row(transformation.at(v).vi);
            tmpV.row(v * 2 + 1) = V.row(transformation.at(v).vj);
            tmpE.row(v) << v * 2 + 0, v * 2 + 1;
        }
        // igl::writePLY("pairs.ply", tmpV, tmpF, tmpE);
    }

    // (Sec. 3) Clustering
    const double BBDiagonalLength = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
    const int minModeCount = 5;
    const int maxModeCount = 20;
    std::vector<Transformation> modes;
    std::vector<std::vector<int>> indices;
    clusteringTransformation(transformation, BBDiagonalLength, minModeCount, maxModeCount, modes, indices);
    // std::cout << modes.size() << " modes are found" << std::endl;
    if (modes.size() <= 0 || indices.at(0).size() <= 0)
    {
        return false;
    }

    // (Sec. 4) Verification
    Eigen::Matrix<double, 1, Eigen::Dynamic> point, normal;
    estimateBestPlane(V, F, transformation, indices.at(0), point, normal);

    // move to center
    const Eigen::Quaternion<double> q = Eigen::Quaternion<double>::FromTwoVectors(Eigen::Vector3d::UnitX(), normal);

    rotMatrix = q.normalized().toRotationMatrix();

    VRaw *= rotMatrix;
    VRaw /= ratio;
    translateXYZ = -(VRaw.colwise().maxCoeff() + VRaw.colwise().minCoeff()) * 0.5;

    point /= ratio;
    point *= rotMatrix;
    translateXYZ(0) = -point(0);

    return true;
}

#endif
