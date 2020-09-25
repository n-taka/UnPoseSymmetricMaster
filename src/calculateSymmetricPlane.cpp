#ifndef calculateSymmetricPlane_CPP
#define calculateSymmetricPlane_CPP

#include "calculateSymmetricPlane.h"

#include "MeanShift.h"
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
        int index;
        Eigen::Matrix<double, 1, 3> xyz;
        Eigen::Matrix<double, 1, 3> c1, c2, n;
        double k1, k2;
    } Signature;

    typedef struct
    {
        double s;
        Eigen::Matrix<double, 1, 3> R;
        Eigen::Matrix<double, 1, 3> t;
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
            s.index = sample;
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

        const int numberOfNeighbors = 4;
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
                if (signQ.index == signP.index)
                {
                    // avoid pair with itself
                    continue;
                }
                if (std::fabs(((signQ.n - signP.n).normalized()).dot((signQ.xyz - signP.xyz).normalized())) < 0.5)
                {
                    // avoid pair with largely different normal
                    continue;
                }

                // calculate transformation
                {
                    Transformation t;
                    t.s = (signQ.k1 / signP.k1 + signQ.k2 / signP.k2) * 0.5;
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

                    transformation.push_back(t);
                    // sigPairs.push_back(std::make_pair(signQ.index, signP.index));
                }
                // calculate transformation (with reflection)
                {
                    Signature signPReflect = signP;
                    // reflect with YZ plane
                    signPReflect.xyz(0) *= -1;
                    signPReflect.c1(0) *= -1;
                    signPReflect.c2(0) *= -1;
                    signPReflect.n(0) *= -1;
                    // make sure (c1, c2, n) form right-hand coordinate
                    signPReflect.c1 *= -1;

                    Transformation t;
                    t.s = (signQ.k1 / signPReflect.k1 + signQ.k2 / signPReflect.k2) * 0.5;
                    const Eigen::Quaternion<double> q0 = Eigen::Quaternion<double>::FromTwoVectors(signQ.n, signPReflect.n);
                    Eigen::Matrix<double, 1, 3> rotatedQC1;
                    Eigen::Matrix<double, 1, 3> rotatedQC2;
                    // Eigen::Matrix<double, 1, 3> rotatedQN;
                    rotatedQC1 = q0 * signQ.c1;
                    rotatedQC2 = q0 * signQ.c2;
                    // rotatedQN = q0 * signQ.n;

                    const Eigen::Quaternion<double> q1 = Eigen::Quaternion<double>::FromTwoVectors(rotatedQC1, signPReflect.c1);
                    const Eigen::Quaternion<double> q2 = Eigen::Quaternion<double>::FromTwoVectors(rotatedQC1, -signPReflect.c1);
                    // q.w = cos(theta/2)
                    const Eigen::Matrix<double, 3, 3> R = ((((q1.w() > q2.w())) ? q1 : q2) * q0).normalized().toRotationMatrix();
                    t.R = R.eulerAngles(0, 1, 2).transpose();
                    t.t = signPReflect.xyz - t.s * (R * signQ.xyz.transpose()).transpose();

                    transformation.push_back(t);
                    // sigPairs.push_back(std::make_pair(signQ.index, signP.index));
                }
            }
        }
    }

    double kernel(double distance, double kernel_bandwidth)
    {
        double temp = exp(-1.0 / 2.0 * (distance * distance) / (kernel_bandwidth * kernel_bandwidth));
        return temp;
    }

    void clusteringTransformation(
        const std::vector<Transformation> &transformation,
        double &BBDiagonalLength)
    {
        const double beta1 = 0.01;
        const double beta2 = 1.0 / (M_PI * M_PI);
        const double beta3 = 4.0 / (BBDiagonalLength * BBDiagonalLength);
        const double windowSize = 3.0;

        // double kernel_bandwidth = 3;
        //     vector<double> weights = {beta_1, beta_2, beta_2, beta_2, beta_3, beta_3, beta_3, 0, 0};
        //     MeanShift *msp = new MeanShift(epanechnikov_kernel, weights);

        //     vector<vector<double> > points;
        //     Transformation::to_points(transf_space, points);
        //     vector<Cluster> clusters = msp->cluster(points, kernel_bandwidth);
        //     for (Cluster &cluster : clusters) {
        //         vector<Transformation> cluster_transf;
        //         for (vector<double> &point : cluster.original_points)
        //             cluster_transf.push_back(Transformation(point));
        //         clusters_transf.push_back(cluster_transf);
        //     }

        auto norm = [&](const Transformation &t) {
            double val = 0.0;
            val += beta1 * t.s * t.s;
            val += beta2 * t.R.squaredNorm();
            val += beta3 * t.t.squaredNorm();
            return val;
        };

        auto kernel = [&](const Transformation &t) {
            // see https://en.wikipedia.org/wiki/Volume_of_an_n-ball
            const double unit7DSphereVol = 16.0 * M_PI * M_PI * M_PI / 105.0;
            const double normT = norm(T);
            if (normT <= 1)
            {
                return 0.5 * (1.0 / unit7DSphereVol) * (7.0 + 2.0) * (1.0 - norm(t));
            }
            else
            {
                return 0.0;
            }
        };
    }

} // namespace

template <typename Scalar, typename Index>
void calculateSymmetricPlane(
    const Mesh<Scalar, Index> &meshIn,
    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &point,
    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &normal)
{
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> V;
    Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> F, F_;
    igl::list_to_matrix(meshIn.V, V);
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
    F.resize(maskedF.size(), 3);
    for (int f = 0; f < F.rows(); ++f)
    {
        F.row(f) = F_.row(maskedF.at(f));
    }
    Eigen::Matrix<int, Eigen::Dynamic, 1> I;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> V_ = V;
    F_ = F;
    igl::remove_unreferenced(V_, F_, V, F, I);
    igl::writeOBJ("masked.obj", V, F);

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
    std::cout << "Use " << signature.size() << " samples" << std::endl;

    ////
    // (Sec. 2.1) Point Pruning
    ////
    prunePoints(signature, 0.75);
    std::cout << "Use " << signature.size() << " samples (after pruning)" << std::endl;

    ////
    // (Sec. 2.2) Pairing
    ////
    // std::vector<std::pair<int, int>> sigPairs;
    std::vector<Transformation> transformation;
    pairPoints(signature, 100, transformation);
    std::cout << "Use " << transformation.size() << " pairs (incl. reflection)" << std::endl;

    // {
    //     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmpV;
    //     Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> tmpF, tmpE;
    //     tmpV.resize(sigPairs.size() * 2, 3);
    //     tmpE.resize(sigPairs.size(), 2);
    //     for (int v = 0; v < sigPairs.size(); ++v)
    //     {
    //         tmpV.row(v * 2 + 0) = V.row(sigPairs.at(v).first);
    //         tmpV.row(v * 2 + 1) = V.row(sigPairs.at(v).second);
    //         tmpE.row(v) << v * 2 + 0, v * 2 + 1;
    //     }
    //     std::cout << std::endl;
    //     igl::writePLY("pairs.ply", tmpV, tmpF, tmpE);
    // }

    // (Sec. 3) Clustering
    const double BBDiagonalLength = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
    clusteringTransformation(transformation, BBDiagonalLength);
    /// TODO

    // (Sec. 4) Verification
    /// TODO

    // (Sec. 2.2)display patches or write to file
    /// TODO
}

#endif
