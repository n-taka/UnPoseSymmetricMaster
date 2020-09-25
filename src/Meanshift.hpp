#pragma once
// https://gist.github.com/ialhashim/b29e5455333aa6ae0071#file-meanshift-hpp

/* Meanshift non-parametric mode estimator.
 *
 * Code adapted from : Sebastian Nowozin <nowozin@gmail.com>
 */

#define _USE_MATH_DEFINES
#include <math.h>

#include <vector>
#include <map>
#include <iostream>

//#include <gmm/gmm.h>
#include <Eigen/Core>

namespace gmm
{
    template <typename V>
    void clear(V &v)
    {
        for (size_t i = 0; i < v.size(); i++)
            v[i] = 0;
    }

    template <typename V, typename W>
    void copy(const V &v, W &w)
    {
        for (size_t i = 0; i < v.size(); i++)
            w[i] = v[i];
    }

    template <typename V, typename W>
    void add(const V &v, W &w)
    {
        for (size_t i = 0; i < w.size(); i++)
            w[i] += v[i];
    }

    template <typename V>
    V scaled(const V &v, double s)
    {
        V a = v;
        for (size_t i = 0; i < a.size(); i++)
            a[i] *= s;
        return a;
    }

    template <typename V>
    void scale(V &v, double s)
    {
        for (size_t i = 0; i < v.size(); i++)
            v[i] *= s;
    }

    template <typename V, typename W>
    double vect_dist2(V &v, W &w, const double &beta1, const double &beta2, const double &beta3)
    {
        Eigen::VectorXd a = Eigen::Map<Eigen::VectorXd>(&v[0], v.size());
        Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(&w[0], w.size());
        Eigen::VectorXd diff = a - b;
        double val = 0.0;
        val += beta1 * diff.block(0, 0, 1, 1).squaredNorm();
        val += beta2 * diff.block(1, 0, 1, 3).squaredNorm();
        val += beta3 * diff.block(4, 0, 1, 3).squaredNorm();
        return val;
        // return (a - b).lpNorm<2>();
    }

    template <typename V>
    double vect_norm2(V &v, const double &beta1, const double &beta2, const double &beta3)
    {
        Eigen::VectorXd a = Eigen::Map<Eigen::VectorXd>(&v[0], v.size());
        double val = 0.0;
        val += beta1 * a.block(0, 0, 1, 1).squaredNorm();
        val += beta2 * a.block(1, 0, 1, 3).squaredNorm();
        val += beta3 * a.block(4, 0, 1, 3).squaredNorm();
        return val;
        // return Eigen::Map<Eigen::VectorXd>(&v[0], v.size()).norm();
    }
} // namespace gmm

namespace clustering
{

    /* We use the notation of the description of Mean Shift in
 *   [Comaniciu2000], Dorin Comaniciu, Peter Meer,
 *   "Mean Shift: A Robust Approach toward Feature Space Analysis"
 *
 * The implementation is naive and does not exploit fast nearest neighbor
 * lookups or triangle inequalities to speed up the mean shift procedure.
 * Therefore, it is only suitable for low-dimensional input spaces (say, <=
 * 10) with relatively few samples (say, < 1e6).
 */
    class Meanshift
    {
    public:
        /* kernel_type selects the kernel profile:
     *   0 for the Epanechnikov kernel profile
     *     k_E(x) = (1 - x) if (0 <= x <= 1), 0 otherwise.
     *   1 for the truncated multivariate normal kernel profile
     *     k_N(x) = exp(-0.5 * x)
     * kernel_bandwidth: The positive bandwidth parameter.
     * mode_tolerance: mode matching tolerance.  Modes which have a L2
     *   distance closer than this value will be treated as being the same.
     */
        Meanshift(int kernel_type, double kernel_bandwidth,
                  int dim, double mode_tolerance, double beta1, double beta2, double beta3);

        /* Find modes of an empirical distribution X.
     *
     * X: N vectors of size M, representing N samples from the distribution.
     * modes: Output, will be allocated properly.  Return all modes found.
     * indexmap: N-vector of absolute mode indices.  Each sample point is
     *   assigned to one mode.  The vector must already be properly sized.
     * procedurecount: The mean shift procedure will be run at least this many
     *   times, sampling randomly from X as initialization.  If zero, all
     *   samples from X are used as initializations.
     */
        void FindModes(const std::vector<std::vector<double>> &X,
                       std::vector<std::vector<double>> &modes,
                       std::vector<int> &indexmap,
                       unsigned int procedure_count = 0) const;

        /* Mean Shift Procedure, starting from mode, perform mean shift on the
     * distribution empirically sampled in X.
     *
     * X: N vectors of size M, representing N samples from the distribution.
     * mode: Starting point (for example a point from X).  The result will be
     *    given in mode.
     * visited: N-vector of indicator variables.  If the mean shift procedure
     *    passes through a point in X, the corresponding index in visited will
     *    be set to one.
     *
     * Return the number of iterations used.
     */
        int MeanshiftProcedure(const std::vector<std::vector<double>> &X,
                               std::vector<double> &mode, std::vector<unsigned int> &visited) const;

    private:
        int dim;         // input space dimension
        int kernel_type; // 0: Epanechnikov, 1: truncated multivariate normal
        double kernel_c; // constant normalization factor
        double kernel_bandwidth;
        double mode_tolerance;
        double beta1;
        double beta2;
        double beta3;
    };

    Meanshift::Meanshift(int kernel_type, double kernel_bandwidth,
                         int dim, double mode_tolerance, double beta1, double beta2, double beta3)
        : dim(dim), kernel_type(kernel_type),
          kernel_bandwidth(kernel_bandwidth),
          mode_tolerance(mode_tolerance),
          beta1(beta1),
          beta2(beta2),
          beta3(beta3)
    {
        assert(kernel_type >= 0 && kernel_type <= 1);
        assert(kernel_bandwidth > 0.0);
        assert(dim > 0);

        // Compute normalization constant
        if (kernel_type == 0)
        {
            // see https://en.wikipedia.org/wiki/Volume_of_an_n-ball
            const double unit7DSphereVol = 16.0 * M_PI * M_PI * M_PI / 105.0;
            kernel_c = 0.5 * (static_cast<double>(dim) + 2.0) / unit7DSphereVol;
            // TODO: replace 1.234 with the volume of the d-dimensional unit
            // sphere FIXME
        }
        else if (kernel_type == 1)
        {
            kernel_c = pow(2.0 * M_PI, -static_cast<double>(dim) / 2.0);
        }
    }

    void Meanshift::FindModes(const std::vector<std::vector<double>> &X,
                              std::vector<std::vector<double>> &modes,
                              std::vector<int> &indexmap,
                              unsigned int procedure_count) const
    {
        //assert(indexmap.size() == X.size());
        indexmap.clear();
        indexmap.resize(X.size(), 0);

        modes.clear();

        if (procedure_count != 0 && procedure_count >= X.size())
            procedure_count = 0;

        // Perform dense mode finding: for each sample in X, perform the mean
        // shift procedure
        std::vector<double> mode(X[0].size());
        std::vector<unsigned int> visited(X.size());
        // [mode_index][sample_index] = number of times sample_index was visited
        //   when arriving at the mode with the mode_index.
        std::vector<std::map<unsigned int, unsigned int>> visit_counts;
        for (unsigned int si = 0; (procedure_count == 0 && si < X.size()) || (procedure_count != 0 && si < procedure_count); ++si)
        {
            if (procedure_count == 0)
            {
                gmm::copy(X[si], mode);
            }
            else
            {
                // Pick a random one from X
                unsigned sample_index = rand() % X.size();
                gmm::copy(X[sample_index], mode);
            }
            // std::cout << "Sample " << si << " of "
            //           << (procedure_count == 0 ? X.size() : procedure_count)
            //           << std::endl;
            gmm::clear(visited);
            MeanshiftProcedure(X, mode, visited);

            // Identify whether this is a novel mode or a known one
            bool found = false;
            unsigned int M_cur = 0;
            for (std::vector<std::vector<double>>::iterator Mi =
                     modes.begin();
                 found == false && Mi != modes.end(); ++Mi)
            {
                if (gmm::vect_dist2(*Mi, mode, beta1, beta2, beta3) < mode_tolerance)
                {
                    M_cur = Mi - modes.begin();
                    found = true;
                }
            }
            if (found == false)
            {
                // Add novel mode
                modes.push_back(mode);
                // TODO: remove this debug output
                // std::cout << modes.size() << " mode"
                //           << (modes.size() >= 2 ? "s" : "") << std::endl;

                // Add a new mapping of which samples have been visited while
                // approaching the novel mode.
                std::map<unsigned int, unsigned int> mode_visits;
                for (std::vector<unsigned int>::const_iterator vi =
                         visited.begin();
                     vi != visited.end(); ++vi)
                {
                    if (*vi != 0)
                        mode_visits[vi - visited.begin()] = 1;
                }
                visit_counts.push_back(mode_visits);
            }
            else
            {
                // The mode has been known, but we maybe crossed old and new
                // samples.  Update the counts.
                for (std::vector<unsigned int>::const_iterator vi =
                         visited.begin();
                     vi != visited.end(); ++vi)
                {
                    if (*vi != 0)
                        visit_counts[M_cur][vi - visited.begin()] += 1;
                }
            }
        }
#ifdef DEBUG
        std::cout << "Found " << modes.size() << " modes." << std::endl;
#endif

        // Perform index mapping: each sample gets assigned to one mode index.
        unsigned int unmapped_count = 0;
        for (unsigned int sample_index = 0;
             sample_index < X.size(); ++sample_index)
        {
            // Find mode index with highest count
            unsigned int maximum_count = 0;
            int maximum_mode_index = -1;
            for (unsigned int mode_index = 0;
                 mode_index < modes.size(); ++mode_index)
            {
                if (visit_counts[mode_index].count(sample_index) == 0)
                    continue;

                unsigned int count_cur = visit_counts[mode_index][sample_index];
                // std::cout << "  mode " << mode_index << " visited "
                //           << "sample " << sample_index << " for "
                //           << count_cur << " times." << std::endl;
                if (count_cur > maximum_count)
                {
                    maximum_mode_index = mode_index;
                    maximum_count = count_cur;
                }
            }
            if (maximum_mode_index == -1)
                unmapped_count += 1;
            indexmap[sample_index] = maximum_mode_index;
        }
        // std::cout << unmapped_count << " unmapped samples." << std::endl;
    }

    int Meanshift::MeanshiftProcedure(const std::vector<std::vector<double>> &X,
                                      std::vector<double> &mode, std::vector<unsigned int> &visited) const
    {
        assert(X.size() > 0);
        assert(visited.size() == X.size());
        assert(X[0].size() == mode.size());
        std::vector<double> meanshift(mode.size(), 0);
        std::vector<double> meanshift_cur(mode.size(), 0);
        int iter = 0;
        bool converged = false;
        do
        {
            // Compute the mean shift vector
            gmm::clear(meanshift);
            double denominator = 0.0;
            for (std::vector<std::vector<double>>::const_iterator Xi = X.begin();
                 Xi != X.end(); ++Xi)
            {
                gmm::copy(mode, meanshift_cur);
                gmm::add(gmm::scaled(*Xi, -1.0), meanshift_cur);
                gmm::scale(meanshift_cur, 1.0 / kernel_bandwidth);
                double weight_cur = gmm::vect_norm2(meanshift_cur, beta1, beta2, beta3);
                if (kernel_type == 0)
                {
                    // Epanechnikov kernel is the shadow of the uniform kernel
                    if (weight_cur <= 1.0)
                        weight_cur = 1.0;
                    else
                        weight_cur = 0.0;
                }
                else if (kernel_type == 1)
                {
                    // Multivariate normal kernel is the shadow of itself
                    weight_cur = kernel_c * exp(-0.5 * weight_cur);

                    // Truncate
                    if (weight_cur < 1e-2)
                        weight_cur = 0.0;
                }
                if (weight_cur >= 1e-6)
                    visited[Xi - X.begin()] = 1;

                gmm::add(gmm::scaled(*Xi, weight_cur), meanshift);
                denominator += weight_cur;
            }

            gmm::scale(meanshift, 1.0 / denominator);
            double distance_moved = gmm::vect_dist2(meanshift, mode, beta1, beta2, beta3);
            gmm::copy(meanshift, mode);

#ifdef DEBUG
            std::cout << "iter " << iter << ", moved "
                      << distance_moved << std::endl;
#endif
            iter += 1;
            if (distance_moved <= 1e-8)
                converged = true;
        } while (converged == false);

        return (iter);
    }

} // namespace clustering
