//
// Created by klaus on 20.01.19.
//

#ifndef LBM_LBM_H
#define LBM_LBM_H

#include <array>
#include <unordered_set>

#include "mesh.h"
#include "utils.h"
#include "eigen.h"
#include "logger.h"

class LatticeBoltzmann2D
{
public:
   LatticeBoltzmann2D(Mesh<2, 2>* mesh, const std::vector<Mesh<2, 1>*>& boundaries, double tau, double c_s);

   void calculate(double dt, size_t nsteps);

private:
   static constexpr double epsilon = 1.0e-10;
   Mesh<2, 2>* mesh;
   Attribute<std::vector<std::array<double, 9>>>* f;
   Attribute<std::vector<std::array<double, 9>>>* fe;
   Attribute<std::vector<Eigen::Vector2d>>* u;
   Attribute<std::vector<double>>* rho;
   std::vector<Eigen::Matrix<double, 9, Eigen::Dynamic>> S;
   std::vector<std::vector<double>> C;
   std::vector<std::vector<int>> P_K;
   std::array<Eigen::Vector2d, 9> c_i;
   std::array<Eigen::Vector2d, 9> betac_i;
   std::array<double, 9> sqrbetac_i;
   std::array<double, 9> w_i;
   std::vector<std::vector<uint8_t>> infacing_boundary_c;
   std::unordered_map<std::size_t, uint8_t> known_boundary_f;
   double c_s;
   double beta;
   double tau;

   void calculateForTime(double dttau, std::size_t t, ProgressBar &progress);

   void preparingForTime(std::size_t t, ProgressBar &progress);

   double getEquilibriumF_i(std::size_t i, std::size_t t) const;

   static double calcPivot(const Eigen::Vector2d &P,
                           const Eigen::Vector2d &P1,
                           const Eigen::Vector2d &P2,
                           const Eigen::Vector2d &C,
                           Eigen::Vector2d &N1,
                           Eigen::Vector2d &N2);

   inline static double triangleArea(const Eigen::Vector2d &A,
                                     const Eigen::Vector2d &B,
                                     const Eigen::Vector2d &C);

   static double orientation(const Eigen::Vector2d &p,
                             const Eigen::Vector2d &q,
                             const Eigen::Vector2d &i);

   template<class T, uint SimplexDim>
   static T interpolateVertexAttributeOnSimplex(Simplex<2, SimplexDim>* simplex, Attribute<T>* attr);
};


#endif //LBM_LBM_H
