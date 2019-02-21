//
// Created by klaus on 20.01.19.
//

#include "eigen.h"
#include "lbm.h"
#include "logger.h"

using namespace std;
using namespace Eigen;

double LatticeBoltzmann2D::triangleArea(const Vector2d& A, const Vector2d& B, const Vector2d& C)
{
   return abs(A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1])) * 0.5;
}

template<class T, uint SimplexDim>
T LatticeBoltzmann2D::interpolateVertexAttributeOnSimplex(Simplex<2, SimplexDim>* simplex,
        Attribute<T>* attr)
{
   T val = 0;
   for (int i = 0; i < (SimplexDim + 1); i++)
      val += attr->getValue((*simplex)[i]);
   return val / ((T) (SimplexDim + 1));
}

double LatticeBoltzmann2D::getEquilibriumF_i(size_t i, size_t t) const
{
   double rho = 0.0;
}

void LatticeBoltzmann2D::preparingForTime(size_t t, ProgressBar &progress)
{
   LOG_T(INFO) << "[Step " << t + 1 << "] Preparing ..." << LogFlags::ENDL;
   progress.start(mesh->getNumVertices());
#ifdef _OPENMP
#pragma omp parallel for shared(progress)
#endif
   for (size_t iP = 0; iP < mesh->getNumVertices(); iP++) {
      vector<Vector2d> &u_P = u->getValue(iP);
      vector<double> &rho_P = rho->getValue(iP);
      vector<array<double, 9>> &fe_P = fe->getValue(iP);
      const vector<array<double, 9>> &f_P = f->getValue(iP);
      for (uint8_t i = 0; i < 9; i++) {
         rho_P[t] += f_P[t][i];
         u_P[t] += c_i[i] * f_P[t][i];
      }
      u_P[t] /= rho_P[t];
      for (uint8_t i = 0; i < 9; i++)
         fe_P[t][i] = rho_P[t] * w_i[i] * (1.0 + betac_i[i].dot(u_P[t]) +
                                           0.5 * (sqrbetac_i[i] - 1.0) / (u_P[t].squaredNorm()));
      progress.update(iP);
   }
   progress.stop();
}

void LatticeBoltzmann2D::calculateForTime(double dt, size_t t, ProgressBar &progress)
{
   preparingForTime(t, progress);
   LOG_T(INFO) << "[Step " << t + 1 << "] Calculating ..." << LogFlags::ENDL;
   const double dttau = dt / tau;
   progress.start(mesh->getNumVertices());
#ifdef _OPENMP
#pragma omp parallel for shared(progress)
#endif
   for (size_t iP = 0; iP < mesh->getNumVertices(); iP++) {
      const auto &S_ik = S[iP];
      const vector<double> &C_k = C[iP];
      const vector<double> &rho_P = rho->getValue(iP);
      const vector<Vector2d> &u_P = u->getValue(iP);
      const vector<array<double, 9>> &fe_P = fe->getValue(iP);
      vector<array<double, 9>> &f_P = f->getValue(iP);
      for (uint8_t i = 0; i < 9; i++)
         f_P[t + 1][i] -= dttau / 3.0 * (f_P[t][i] - fe_P[t][i]);
      for (int iK = 0; iK < C_k.size(); iK++) {
         const vector<array<double, 9>> &f_PK = f->getValue(P_K[iP][iK]);
         const vector<array<double, 9>> &fe_PK = fe->getValue(P_K[iP][iK]);
         const vector<Vector2d> &u_PK = u->getValue(P_K[iP][iK]);
         const vector<double> &rho_PK = rho->getValue(P_K[iP][iK]);
         for (uint8_t i = 0; i < 9; i++) {
            f_P[t + 1][i] += dt * S_ik(i, iK) * f_PK[t][i];
            if (i == iK)
               f_P[t + 1][i] -= dttau * C_k[iK] * (f_PK[t][i] - fe_PK[t][i]);
         }
      }
      progress.update(iP);
   }
   progress.stop();
}

void LatticeBoltzmann2D::calculate(double dt, size_t nsteps)
{
   ProgressBar progress;
   for (size_t iP = 0; iP < mesh->getNumVertices(); iP++) {
      f->getValue(iP).resize(nsteps, {1.0});
      fe->getValue(iP).resize(nsteps);
      u->getValue(iP).resize(nsteps, Vector2d::Zero());
      rho->getValue(iP).resize(nsteps, 0.0);
   }
   for (size_t t = 0; t < nsteps - 1; t++)
      calculateForTime(dt, t, progress);
   preparingForTime(nsteps - 1, progress);
}

double LatticeBoltzmann2D::calcPivot(const Vector2d &P,
                                    const Vector2d &E1,
                                    const Vector2d &E2,
                                    const Vector2d &C,
                                    Vector2d &N1,
                                    Vector2d &N2)
{
   const Vector2d CE1 = E1 - C;
   const Vector2d CE2 = E2 - C;
   N1 << CE1(1), CE1(0);
   N2 << CE2(1), CE2(0);
   if (((E1(0) - P(0)) * (C(1) - P(1)) - (E1(1) - P(1)) * (C(0) - P(0))) > 0) {
      N1(0) = -N1(0);
      N2(1) = -N2(1);
   } else {
      N1(1) = -N1(1);
      N2(0) = -N2(0);
   }
   return triangleArea(P, E1, C) + triangleArea(P, E2, C);
}

double LatticeBoltzmann2D::orientation(const Vector2d &p,
                                       const Vector2d &q,
                                       const Vector2d &i)
{
   const Vector2d pq = q - p;
   const Vector2d pi = i - p;
   const double o = pq(0) * pi(1) - pq(1) * pi(0);
   if (abs(o) > LatticeBoltzmann2D::epsilon)
      return o;
   return 0.0;
}

LatticeBoltzmann2D::LatticeBoltzmann2D(Mesh<2, 2>* mesh, const vector<Mesh<2, 1>*>& boundaries, double tau, double c_s)
   : mesh(mesh), tau(tau), c_s(c_s), beta(1.0 / (c_s * c_s))
{
   auto E_ij = mesh->getOrCreateAttributeOnVertex<vector<Vector2d>>("E_ij");
   auto C_ij = mesh->getOrCreateAttributeOnVertex<vector<Matrix<double, 2, 1>>>("C_ij");
   auto C_i = mesh->getOrCreateAttributeOnFace<Vector2d>("C_i");
   auto E_i = mesh->getOrCreateAttributeOnEdge<Vector2d>("E_i");
   auto AttrN_K = mesh->getOrCreateAttributeOnVertex<vector<Vector2d>>("N_K");
   auto AttrP_K = mesh->getOrCreateAttributeOnVertex<vector<int>>("P_K");
   f = mesh->getOrCreateAttributeOnVertex<vector<array<double, 9>>>("f");
   fe = mesh->getOrCreateAttributeOnVertex<vector<array<double, 9>>>("fe");
   u = mesh->getOrCreateAttributeOnVertex<vector<Vector2d>>("u");
   rho = mesh->getOrCreateAttributeOnVertex<vector<double>>("rho");

   P_K.resize(mesh->getNumVertices());
   vector<vector<Vector2d>> N_K(mesh->getNumVertices());
   vector<vector<double>> V_K(mesh->getNumVertices());
   vector<double> V_P(mesh->getNumVertices(), 0.0);
   array<Vector2d, 3> E;
   array<Vector2d, 3> P;
   array<int, 3> Pidx;
   LOG_T(INFO) << "Calculate neighbour relations ..." << LogFlags::ENDL;
   ProgressBar progress;
   progress.start(mesh->getNumFaces());

   infacing_boundary_c.resize(boundaries.size());
   for (size_t ib = 0; ib < boundaries.size(); ib++) {
      Mesh<2, 1>* boundary = boundaries[ib];
      auto Nb = boundary->getOrCreateAttributeOnEdge<Vector2d>("Nb");
      infacing_boundary_c[ib].resize(boundary->getNumEdges(), 0);
      for (size_t ie = 0; ie < boundary->getNumEdges(); ie++) {
         Edge<2> *edge = boundary->getEdge(ie);
         Face<2>* face = mesh->getFacesOfEdge(edge).first->second;
         const Vector2d& B1 = boundary->getPoint((*edge)[0]);
         const Vector2d& B2 = boundary->getPoint((*edge)[1]);
         Vector2d P;
         for (size_t iv = 0; iv < face->getNumVertices(); iv++) {
            if (((*face)[iv] != (*edge)[0]) && ((*face)[iv] != (*edge)[1])) {
               P = mesh->getPoint((*face)[iv]);
               break;
            }
         }
         const Vector2d e = (B2 - B1);
         const Vector2d E = B1 + e * 0.5;
         const Vector2d C = P + (E - P) * 2.0 / 3.0;
         const Vector2d t = C - E;
         Vector2d n(e[1], -e[0]);
         if (n.dot(t) > 0.0)
            n *= -1.0;
         for (uint16_t i = 1; i < c_i.size(); i++)
            infacing_boundary_c[ib][ie] |= (1 << (i - 1));
         Nb->setValue(ie, n);
      }
   }
   for (size_t i = 0; i < mesh->getNumFaces(); i++) {
      // Calculate face center point
      const Face<2>* face = mesh->getFace(i);
      Pidx[0] = (*face)[0];
      P[0] = mesh->getPoint(Pidx[0]);
      Pidx[1] = (*face)[1];
      P[1] = mesh->getPoint(Pidx[1]);
      Pidx[2] = (*face)[2];
      P[2] = mesh->getPoint(Pidx[2]);
      E[0] = P[0] + (P[1] - P[0]) * 0.5;
      const Vector2d C = P[2] + (E[0] - P[2]) * 2.0 / 3.0;
      C_i->setValue(i, C);
      C_ij->getValue(Pidx[0]).push_back(C);
      C_ij->getValue(Pidx[1]).push_back(C);
      C_ij->getValue(Pidx[2]).push_back(C);

      // Calculate edge center points
      E_i->setValue(mesh->getEdge(Pidx[0], Pidx[1])->getID(), E[0]);
      E_ij->getValue(Pidx[0]).push_back(E[0]);
      E_ij->getValue(Pidx[1]).push_back(E[0]);
      E[1] = P[0] + (P[2] - P[0]) * 0.5;
      E_i->setValue(mesh->getEdge(Pidx[0], Pidx[2])->getID(), E[1]);
      E_ij->getValue(Pidx[0]).push_back(E[1]);
      E_ij->getValue(Pidx[2]).push_back(E[1]);
      E[2] = P[1] + (P[2] - P[1]) * 0.5;
      E_i->setValue(mesh->getEdge(Pidx[1], Pidx[2])->getID(), E[2]);
      E_ij->getValue(Pidx[1]).push_back(E[2]);
      E_ij->getValue(Pidx[2]).push_back(E[2]);

      const array<int, 3> iE1 = {0, 0, 2};
      const array<int, 3> iE2 = {1, 2, 1};
      const array<int, 6> iP_K = {1, 2, 0, 2, 0, 1};
      Vector2d N1, N2;
      for (int iP = 0; iP < 3; iP++) {
         const double V = calcPivot(P[iP], E[iE1[iP]], E[iE2[iP]], C, N1, N2);
         V_P[Pidx[iP]] += V;
         Vector2d N = 5.0 / 12.0 * N1 + 2.0 / 12.0 * N2;
         for (int j = 0; j < 2; j++) {
            const auto it = find(P_K[Pidx[iP]].begin(), P_K[Pidx[iP]].end(), Pidx[iP_K[2 * iP + j]]);
            if (it == P_K[Pidx[iP]].end()) {
               P_K[Pidx[iP]].push_back(Pidx[iP_K[2 * iP  + j]]);
               AttrP_K->getValue(Pidx[iP]).push_back(Pidx[iP_K[2 * iP  + j]]);
               N_K[Pidx[iP]].push_back(N);
               AttrN_K->getValue(Pidx[iP]).push_back(N);
               V_K[Pidx[iP]].push_back(V);
            } else {
               const size_t idx = distance(P_K[Pidx[iP]].begin(), it);
               N_K[Pidx[iP]][idx] += N;
               AttrN_K->getValue(Pidx[iP])[idx] = N_K[Pidx[iP]][idx];
               V_K[Pidx[iP]][idx] += V;
            }
            N = 5.0 / 12.0 * N2 + 2.0 / 12.0 * N1;
         }
      }
      progress.update(i);
   }
   progress.stop();

   const double c = 1.0 / sqrt(2.0);
   c_i[0] << 0.0, 0.0;
   c_i[1] << 1.0, 0.0;
   c_i[2] << 0.0, 1.0;
   c_i[3] << -1.0, 0.0;
   c_i[4] << 0.0, -1.0;
   c_i[5] << c, c;
   c_i[6] << -c, c;
   c_i[7] << -c, -c;
   c_i[8] << c, -c;
   const double w1 = 1.0 / 9.0;
   const double w2 = 1.0 / 36.0;
   w_i = {4.0 / 9.0, w1, w1, w1, w1, w2, w2, w2, w2};
   for (size_t i = 0; i < 9; i++) {
      betac_i[i] = c_i[i] / c_s;
      sqrbetac_i[i] = betac_i[i].squaredNorm();
   }

   LOG_T(INFO) << "Calculate streaming and collision matrices ..." << LogFlags::ENDL;
   progress.start(mesh->getNumVertices());
   S.resize(mesh->getNumVertices());
   C.resize(mesh->getNumVertices());
   for (size_t iv = 0; iv < mesh->getNumVertices(); iv++) {
      const int n = P_K[iv].size();
      Matrix<double, 9, Dynamic> S_ik(9, n);
      C[iv].resize(n);
      const double V = V_P[iv];
      for (int i = 0; i < 9; i++) {
         for (int k = 0; k < n; k++) {
            S_ik(i, k) =  c_i[i].dot(N_K[iv][k]) / V;
            if (i == 0)
               C[iv][k] = V_K[iv][k] / (3.0 * V);
         }
      }
      S[iv] = S_ik;
      progress.update(iv);
   }
   progress.stop();
}