//
// Created by klaus on 20.01.19.
//

#include <unordered_set>

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

void LatticeBoltzmann2D::init()
{
   const size_t nvertices = mesh->getNumVertices();

   LOG_T(INFO) << "Initializing ..." << LogFlags::ENDL;
   ProgressBar progress;
   progress.start(nvertices);

   const vector<Vector2d>& u_0 = u->getValue(0);
   if (f->getSize() == 0)
      f->addValue(vector<Matrix<double, 9, 1>>(nvertices, Matrix<double, 9, 1>::Zero()));
   vector<Matrix<double, 9, 1>> &f_0 = f->getValue(0);
   if (rho->getSize() == 0)
      rho->addValue(vector<double>(nvertices, 1.0));
   vector<double>& rho_0 = rho->getValue(0);
   if (fe->getSize() == 0)
      fe->addValue(vector<Matrix<double, 9, 1>>(nvertices));
   vector<Matrix<double, 9, 1>>& fe_0 = fe->getValue(0);
   const double epsilon = 1.0e-6;
   size_t done = 0;
#ifdef _OPENMP
#pragma omp parallel for shared(progress), schedule(dynamic)
#endif
   for (size_t iP = 0; iP < nvertices; iP++) {
      const Vector2d &u = u_0[iP];
      const Matrix<double, 9, 1> cu = betac_i * u;
      const Matrix<double, 9, 1> cusqr = 0.5 * (cu.cwiseProduct(cu).array() - u.squaredNorm() / beta).matrix();
      double rho_pre = rho_0[iP];
      Matrix<double, 9, 1> &f = f_0[iP];
      double& rho = rho_0[iP];
      Matrix<double, 9, 1> &fe = fe_0[iP];
      /*do {
         rho_pre = rho;
         fe = rho * w_i.cwiseProduct(Matrix<double, 9, 1>::Identity() + cu + cusqr);
         const Matrix<double, 9, 1> fc = f + tau * (fe - f);
         f = fc + 1.0 / tau * (fe - fc)
         rho = fc.sum();
      } while (fabs(rho- rho_pre) < epsilon);*/
#pragma omp critical
      {
         done++;
         progress.update(done);
      }
   }
   progress.stop();
}

void LatticeBoltzmann2D::update()
{
   const size_t t = f->getSize() - 1;
   const size_t nvertices = mesh->getNumVertices();

   LOG_T(INFO) << "[Step " << t + 1 << "] Preparing ..." << LogFlags::ENDL;
   ProgressBar progress;
   progress.start(nvertices);

   const vector<Matrix<double, 9, 1>> &f_t = f->getValue(t);
   u->addValue(vector<Vector2d>(nvertices, Vector2d::Zero()));
   vector<Vector2d>& u_t = u->getValue(t);
   rho->addValue(vector<double>(nvertices));
   vector<double>& rho_t = rho->getValue(t);
   fe->addValue(vector<Matrix<double, 9, 1>>(nvertices));
   vector<Matrix<double, 9, 1>>& fe_t = fe->getValue(t);
#ifdef _OPENMP
#pragma omp parallel for shared(progress)
#endif
   for (size_t iP = 0; iP < nvertices; iP++) {
      const Matrix<double, 9, 1> &f = f_t[iP];
      const double rho = f.sum();
      rho_t[iP] = rho;
      u_t[iP] = c_i * f / rho;
      const Vector2d &u = u_t[iP];
      const Matrix<double, 9, 1> cu = betac_i * u;
      const Matrix<double, 9, 1> cusqr = 0.5 * (cu.cwiseProduct(cu).array() - u.squaredNorm() / beta).matrix();
      fe_t[iP] = rho * w_i.cwiseProduct(Matrix<double, 9, 1>::Identity() + cu + cusqr);
#pragma omp critical
      progress.update(iP);
   }
   progress.stop();
}

void LatticeBoltzmann2D::calculateNextTime(double dt)
{
   const size_t t = f->getSize() - 1;
   const size_t nvertices = mesh->getNumVertices();

   LOG_T(INFO) << "[Step " << t + 1 << "] Calculating ..." << LogFlags::ENDL;
   ProgressBar progress;
   progress.start(nvertices);

   const double dttau = dt / tau;
   const vector<Matrix<double, 9, 1>>& f_t = f->getValue(t);
   const vector<Matrix<double, 9, 1>>& fe_t = fe->getValue(t);
   f->addValue(vector<Matrix<double, 9, 1>>(nvertices));
   vector<Matrix<double, 9, 1>>& f_t2 = f->getValue(t + 1);

#ifdef _OPENMP
#pragma omp parallel for shared(progress)
#endif
   for (size_t iP = 0; iP < nvertices; iP++) {
      const auto boundaryit = infacing_boundary_c.find(iP);
      const bool isboundary = boundaryit != infacing_boundary_c.end();
      const auto &S_ik = S[iP];
      const vector<double>& C_k = C[iP];
      const Matrix<double, 9, 1>& f = f_t[iP];
      Matrix<double, 9, 1>& f2 = f_t2[iP];
      f2 = f;
      for (int iK = -1; iK < C_k.size(); iK++) {
         size_t iv = iP;
         Matrix<double, 9, 1> C = Matrix<double, 9, 1>::Zero();
         C(0) = 1.0 / 3.0;
         Matrix<double, 9, 1> S = Matrix<double, 9, 1>::Zero();
         if (iK >= 0) {
            iv = P_K[iP][iK];
            if (iK < 9)
               C(iK) = C_k[iK];
            S = S_ik.col(iK);
         }
         const Matrix<double, 9, 1> &fe = fe_t[iv];
         const Matrix<double, 9, 1> &f_K = f_t[iv];
         if (!isboundary)
            f2 += dt*S.cwiseProduct(f_K) - dttau*C.cwiseProduct(f_K - fe);
         else {
            const uint8_t infacings = boundaryit->second;
            for (uint8_t i = 0; i <= 8; i += 2) {
               if (i > 0) {
                  const uint8_t mask_odd  = (1 << (i - 2));
                  const uint8_t mask_even = (1 << (i - 1));
                  if ((infacings & mask_odd) && !(infacings & mask_even)) {
                     f2(i) += f(i) + dt*S(i) * f_K(i) - dttau*C(i) * (f_K(i) - fe(i));
                     f2(i - 1) = f2(i);
                  } else if (!(infacings & mask_odd) && (infacings & mask_even)) {
                     f2(i - 1) += f(i - 1) + dt*S(i - 1) * f_K(i - 1) - dttau*C(i - 1) * (f_K(i - 1) - fe(i - 1));
                     f2(i) = f2(i - 1);
                  } else {
                     const double fi1 = dt*S(i) * f_K(i) - dttau*C(i) * (f_K(i) - fe(i));
                     const double fi2 = dt*S(i - 1) * f_K(i - 1) - dttau*C(i - 1) * (f_K(i - 1) - fe(i - 1));
                     const double fi = 0.5 * (f(i) + fi1 + f(i - 1) + fi2);
                     f2(i) = fi;
                     f2(i - 1) = fi;
                  }
               } else
                  f2(i) += dt*S(i) * f_K(i) - dttau*C(i) * (f_K(i) - fe(i));
            }
         }
      }
#pragma omp critical
      progress.update(iP);
   }
   progress.stop();
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

LatticeBoltzmann2D::LatticeBoltzmann2D(Mesh<2, 2>* mesh, const vector<Mesh<2, 1>*>& boundaries, double tau, double c_s)
   : mesh(mesh), tau(tau), c_s(c_s), beta(1.0 / (c_s * c_s))
{
   auto E_ij = mesh->getOrCreateAttributeOnVertex<vector<Vector2d>>("E_ij");
   auto C_ij = mesh->getOrCreateAttributeOnVertex<vector<Matrix<double, 2, 1>>>("C_ij");
   auto C_i = mesh->getOrCreateAttributeOnFace<Vector2d>("C_i");
   auto E_i = mesh->getOrCreateAttributeOnEdge<Vector2d>("E_i");
   auto AttrN_K = mesh->getOrCreateAttributeOnVertex<vector<Vector2d>>("N_K");
   auto AttrP_K = mesh->getOrCreateAttributeOnVertex<vector<int>>("P_K");

   f = mesh->getOrCreateListAttributeOnVertex<vector<Matrix<double, 9, 1>>>("f");
   fe = mesh->getOrCreateListAttributeOnVertex<vector<Matrix<double, 9, 1>>>("fe");
   u = mesh->getOrCreateListAttributeOnVertex<vector<Vector2d>>("u");
   rho = mesh->getOrCreateListAttributeOnVertex<vector<double>>("rho");

   const double c = 1.0 / sqrt(2.0);
   c_i.col(0) <<  0.0,  0.0;
   c_i.col(1) <<  1.0,  0.0;
   c_i.col(2) <<  0.0,  1.0;
   c_i.col(3) << -1.0,  0.0;
   c_i.col(4) <<  0.0, -1.0;
   c_i.col(5) <<  c,  c;
   c_i.col(6) << -c,  c;
   c_i.col(7) << -c, -c;
   c_i.col(8) <<  c, -c;

   const double w1 = 1.0 / 9.0;
   const double w2 = 1.0 / 36.0;
   w_i << 4.0 / 9.0, w1, w1, w1, w1, w2, w2, w2, w2;
   betac_i = (c_i * beta).transpose();

   const size_t nvertices = mesh->getNumVertices();

   LOG_T(INFO) << "Calculate boundary conditions ..." << LogFlags::ENDL;
   unordered_map<size_t, Vector2d> vertex_normals;
   unordered_set<size_t> processed_edges;
   for (size_t ib = 0; ib < boundaries.size(); ib++) {
      Mesh<2, 1>* boundary = boundaries[ib];
      auto Nb = boundary->getOrCreateAttributeOnVertex<Vector2d>("Nb");
      auto infaceing = boundary->getOrCreateAttributeOnVertex<uint8_t>("infaceing");
      for (size_t ie = 0; ie < boundary->getNumEdges(); ie++) {
         Edge<2> *edge = boundary->getEdge(ie);
         if (processed_edges.find(edge->getID()) != processed_edges.end())
            continue;
         processed_edges.insert(edge->getID());
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
         for (size_t iv = 0; iv < 2; iv++) {
            const size_t vid = (*edge)[iv];
            const auto nit = vertex_normals.find(vid);
            if (nit != vertex_normals.end()) {
               Vector2d vnormal = n + nit->second;
               vnormal /= vnormal.norm();
               for (uint8_t i = 1; i < c_i.size(); i++)
                  if (vnormal.dot(c_i.col(i)) < 0)
                     infacing_boundary_c[vid] |= (1 << (i - 1));
               const size_t vidx = boundary->getVertexIdx(vid);
               Nb->setValue(vidx, vnormal);
               infaceing->setValue(vidx, infacing_boundary_c[vid]);
            } else
               vertex_normals[vid] = n;
         }
      }
   }

   LOG_T(INFO) << "Calculate neighbour relations ..." << LogFlags::ENDL;
   ProgressBar progress;
   progress.start(mesh->getNumFaces());
   P_K.resize(nvertices);
   vector<vector<Vector2d>> N_K(nvertices);
   vector<vector<double>> V_K(nvertices);
   vector<double> V_P(nvertices, 0.0);
   array<Vector2d, 3> E;
   array<Vector2d, 3> P;
   array<int, 3> Pidx;
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

   LOG_T(INFO) << "Calculate streaming and collision matrices ..." << LogFlags::ENDL;
   progress.start(nvertices);
   S.resize(nvertices);
   C.resize(nvertices);
   for (size_t iv = 0; iv < nvertices; iv++) {
      const int n = P_K[iv].size();
      Matrix<double, 9, Dynamic> S_ik(9, n);
      C[iv].resize(n);
      const double V = V_P[iv];
      for (int i = 0; i < 9; i++) {
         for (int k = 0; k < n; k++) {
            S_ik(i, k) =  c_i.col(i).dot(N_K[iv][k]) / V;
            if (i == 0)
               C[iv][k] = V_K[iv][k] / (3.0 * V);
         }
      }
      S[iv] = S_ik;
      progress.update(iv);
   }
   progress.stop();
}