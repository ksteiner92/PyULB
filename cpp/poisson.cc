//
// Created by klaus on 02.03.19.
//

#include <vector>
#include <array>

#include "poisson.h"
#include "eigen.h"
#include "logger.h"

using namespace std;
using namespace Eigen;

template<uint Dim, uint TopDim>
Poisson<Dim, TopDim>::Poisson(IMesh* mesh)
: mesh(mesh)
{
   LOG_T(INFO) << "[Poisson] Pre-processing ..." << LogFlags::ENDL;
   ProgressBar progress;
   const size_t nvertices = mesh->getNumVertices();
   K.resize(nvertices, nvertices);
   Mesh<Dim, 0>* mesh1D = dynamic_cast<Mesh<Dim, 0>*>(mesh);
   progress.start(mesh->getNumBodies());
   for (size_t ic = 0; ic < mesh->getNumBodies(); ic++) {
      MeshElement* body = mesh->getBody(ic);
      vector<Matrix<double, Dim, 1>> p(body->getNumVertices());
      for (uint8_t i = 0; i < body->getNumVertices(); i++)
         p[i] = mesh1D->getPoint((*body)[i]);
      array<double, 3> b;
      array<double, 3> c;
      for (uint8_t i = 0; i < 3; i++) {
         b[i] = p[(i + 1) % 3](1) - p[(i + 2) % 3](1);
         c[i] = p[(i + 2) % 3](0) - p[(i + 1) % 3](0);
      }
      const double Omega = 2 * fabs(p[1](0) * p[2](1) + p[0](0) * p[1](1) + p[0](1) * p[2](0) -
                                p[1](0) * p[0](1) - p[2](0) * p[1](1) - p[2](1) * p[0](0));
      for (uint8_t i = 0; i < 3; i++) {
         for (uint8_t j = i; j < 3; j++) {
            const double Kij = (b[j] * b[i] + c[j] * c[i]) / Omega;
            const size_t k = min((*body)[i], (*body)[j]);
            const size_t l = max((*body)[i], (*body)[j]);
            K.coeffRef(k, l) += Kij;
         }
      }
      progress.update(ic);
   }
   progress.stop();
}

template<uint Dim, uint TopDim>
double Poisson<Dim, TopDim>::orientation(const Matrix<double, Dim, 1>& p,
                          const Matrix<double, Dim, 1>& q,
                          const Matrix<double, Dim, 1>& i)
{
   const Matrix<double, Dim, 1> pq = q - p;
   const Matrix<double, Dim, 1> pi = i - p;
   const double o = pq(0) * pi(1) - pq(1) * pi(0);
   if (abs(o) > 1.0e-16)
      return o;
   return 0.0;
}

template<uint Dim, uint TopDim>
void Poisson<Dim, TopDim>::grad(const string& u, const string& w) const
{
   Mesh<Dim, 0>* mesh0D = dynamic_cast<Mesh<Dim, 0>*>(mesh);
   auto phi = mesh0D->template getAttributeOnVertex<double>(u);
   if (phi == nullptr) {
      stringstream ss;
      ss << "Attribute '" << u << "' does not exist on vertices";
      throw invalid_argument(ss.str());
   }
   //auto E = mesh1D->template getOrCreateAttributeOnVertex<Matrix<double, Dim, 1>>(w);
   auto E = reinterpret_cast<Attribute<Matrix<double, Dim, 1>>*>(mesh->getAttributeOnBody(w));
   if (E == nullptr) {
      stringstream ss;
      ss << "Attribute '" << w << "' does not exist on mesh bodies";
      throw invalid_argument(ss.str());
   }
   LOG_T(INFO) << "[Poisson] Calculating gradient..." << LogFlags::ENDL;
   ProgressBar progress;
   progress.start(mesh->getNumBodies());
   //vector<double> vvol(mesh->getNumVertices(), 0.0);
   //vector<Matrix<double, Dim, 1>> gradient(mesh1D->getNumVertices(), Matrix<double, Dim, 1>::Zero());
   for (size_t ib = 0; ib < mesh->getNumBodies(); ib++) {
      MeshElement* body = mesh->getBody(ib);
      const size_t nvertices = body->getNumVertices();
      vector<Matrix<double, Dim, 1>> p(nvertices);
      for (uint8_t i = 0; i < nvertices; i++)
         p[i] = mesh0D->getPoint((*body)[i]);
      const double Omega = 0.5 * (p[1](0) * p[2](1) + p[0](0) * p[1](1) + p[0](1) * p[2](0) -
                                p[1](0) * p[0](1) - p[2](0) * p[1](1) - p[2](1) * p[0](0));
      E->getValue(ib) = Matrix<double, Dim, 1>::Zero();
      for (int i = 0; i < nvertices; i++) {
         E->getValue(ib)(0) += (p[(i + 1) % 3](1) - p[(i + 2) % 3](1)) * phi->getValue((*body)[i]) / (2.0 * Omega);
         E->getValue(ib)(1) += (p[(i + 2) % 3](0) - p[(i + 1) % 3](0)) * phi->getValue((*body)[i]) / (2.0 * Omega);
      }
      progress.update(ib);
   }
   progress.stop();
}

template<uint Dim, uint TopDim>
void Poisson<Dim, TopDim>::solve(const string& phi_str,
                                 const string& rho_str,
                                 const vector<size_t>& tmp)
{
   LOG_T(INFO) << "[Poisson] Assembling ..." << LogFlags::ENDL;
   //auto rho_attr = mesh1D->template getOrCreateAttributeOnVertex<double>(rho_str);
   Mesh<Dim, 1>* mesh1D = dynamic_cast<Mesh<Dim, 1>*>(mesh);
   auto phi_attr = mesh1D->template getOrCreateAttributeOnVertex<double>(phi_str);
   vector<size_t> dirichlet(tmp.size());
   copy(tmp.begin(), tmp.end(), dirichlet.begin());
   sort(dirichlet.begin(), dirichlet.end());
   const size_t nvertices = mesh1D->getNumVertices();
   VectorXd F = VectorXd::Zero(nvertices);
   SparseMatrix<double, RowMajor> A = K;
   for (size_t i = 0; i < dirichlet.size(); i++) {
      const size_t dv = dirichlet[i];
      for (SparseMatrix<double, RowMajor>::InnerIterator it(A, dv); it; ++it)
         it.valueRef() = it.col() == dv ? 1.0 : 0.0;
      const double phi = phi_attr->getValue(dv);
      F(dv) = phi;
      const auto eov = mesh1D->getEdgesOfVertex(dv);
      for (auto it = eov.first; it != eov.second; it++) {
         Edge<Dim>* edge = it->second;
         const size_t neighbor = (*edge)[0] == dv ? (*edge)[1] : (*edge)[0];
         const size_t k = min(neighbor, dv);
         const size_t l = max(neighbor, dv);
         A.coeffRef(k, l) = 0.0;
         if (!binary_search(dirichlet.begin(), dirichlet.end(), neighbor))
            F(neighbor) -= K.coeff(k, l) * phi;
      }
   }
   LOG_T(INFO) << "[Poisson] Solving ..." << LogFlags::ENDL;
   //Matrix<double, Dynamic, Dynamic, RowMajor> res = A.toDense();
   //copy(res.data(), res.data() + res.rows() * res.cols(), back_inserter(M));
   //copy(F.data(), F.data() + F.rows() * F.cols(), back_inserter(N));
   ConjugateGradient<SparseMatrix<double, RowMajor>, Eigen::Upper> solver;
   VectorXd phi = solver.compute(A).solve(F);
   if(solver.info()!= Success)
      throw runtime_error("Could not solve Poisson system");
   for (size_t i = 0; i < nvertices; i++)
      phi_attr->setValue(i, phi(i));
}

template class Poisson<1, 1>;
template class Poisson<2, 1>;
template class Poisson<2, 2>;
template class Poisson<3, 1>;
template class Poisson<3, 2>;
template class Poisson<3, 3>;