//
// Created by klaus on 06.01.19.
//

#ifndef LBM_POINT_H
#define LBM_POINT_H

#include <array>
#include <unordered_map>
#include <memory>
#include <algorithm>

#include "attributes.h"
#include "eigen.h"

template<uint Dim, uint TopDim> class Vertex;
template<uint Dim, uint TopDim> class Edge;
template<uint Dim, uint TopDim> class Face;
template<uint Dim, uint TopDim> class Mesh;

class MeshElement
{
public:
   virtual std::size_t size() const = 0;
};

template<uint SimplexDim, uint Dim, uint TopDim>
class Simplex : public MeshElement
{
   static_assert(SimplexDim <= Dim, "Simplex dimension has to be smaller or equal than dimension");
   static_assert(Dim <= 3, "Dimension can only be 1, 2 or 3");
   static_assert(TopDim <= Dim, "Topological dimension must be smaller or equal than dimension");

private:
   template<bool, class TrueT = void, class FalseT = void>
   struct enable_if_else {};

   template<class TrueT, class FalseT>
   struct enable_if_else<true, TrueT, FalseT>
   {
      typedef TrueT type;
   };

   template<class TrueT, class FalseT>
   struct enable_if_else<false, TrueT, FalseT>
   {
      typedef FalseT type;
   };

   template< bool B, class TrueT = void, class FalseT = void>
   using enable_if_else_t = typename enable_if_else<B, TrueT, FalseT>::type;

public:
   typedef enable_if_else_t<SimplexDim == 1, Vertex<Dim, TopDim>*,
           enable_if_else_t<SimplexDim == 2, Edge<Dim, TopDim>*,
           enable_if_else_t<SimplexDim == 3, Face<Dim, TopDim>*,
           int > > >  FacetElementType;

   typedef enable_if_else_t<SimplexDim == 1, std::array<FacetElementType, 2>,
           enable_if_else_t<SimplexDim == 2, std::array<FacetElementType, 3>,
           enable_if_else_t<SimplexDim == 3, std::array<FacetElementType, 4>,
           int > > >  FacetsType;

   typedef std::array<int, SimplexDim + 1> FacetPointListType;

private:
   template<class T, uint I>
   struct GetPoints
   {
      static constexpr int get(
              FacetPointListType &facetPoints,
              const typename Simplex<I, Dim, TopDim>::FacetsType &facets,
              int idx)
      {
         for (uint i = 0; i < facets.size(); i++)
            idx = GetPoints<T, I - 1>::get(facetPoints, facets[i]->getFacets(), idx);
         return idx;
      }
   };

   template<class T>
   struct GetPoints<T, 0>
   {
      static constexpr int get(
              FacetPointListType &facetPoints,
              const typename Simplex<0, Dim, TopDim>::FacetsType &facets,
              int idx)
      {
         if (!std::any_of(std::begin(facetPoints),
                 std::end(facetPoints), [&](int i) { return i == facets;})) {
            facetPoints[idx] = facets;
            return idx + 1;
         }
         return idx;
      }
   };

public:
   Simplex() {}

   Simplex(Mesh<Dim, TopDim>* mesh, const FacetsType& facets, int id)
      : mesh(mesh), facets(facets), id(id), pts({-1})
   {
      GetPoints<void, SimplexDim>::get(pts, facets, 0);
   }

   int getID() const
   {
      return id;
   }

   int operator[](std::size_t idx) const
   {
      if (idx >= SimplexDim + 1) {
         std::cout << "Simplex operator[]: Index out of bounds" << std::endl;
         std::cout.flush();
         throw std::out_of_range("Index out of range");
      }
      return pts[idx];
   }

   const Eigen::Matrix<double, Dim, 1>& getPoint(std::size_t idx) const
   {
      if (idx >= SimplexDim + 1) {
         std::cout << "Simplex getPoint: Index out of bounds" << std::endl;
         std::cout.flush();
         throw std::out_of_range("Index out of range");
      }
      return (*mesh->points)[pts[idx]];
   }

   FacetElementType getFacet(std::size_t idx)
   {
      if (idx >= facets.size()) {
         std::cout << "Simplex getFacet: Index out of bounds" << std::endl;
         std::cout.flush();
         throw std::out_of_range("Index out of range");
      }
      return facets[idx];
   }

   FacetsType& getFacets()
   {
      return facets;
   }

   std::size_t size() const override
   {
      return 0;
   }

protected:
   int id;
   FacetsType facets;
   FacetPointListType pts;
   Eigen::Matrix<double, Dim, 1> center;
   Mesh<Dim, TopDim>* mesh;

};

template<uint Dim, uint TopDim = Dim>
class Vertex : public Simplex<0, Dim, TopDim>
{
public:
   Vertex() {}

   Vertex(Mesh<Dim, TopDim>* mesh, const typename Vertex::FacetsType& hull, int id)
      : Simplex<0, Dim, TopDim>(mesh, hull, id)
   {
   }
};

template<uint Dim, uint TopDim = Dim>
class Edge : public Simplex<1, Dim, TopDim>
{
public:
   Edge() {}

   Edge(Mesh<Dim, TopDim>* mesh, const typename Edge::FacetsType& hull, int id)
           : Simplex<1, Dim, TopDim>(mesh, hull, id)
   {
   }

};

template<uint Dim, uint TopDim = Dim>
class Face : public Simplex<2, Dim, TopDim>
{
   static_assert(Dim >= 2, "Face can only exist on 2D or 3D grid");

public:
   Face() {}

   Face(Mesh<Dim, TopDim>* mesh, const typename Face::FacetsType& hull, int id)
           : Simplex<2, Dim, TopDim>(mesh, hull, id)
   {
   }

};

template<uint Dim>
class Cell : public Simplex<3, Dim, Dim>
{
   static_assert(Dim == 3, "Cell can only exist on a 3D grid");

public:
   Cell() {}

   Cell(Mesh<Dim, Dim>* mesh, const typename Cell::FacetsType& hull, int id)
      : Simplex<3, Dim, Dim>(mesh, hull, id)
   {
   }

};

class IMesh
{
public:
   template<class T>
   static inline Attribute<T>* getOrCreateAttribute(
           const std::string &name,
           std::unordered_map<std::string, std::unique_ptr<BaseAttribute>>& attrlst,
           std::size_t size,
           const T& def)
   {
      const auto it = attrlst.find(name);
      if (it == attrlst.end()) {
         attrlst[name] = std::make_unique<Attribute<T>>(size, def);
         return static_cast<Attribute<T> *>(attrlst[name].get());
      }
      return static_cast<Attribute<T> *>(it->second.get());
   }

   template<class T>
   static inline Attribute<T>* getOrCreateAttribute(
           const std::string &name,
           std::unordered_map<std::string, std::unique_ptr<BaseAttribute>>& attrlst,
           std::size_t size)
   {
      const auto it = attrlst.find(name);
      if (it == attrlst.end()) {
         attrlst[name] = std::make_unique<Attribute<T>>(size);
         return static_cast<Attribute<T> *>(attrlst[name].get());
      }
      return static_cast<Attribute<T> *>(it->second.get());
   }
};

template<uint Dim, uint TopDim = Dim>
class Mesh : public IMesh
{
   static_assert(Dim > 0 && Dim <= 3, "Dimension not supported");
   static_assert(TopDim <= Dim, "Topological dimension not supported");
};

template<uint Dim>
class Mesh<Dim, 1>
{
public:
   Mesh(std::vector<Eigen::Matrix<double, Dim, 1>>* points) : points(points) {}

   Vertex<Dim, Dim>* getVertex(std::size_t idx) const
   {
      if (idx >= vertices.size()) {
         std::cout << "Mesh1D getVertex: Index out of bounds" << std::endl;
         std::cout.flush();
         throw std::out_of_range("Vertex index out of range");
      }
      return vertices[idx].get();
   }

   std::size_t getNumVertices() const
   {
      return vertices.size();
   }

   Edge<Dim, Dim>* getEdge(std::size_t idx) const
   {
      if (idx >= edges.size()) {
         std::cout << "Mesh1D getEdge: Index out of bounds" << std::endl;
         std::cout.flush();
         throw std::out_of_range("Edge index out of range");
      }
      return edges[idx].get();
   }

   std::size_t getNumEdges() const
   {
      return edges.size();
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnVertex(const std::string &name)
   {
      return IMesh::getOrCreateAttribute<T>(name, vertexAttrs, vertices.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnVertex(const std::string &name, const T& def)
   {
      return IMesh::getOrCreateAttribute<T>(name, vertexAttrs, vertices.size(), def);
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnEdge(const std::string &name)
   {
      return IMesh::getOrCreateAttribute<T>(name, edgeAttrs, edges.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnEdge(const std::string &name, const T& def)
   {
      return IMesh::getOrCreateAttribute<T>(name, edgeAttrs, edges.size(), def);
   }

   std::vector<Eigen::Matrix<double, Dim, 1>>* points;
   std::vector<std::unique_ptr<Vertex<Dim>>> vertices;
   std::vector<std::unique_ptr<Edge<Dim>>> edges;
   std::unordered_map<unsigned long long int, Edge<Dim>*> edge2edge;
   std::unordered_multimap<int, Edge<Dim>*> point2edge;

private:
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> vertexAttrs;
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> edgeAttrs;

};

template<uint Dim>
class Mesh<Dim, 2> : public Mesh<Dim, 1>
{
public:
   Mesh(std::vector<Eigen::Matrix<double, Dim, 1>>* points) : Mesh<Dim, 1>(points) {}

   Face<Dim, Dim>* getFace(std::size_t idx) const
   {
      if (idx >= faces.size()) {
         std::cout << "Mesh2D getFace: Index out of bounds" << std::endl;
         std::cout.flush();
         throw std::out_of_range("Face index out of range");
      }
      return faces[idx].get();
   }

   std::size_t getNumFaces() const
   {
      return faces.size();
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnFace(const std::string &name)
   {
      return IMesh::getOrCreateAttribute<T>(name, facesAttrs, faces.size());
   }

   /*template<class T>
   Attribute<T>* getOrCreateAttributeOnFace(const std::string &name, const T& def)
   {
      return IMesh::getOrCreateAttribute<T>(name, facesAttrs, faces.size(), def);
   }*/

   std::vector<std::shared_ptr<Face<Dim>>> faces;
   std::unordered_map<unsigned long long int, Face<Dim>*> face2face;
   std::unordered_multimap<unsigned long long int, Face<Dim>*> edge2face;
   std::unordered_multimap<int, Face<Dim>*> point2face;

private:
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> facesAttrs;

};

template<uint Dim>
class Mesh<Dim, 3> : public Mesh<Dim, 2>
{
public:
   Mesh(std::vector<Eigen::Matrix<double, Dim, 1>>* points) : Mesh<Dim, 2>(points) {}

   Cell<Dim>* getBody(std::size_t idx) const
   {
      if (idx >= bodies.size()) {
         std::cout << "Mesh3D getBody: Index out of bounds" << std::endl;
         std::cout.flush();
         throw std::out_of_range("Cell index out of range");
      }
      return bodies[idx].get();
   }

   std::size_t getNumBodies() const
   {
      return bodies.size();
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnBody(const std::string &name)
   {
      return IMesh::getOrCreateAttribute<T>(name, bodiesAttrs, bodies.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnBody(const std::string &name, const T& def)
   {
      return IMesh::getOrCreateAttribute<T>(name, bodiesAttrs, bodies.size(), def);
   }

   std::vector<std::unique_ptr<Cell<3>>> bodies;
   std::unordered_multimap<unsigned long long int, Cell<3>*> point2body;
   std::unordered_multimap<unsigned long long int, Cell<3>*> edge2body;
   std::unordered_multimap<unsigned long long int, Cell<3>*> face2body;

private:
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> bodiesAttrs;

};

template<uint Dim>
class MeshGenerator
{
public:
   virtual void generate(Mesh<Dim, Dim>& mesh) = 0;
};

#endif //LBM_POINT_H
