//
// Created by klaus on 06.01.19.
//

#ifndef LBM_POINT_H
#define LBM_POINT_H

#include <array>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <algorithm>

#include "attributes.h"
#include "eigen.h"
#include "utils.h"

template<uint Dim>
class Vertex;

template<uint Dim>
class Edge;

template<uint Dim>
class Face;

template<uint Dim>
class Grid;

template<uint Dim, uint TopDim>
class Mesh;

template<uint Dim, uint TopDim>
class IMesher;

class MeshElement
{
public:
   virtual std::size_t size() const = 0;
};

template<uint Dim, uint SimplexDim>
class Simplex : public MeshElement
{
   static_assert(SimplexDim <= Dim, "Simplex dimension has to be smaller or equal than dimension");
   static_assert(Dim <= 3, "Dimension can only be 1, 2 or 3");

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
   typedef enable_if_else_t<SimplexDim == 1, Vertex<Dim>*,
           enable_if_else_t<SimplexDim == 2, Edge<Dim>*,
           enable_if_else_t<SimplexDim == 3, Face<Dim>*,
           int > > >  FacetElementType;

   typedef enable_if_else_t<SimplexDim == 1, std::array<FacetElementType, 2>,
           enable_if_else_t<SimplexDim == 2, std::array<FacetElementType, 3>,
           enable_if_else_t<SimplexDim == 3, std::array<FacetElementType, 4>,
           std::size_t > > >  FacetsType;

   typedef std::array<std::size_t, SimplexDim + 1> FacetPointListType;

private:
   template<class T, uint I>
   struct GetPoints
   {
      static constexpr int get(
              FacetPointListType &facetPoints,
              const typename Simplex<Dim, I>::FacetsType &facets,
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
              const typename Simplex<Dim, 0>::FacetsType &facets,
              std::size_t idx)
      {
         if (!std::any_of(std::begin(facetPoints),
                 std::end(facetPoints), [&](std::size_t i) { return i == facets;})) {
            facetPoints[idx] = facets;
            return idx + 1;
         }
         return idx;
      }
   };

public:
   Simplex() {}

   Simplex(const FacetsType& facets, std::size_t id)
      : facets(facets), id(id), pts({0})
   {
      GetPoints<void, SimplexDim>::get(pts, facets, 0);
   }

   std::size_t getID() const
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
   std::size_t id;
   FacetsType facets;
   FacetPointListType pts;

};

template<uint Dim>
class Vertex : public Simplex<Dim, 0>
{
public:
   Vertex() {}

   Vertex(const typename Vertex::FacetsType& facets, std::size_t id)
      : Simplex<Dim, 0>(facets, id)
   {}
};

template<uint Dim>
class Edge : public Simplex<Dim, 1>
{
public:
   Edge() {}

   Edge(const typename Edge::FacetsType& facets, std::size_t id)
           : Simplex<Dim, 1>(facets, id)
   {}

};

template<uint Dim>
class Face : public Simplex<Dim, 2>
{
   static_assert(Dim >= 2, "Face can only exist on 2D or 3D grid");

public:
   Face() {}

   Face(const typename Face::FacetsType& facets, std::size_t id)
           : Simplex<Dim, 2>(facets, id)
   {}

};

class Cell : public Simplex<3, 3>
{
public:
   Cell() {}

   Cell(const typename Cell::FacetsType& facets, std::size_t id)
      : Simplex<3, 3>(facets, id)
   {}

};

template<uint Dim, uint TopDim=Dim>
class Mesh
{};

template<uint Dim>
class Mesh<Dim, 0>
{
   static_assert(Dim >= 0 && Dim <= 3, "Dimension not supported");

public:
   Mesh(std::vector<Eigen::Matrix<double, Dim, 1>>* points);

   template<uint TopDim>
   Mesh(Mesh<Dim, TopDim>* mesh) : points(mesh->points)
   {
      static_assert(TopDim >= 0, "Dimension mismatch");
      vertices = mesh->vertices;
   }

   const Eigen::Matrix<double, Dim, 1>& getPoint(std::size_t idx) const;

   std::vector<Eigen::Matrix<double, Dim, 1>>* getPoints() const;

   std::size_t getNumPoints() const;

   Vertex<Dim>* getVertex(std::size_t idx) const;

   Vertex<Dim>* getVertexByID(std::size_t id) const;

   Vertex<Dim>* getOrCreateVertexByID(std::size_t id);

   std::size_t getNumVertices() const;

   template<class T>
   Attribute<T>* getOrCreateAttributeOnVertex(const std::string &name)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, vertexAttrs, vertices->size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnVertex(const std::string &name, const T& def)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, vertexAttrs, vertices->size(), def);
   }

protected:
   template<class T>
   static inline Attribute<T>* getOrCreateAttribute(const std::string &name,
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
   static inline Attribute<T>* getOrCreateAttribute(const std::string &name,
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

   template<class T>
   typename std::vector<T>::iterator insertSorted(std::vector<T>& vec, const T& item)
   {
      return vec.insert(std::upper_bound(vec.begin(), vec.end(), item), item);
   }

protected:
   std::vector<Eigen::Matrix<double, Dim, 1>>* points;
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> vertexAttrs;
   std::unique_ptr<std::vector<std::unique_ptr<Vertex<Dim>>>> vertices_owner;
   std::vector<std::unique_ptr<Vertex<Dim>>>* vertices;
   std::vector<std::size_t> refvertices;

};

template<uint Dim>
class Mesh<Dim, 1> : public Mesh<Dim, 0>
{
   static_assert(Dim >= 1, "Topological dimension not supported");

public:
   Mesh(std::vector<Eigen::Matrix<double, Dim, 1>>* points);

   template<uint TopDim>
   Mesh(Mesh<Dim, TopDim>* mesh) : Mesh<Dim, 0>(mesh)
   {
      static_assert(TopDim >= 1, "Dimension mismatch");
      edges = mesh->edges;
      vhash2id = mesh->vhash2id;
      hull = std::make_unique<Mesh<Dim, 0>>(this);
   }

   Edge<Dim>* getEdge(std::size_t idx) const;

   Edge<Dim>* getEdgeByID(std::size_t id) const;

   Edge<Dim>* getOrCreateEdgeByID(std::size_t id);

   Edge<Dim>* getEdge(std::size_t vid1, std::size_t vid2) const;

   Edge<Dim>* getOrCreateEdge(std::size_t vid1, std::size_t vid2);

   std::size_t getNumEdges() const;

   Mesh<Dim, 0>* getHull() const;

   template<class T>
   Attribute<T>* getOrCreateAttributeOnEdge(const std::string &name)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, edgeAttrs, edges->size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnEdge(const std::string &name, const T& def)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, edgeAttrs, edges->size(), def);
   }

protected:
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> edgeAttrs;
   std::unique_ptr<std::vector<std::unique_ptr<Edge<Dim>>>> edges_owner;
   std::unique_ptr<std::unordered_map<ullong, std::size_t>> vhash2id_owner;
   std::vector<std::unique_ptr<Edge<Dim>>>* edges;
   std::unordered_map<ullong, std::size_t>* vhash2id;
   std::vector<std::size_t> refedges;
   std::unordered_multimap<std::size_t, Edge<Dim>*> vertex2edge;
   LineHash linhash;
   std::unique_ptr<Mesh<Dim, 0>> hull;

};

template<uint Dim>
class Mesh<Dim, 2> : public Mesh<Dim, 1>
{
   static_assert(Dim >= 2, "Topological dimension not supported");

public:
   Mesh(std::vector<Eigen::Matrix<double, Dim, 1>>* points);

   template<uint TopDim>
   Mesh(Mesh<Dim, TopDim>* mesh) : Mesh<Dim, 1>(mesh)
   {
      static_assert(TopDim >= 2, "Dimension mismatch");
      faces = mesh->faces;
      ehash2id = mesh->ehash2id;
      hull = std::make_unique<Mesh<Dim, 1>>(this);
   }

   Face<Dim>* getFace(std::size_t idx) const;

   Face<Dim>* getFaceByID(std::size_t id) const;

   Face<Dim>* getFace(std::size_t eid1, std::size_t eid2, std::size_t eid3) const;

   Face<Dim>* getOrCreateFace(std::size_t eid1, std::size_t eid2, std::size_t eid3);

   Face<Dim>* getOrCreateFace(Edge<Dim>* e1, Edge<Dim>* e2, Edge<Dim>* e3);

   std::size_t getNumFaces() const;

   Mesh<Dim, 1>* getHull() const;

   template<class T>
   Attribute<T>* getOrCreateAttributeOnFace(const std::string &name)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, facesAttrs, faces->size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnFace(const std::string &name, const T& def)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, facesAttrs, faces->size(), def);
   }

protected:
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> facesAttrs;
   std::unique_ptr<std::vector<std::unique_ptr<Face<Dim>>>> faces_owner;
   std::unique_ptr<std::unordered_map<ullong, std::size_t>> ehash2id_owner;
   std::vector<std::unique_ptr<Face<Dim>>>* faces;
   std::vector<std::size_t> reffaces;
   std::unordered_map<ullong, std::size_t>* ehash2id;
   /*std::unordered_multimap<ullong, Face<Dim>*> point2face;
   std::unordered_multimap<ullong, Face<Dim>*> edge2face;*/
   TriangleHash trihash;
   std::unique_ptr<Mesh<Dim, 1>> hull;

};

template<>
class Mesh<3, 3> : public Mesh<3, 2>
{
public:
   Mesh(std::vector<Eigen::Matrix<double, 3, 1>>* points);

   Cell* getCell(std::size_t idx) const;

   std::size_t getNumCells() const;

   template<class T>
   Attribute<T>* getOrCreateAttributeOnBody(const std::string &name)
   {
      return Mesh<3, 0>::template getOrCreateAttribute<T>(name, bodiesAttrs, cells.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnBody(const std::string &name, const T& def)
   {
      return Mesh<3, 0>::template getOrCreateAttribute<T>(name, bodiesAttrs, cells.size(), def);
   }

protected:
   std::vector<std::unique_ptr<Cell>> cells;
   std::unordered_multimap<unsigned long long int, Cell*> point2body;
   std::unordered_multimap<unsigned long long int, Cell*> edge2body;
   std::unordered_multimap<unsigned long long int, Cell*> face2body;
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> bodiesAttrs;
   Mesh<3, 2>* hull;

};

template<uint Dim, uint TopDim>
class HullMesh : public Mesh<Dim, TopDim - 1>
{
public:
   HullMesh(Mesh<Dim, TopDim>& mesh);



};

template<uint Dim, uint TopDim>
class IMesher
{
public:
   virtual void generate(Mesh<Dim, TopDim>& mesh) = 0;
};

#endif //LBM_POINT_H
