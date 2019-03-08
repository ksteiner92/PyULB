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
   virtual std::size_t getNumVertices() const = 0;

   virtual std::size_t getID() const = 0;

   virtual std::size_t operator[](std::size_t idx) const = 0;
};

template<uint Dim, int SimplexDim>
struct SimplexTraits
{};

template<uint Dim>
struct SimplexTraits<Dim, 0>
{
   typedef int FacetElementType;
   typedef std::size_t FacetsType;
};

template<uint Dim>
struct SimplexTraits<Dim, 1>
{
   typedef Vertex<Dim>* FacetElementType;
   typedef std::array<FacetElementType, 2> FacetsType;
};

template<uint Dim>
struct SimplexTraits<Dim, 2>
{
   typedef Edge<Dim>* FacetElementType;
   typedef std::array<FacetElementType, 3> FacetsType;
};

template<uint Dim>
struct SimplexTraits<Dim, 3>
{
   typedef Face<Dim>* FacetElementType;
   typedef std::array<FacetElementType, 4> FacetsType;
};

template<uint Dim, int SimplexDim>
class Simplex : public MeshElement
{
   static_assert(SimplexDim <= Dim, "Simplex dimension has to be smaller or equal than dimension");
   static_assert(Dim <= 3, "Dimension can only be 1, 2 or 3");


public:
   typedef typename SimplexTraits<Dim, SimplexDim>::FacetElementType FacetElementType;
   typedef typename SimplexTraits<Dim, SimplexDim>::FacetsType FacetsType;
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

   std::size_t getID() const override
   {
      return id;
   }

   std::size_t operator[](std::size_t idx) const override
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

   std::size_t getNumVertices() const override
   {
      return pts.size();
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

class IMesh
{
public:
   virtual std::size_t getNumVertices() const = 0;

   virtual MeshElement* getVertex(std::size_t idx) const = 0;

   virtual std::size_t getVertexID(std::size_t idx) const = 0;

   virtual std::size_t getVertexIdx(std::size_t id) const = 0;

   virtual MeshElement* getVertexByID(std::size_t id) const = 0;

   virtual std::size_t getNumBodies() const = 0;

   virtual MeshElement* getBody(std::size_t idx) const = 0;

   virtual MeshElement* getBodyByID(std::size_t id) const = 0;

   virtual std::size_t getNumFacets() const = 0;

   virtual MeshElement* getFacet(std::size_t idx) const = 0;

   virtual MeshElement* getFacetByID(std::size_t id) const = 0;

   virtual std::size_t getNumRidges() const = 0;

   virtual MeshElement* getRidge(std::size_t idx) const = 0;

   virtual MeshElement* getRidgeByID(std::size_t id) const = 0;

   virtual std::size_t getNumPeaks() const = 0;

   virtual MeshElement* getPeak(std::size_t idx) const = 0;

   virtual MeshElement* getPeakByID(std::size_t id) const = 0;

   virtual IMesh* getHull() const = 0;

};

template<uint Dim, uint TopDim=Dim>
class Mesh
{};

template<uint Dim>
class Mesh<Dim, 0> : public IMesh
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

   virtual std::size_t getNumBodies() const override;

   virtual MeshElement* getBody(std::size_t idx) const override;

   virtual MeshElement* getBodyByID(std::size_t id) const override;

   virtual std::size_t getNumFacets() const override;

   virtual MeshElement* getFacet(std::size_t idx) const override;

   virtual MeshElement* getFacetByID(std::size_t id) const override;

   virtual std::size_t getNumRidges() const override;

   virtual MeshElement* getRidge(std::size_t idx) const override;

   virtual MeshElement* getRidgeByID(std::size_t id) const override;

   virtual std::size_t getNumPeaks() const override;

   virtual MeshElement* getPeak(std::size_t idx) const override;

   virtual MeshElement* getPeakByID(std::size_t id) const override;

   const Eigen::Matrix<double, Dim, 1>& getPoint(std::size_t idx) const;

   std::vector<Eigen::Matrix<double, Dim, 1>>* getPoints() const;

   std::size_t getNumPoints() const;

   Vertex<Dim>* getVertex(std::size_t idx) const override;

   std::size_t getVertexID(std::size_t idx) const override;

   std::size_t getVertexIdx(std::size_t id) const override;

   Vertex<Dim>* getVertexByID(std::size_t id) const override;

   Vertex<Dim>* getOrCreateVertexByID(std::size_t id);

   std::size_t getNumVertices() const override;

   virtual IMesh* getHull() const override;

   template<class T>
   Attribute<T>* getOrCreateAttributeOnVertex(const std::string &name)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, vertexAttrs, refvertices.size());
   }

   template<class T>
   ListAttribute<T>* getOrCreateListAttributeOnVertex(const std::string &name)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T, ListAttribute>(name, vertexAttrs, refvertices.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnVertex(const std::string &name, const T& def)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, vertexAttrs, refvertices.size(), def);
   }

   template<class T>
   ListAttribute<T>* getOrCreateListAttributeOnVertex(const std::string &name, const T& def)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T, ListAttribute>(name, vertexAttrs, refvertices.size(), def);
   }

protected:
   template<typename T, template<typename> typename AttibuteType = Attribute>
   static inline AttibuteType<T> *getOrCreateAttribute(const std::string &name,
                                                       std::unordered_map<std::string,
                                                               std::unique_ptr<BaseAttribute>> &attrlst,
                                                       std::size_t size,
                                                       const T &def)
   {
      const auto it = attrlst.find(name);
      if (it == attrlst.end()) {
         attrlst[name] = std::make_unique<AttibuteType<T>>(size, def);
         return static_cast<AttibuteType<T> *>(attrlst[name].get());
      }
      return static_cast<AttibuteType<T> *>(it->second.get());
   }

   template<typename T, template<typename> typename AttibuteType = Attribute>
   static inline AttibuteType<T> *getOrCreateAttribute(const std::string &name,
                                                       std::unordered_map<std::string,
                                                               std::unique_ptr<BaseAttribute>> &attrlst,
                                                       std::size_t size)
   {
      const auto it = attrlst.find(name);
      if (it == attrlst.end()) {
         attrlst[name] = std::make_unique<Attribute<T>>(size);
         return static_cast<AttibuteType<T> *>(attrlst[name].get());
      }
      return static_cast<AttibuteType<T> *>(it->second.get());
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
   /**
    *
    * @param points
    */
   Mesh(std::vector<Eigen::Matrix<double, Dim, 1>>* points);

   /**
    *
    * @tparam TopDim The topological dimension of the parent mesh
    * @param mesh The parent mesh
    */
   template<uint TopDim>
   Mesh(Mesh<Dim, TopDim>* mesh) : Mesh<Dim, 0>(mesh)
   {
      static_assert(TopDim >= 1, "Dimension mismatch");
      edges = mesh->edges;
      vhash2id = mesh->vhash2id;
      vertex2edge = mesh->vertex2edge;
      hull = std::make_unique<Mesh<Dim, 0>>(this);
   }

   virtual std::size_t getNumBodies() const override;

   virtual MeshElement* getBody(std::size_t idx) const override;

   virtual MeshElement* getBodyByID(std::size_t id) const override;

   virtual std::size_t getNumFacets() const override;

   virtual MeshElement* getFacet(std::size_t idx) const override;

   virtual MeshElement* getFacetByID(std::size_t id) const override;

   virtual std::size_t getNumRidges() const override;

   virtual MeshElement* getRidge(std::size_t idx) const override;

   virtual MeshElement* getRidgeByID(std::size_t id) const override;

   virtual std::size_t getNumPeaks() const override;

   virtual MeshElement* getPeak(std::size_t idx) const override;

   virtual MeshElement* getPeakByID(std::size_t id) const override;

   /**
    * Returns the i-th edge of this mesh.
    *
    * @param idx The index
    * @return The pointer to the edge if found, otherwise
    * nullptr
    */
   Edge<Dim>* getEdge(std::size_t idx) const;

   /**
    * Returns the edge with the ID @param id
    *
    * @param id The ID of the edge
    * @return The pointer to the edge if found, otherwise
    * nullptr
    */
   Edge<Dim>* getEdgeByID(std::size_t id) const;

   /**
    *
    * @param id
    * @return
    */
   Edge<Dim>* getOrCreateEdgeByID(std::size_t id);

   Edge<Dim>* getEdge(std::size_t vid1, std::size_t vid2) const;

   Edge<Dim>* getOrCreateEdge(std::size_t vid1, std::size_t vid2);

   std::size_t getNumEdges() const;

   virtual IMesh* getHull() const override;

   std::pair<typename std::unordered_multimap<std::size_t, Edge<Dim>*>::const_iterator,
           typename std::unordered_multimap<std::size_t, Edge<Dim>*>::const_iterator> getEdgesOfVertex(std::size_t idx) const;

   template<class T>
   Attribute<T>* getOrCreateAttributeOnEdge(const std::string &name)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, edgeAttrs, refedges.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnEdge(const std::string &name, const T& def)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, edgeAttrs, refedges.size(), def);
   }

protected:
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> edgeAttrs;
   std::unique_ptr<std::vector<std::unique_ptr<Edge<Dim>>>> edges_owner;
   std::unique_ptr<std::unordered_map<ullong, std::size_t>> vhash2id_owner;
   std::vector<std::unique_ptr<Edge<Dim>>>* edges;
   std::unordered_map<ullong, std::size_t>* vhash2id;
   std::vector<std::size_t> refedges;
   std::unique_ptr<std::unordered_multimap<std::size_t, Edge<Dim>*>> vertex2edge_owner;
   std::unordered_multimap<std::size_t, Edge<Dim>*>* vertex2edge;
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
      vertex2face = mesh->vertex2face;
      edge2face = mesh->edge2face;
      hull = std::make_unique<Mesh<Dim, 1>>(this);
   }

   virtual std::size_t getNumBodies() const override;

   virtual MeshElement* getBody(std::size_t idx) const override;

   virtual MeshElement* getBodyByID(std::size_t id) const override;

   virtual std::size_t getNumFacets() const override;

   virtual MeshElement* getFacet(std::size_t idx) const override;

   virtual MeshElement* getFacetByID(std::size_t id) const override;

   virtual std::size_t getNumRidges() const override;

   virtual MeshElement* getRidge(std::size_t idx) const override;

   virtual MeshElement* getRidgeByID(std::size_t id) const override;

   virtual std::size_t getNumPeaks() const override;

   virtual MeshElement* getPeak(std::size_t idx) const override;

   virtual MeshElement* getPeakByID(std::size_t id) const override;

   Face<Dim>* getFace(std::size_t idx) const;

   Face<Dim>* getFaceByID(std::size_t id) const;

   Face<Dim>* getFace(std::size_t eid1, std::size_t eid2, std::size_t eid3) const;

   Face<Dim>* getOrCreateFace(std::size_t eid1, std::size_t eid2, std::size_t eid3);

   Face<Dim>* getOrCreateFace(Edge<Dim>* e1, Edge<Dim>* e2, Edge<Dim>* e3);

   std::size_t getNumFaces() const;

   std::pair<typename std::unordered_multimap<std::size_t, Face<Dim>*>::const_iterator,
   typename std::unordered_multimap<std::size_t, Face<Dim>*>::const_iterator> getFacesOfEdge(Edge<Dim>* edge) const;

   virtual IMesh* getHull() const override;

   template<class T>
   Attribute<T>* getOrCreateAttributeOnFace(const std::string &name)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, facesAttrs, reffaces.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnFace(const std::string &name, const T& def)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, facesAttrs, reffaces.size(), def);
   }

protected:
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> facesAttrs;
   std::unique_ptr<std::vector<std::unique_ptr<Face<Dim>>>> faces_owner;
   std::unique_ptr<std::unordered_map<ullong, std::size_t>> ehash2id_owner;
   std::vector<std::unique_ptr<Face<Dim>>>* faces;
   std::vector<std::size_t> reffaces;
   std::unordered_map<ullong, std::size_t>* ehash2id;
   std::unique_ptr<std::unordered_multimap<std::size_t, Face<Dim>*>> vertex2face_owner;
   std::unordered_multimap<std::size_t, Face<Dim>*>* vertex2face;
   std::unique_ptr<std::unordered_multimap<std::size_t, Face<Dim>*>> edge2face_owner;
   std::unordered_multimap<std::size_t, Face<Dim>*>* edge2face;
   TriangleHash trihash;
   std::unique_ptr<Mesh<Dim, 1>> hull;

};

template<>
class Mesh<3, 3> : public Mesh<3, 2>
{
public:
   Mesh(std::vector<Eigen::Matrix<double, 3, 1>>* points);

   std::size_t getNumBodies() const override;

   MeshElement* getBody(std::size_t idx) const override;

   MeshElement* getBodyByID(std::size_t id) const override;

   std::size_t getNumFacets() const override;

   MeshElement* getFacet(std::size_t idx) const override;

   MeshElement* getFacetByID(std::size_t id) const override;

   std::size_t getNumRidges() const override;

   MeshElement* getRidge(std::size_t idx) const override;

   MeshElement* getRidgeByID(std::size_t id) const override;

   std::size_t getNumPeaks() const override;

   MeshElement* getPeak(std::size_t idx) const override;

   MeshElement* getPeakByID(std::size_t id) const override;

   Cell* getCell(std::size_t idx) const;

   Cell* getCellByID(std::size_t id) const;

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

/*template<uint Dim, uint TopDim>
class Boundary : public Mesh<Dim, TopDim - 1>
{
public:
   template<uint MeshTopDim>
   Boundary(Mesh<Dim, MeshTopDim>* mesh) : Mesh<Dim, TopDim - 1>(mesh)
   {
      typedef Eigen::Matrix<double, Dim, 1> VectorType;
      std::size_t mostleftbottom = 0;
      const VectorType& p = mesh.getPoint(Mesh<Dim, TopDim - 1>::refvertices[0]);
      double x1 = p[0];
      double y1 = p[1];
      for (std::size_t iv = 1; iv < Mesh<Dim, TopDim - 1>::refvertices.size(); iv++) {
         const VectorType& p = mesh.getPoint(Mesh<Dim, TopDim - 1>::refvertices[iv]);
         const double xdiff = p[0] - x1;
         const bool ylower = p[1] < y1;
         if (ylower)
            y1 = p[1];
         if (abs(xdiff) <= epsilon) {
            if (ylower)
               mostleftbottom = i;
         } else if (xdiff < 0.0) {
            x1 = p[0];
            mostleftbottom = i;
         }
      }
      vector<size_t> bvertices(boundary->getNumVertices());
      iota(bvertices.begin(), bvertices.end(), 0);
      swap(bvertices[0], bvertices[mostleftbottom]);
      sort(bvertices.begin() + 1, bvertices.end(), [&p, &maxcolin, &pstart, this](int a, int b) {
         const double o = Delaunay2D::orientation(p, pts[a], pts[b]);
         if (o == 0.0) {
            const double pa = (pts[a] - p).squaredNorm();
            const double pb = (pts[b] - p).squaredNorm();
            const double m = max(pa, pb);
            const bool res = pa < pb;
            if (m > maxcolin) {
               maxcolin = m;
               pstart = res ? b : a;
            }
            return res;
         }
         return o > 0.0;
      });
   }

private:
   static constexpr double epsilon = 1.0e-10;

   static double orientation(const Eigen::Vector2d &p,
                             const Eigen::Vector2d &q,
                             const Eigen::Vector2d &i);


};*/

template<uint Dim, uint TopDim>
class IMesher
{
public:
   virtual void generate(Mesh<Dim, TopDim>& mesh) = 0;
};

#endif //LBM_POINT_H
