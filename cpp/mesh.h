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
#include <sstream>

#include "attributes.h"
#include "eigen.h"
#include "utils.h"

typedef int ID;

class IMesh;

template<uint Dim, uint TopDim>
class Mesh;

class Meshing;

class MeshElement
{
public:
   //MeshElement(IMesh* mesh) : mesh(mesh) {}

   virtual std::size_t getNumVertices() const = 0;

   virtual ID getID() const = 0;

   virtual ID operator[](std::size_t idx) const = 0;

   /*const Eigen::MatrixBase<double>& getPoint(std::size_t idx) const
   {
      mesh.
   }

protected:
   IMesh* mesh;*/

};

template<uint Dim, int SimplexDim>
class Simplex;

#ifndef SWIG
template <uint Dim>
using Vertex = Simplex<Dim, 0>;

template <uint Dim>
using Edge = Simplex<Dim, 1>;

template <uint Dim>
using Face = Simplex<Dim, 2>;

using Cell = Simplex<3, 3>;
#endif

template<uint Dim, int SimplexDim>
struct SimplexTraits
{};

template<uint Dim>
struct SimplexTraits<Dim, 0>
{
   typedef int FacetElementType;
   typedef ID FacetsType;

   static constexpr std::vector<ID>* getVertexList(Mesh<Dim, 0>* mesh)
   {
      return nullptr;
   }
};

template<uint Dim>
struct SimplexTraits<Dim, 1>
{
   typedef Vertex<Dim>* FacetElementType;
   typedef std::array<FacetElementType, 2> FacetsType;

   static constexpr std::vector<ID>* getVertexList(Mesh<Dim, 1>* mesh)
   {
      return mesh->edges_lst;
   }
};

template<uint Dim>
struct SimplexTraits<Dim, 2>
{
   typedef Edge<Dim>* FacetElementType;
   typedef std::array<FacetElementType, 3> FacetsType;

   static constexpr std::vector<ID>* getVertexList(Mesh<Dim, 2>* mesh)
   {
      return mesh->faces_lst;
   }
};

template<uint Dim>
struct SimplexTraits<Dim, 3>
{
   typedef Face<Dim>* FacetElementType;
   typedef std::array<FacetElementType, 4> FacetsType;

   static constexpr std::vector<ID>* getVertexList(Mesh<Dim, 3>* mesh)
   {
      return mesh->cells_lst;
   }
};

template<uint Dim, int SimplexDim>
class Simplex : public MeshElement
{
   static_assert(SimplexDim <= Dim, "Simplex dimension has to be smaller or equal than dimension");
   static_assert(Dim <= 3, "Dimension can only be 1, 2 or 3");


public:
   typedef typename SimplexTraits<Dim, SimplexDim>::FacetElementType FacetElementType;
   typedef typename SimplexTraits<Dim, SimplexDim>::FacetsType FacetsType;
   typedef std::array<ID, SimplexDim + 1> FacetPointListType;

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
   /*Simplex(const FacetsType& facets, ID id)
      :  facets(facets), id(id), pts({0})//, mesh(mesh)
   {
      GetPoints<void, SimplexDim>::get(pts, facets, 0);
   }*/

   Simplex(IMesh* mesh, ID id) : id(id), mesh(mesh), offset(0)
   {
      vertices = SimplexTraits<Dim, SimplexDim>::getVertexList(dynamic_cast<Mesh<Dim, SimplexDim>*>(mesh));
      if (vertices != nullptr) {
         vertices += id * SimplexDim;
         offset = (*vertices)[0];
      }
   }

   ID getID() const override
   {
      return id;
   }

   ID operator[](std::size_t idx) const override
   {
      if (idx >= SimplexDim + 1) {
         std::stringstream ss;
         ss << "Index " << idx << " out of bounds";
         throw std::out_of_range(ss.str());
      }
      if (vertices != nullptr)
         return (*vertices)[idx];
      return id;
   }

   const Eigen::Matrix<double, Dim, 1>& getPoint(std::size_t idx) const
   {
      if (idx >= SimplexDim + 1) {
         std::stringstream ss;
         ss << "Index " << idx << " out of bounds";
         throw std::out_of_range(ss.str());
      }
      return reinterpret_cast<Mesh<Dim, 0>*>(mesh)->getPoint(offset + idx);
   }

   /*FacetElementType getFacet(std::size_t idx)
   {
      if (idx >= facets.size()) {
         std::stringstream ss;
         ss << "Index " << idx << " out of bounds";
         throw std::out_of_range(ss.str());
      }
      return facets[idx];
   }

   FacetsType& getFacets()
   {
      return facets;
   }*/

   std::size_t getNumVertices() const override
   {
      return SimplexDim + 1;
   }

protected:
   ID id;
   std::vector<ID>* vertices;
   std::size_t offset;
   IMesh* mesh;

   //FacetsType facets;
   //FacetPointListType pts;


};

class IMesh
{
public:
#ifndef SWIG
   virtual Eigen::Ref<const Eigen::VectorXd> getPoint(std::size_t idx) = 0;
#endif

   virtual const std::vector<double>& getCoords() const = 0;

   virtual std::vector<double>& getCoords() = 0;

   virtual std::size_t getNumVertices() const = 0;

   virtual MeshElement* getVertex(std::size_t idx) const = 0;

   virtual ID getVertexID(std::size_t idx) const = 0;

   virtual std::size_t getVertexIdx(ID id) const = 0;

   virtual MeshElement* getVertexByID(ID id) const = 0;

   virtual std::size_t getNumBodies() const = 0;

   virtual MeshElement* getBody(std::size_t idx) const = 0;

   virtual MeshElement* getBodyByID(ID id) const = 0;

   virtual std::size_t getNumFacets() const = 0;

   virtual MeshElement* getFacet(std::size_t idx) const = 0;

   virtual MeshElement* getFacetByID(ID id) const = 0;

   virtual std::size_t getNumRidges() const = 0;

   virtual MeshElement* getRidge(std::size_t idx) const = 0;

   virtual MeshElement* getRidgeByID(ID id) const = 0;

   virtual std::size_t getNumPeaks() const = 0;

   virtual MeshElement* getPeak(std::size_t idx) const = 0;

   virtual MeshElement* getPeakByID(ID id) const = 0;

   virtual IMesh* getHull() const = 0;

   virtual BaseAttribute* getAttributeOnBody(const std::string &name) const = 0;

   virtual BaseAttribute* getAttributeOnFacet(const std::string &name) const = 0;

   virtual BaseAttribute* getAttributeOnRidge(const std::string &name) const = 0;

   virtual BaseAttribute* getAttributeOnPeak(const std::string &name) const = 0;

};

template<uint Dim, uint TopDim=Dim>
class Mesh
{};

template<uint Dim>
class Mesh<Dim, 0> : public IMesh
{
   static_assert(Dim >= 0 && Dim <= 3, "Dimension not supported");
   friend SimplexTraits<Dim, 0>;
   friend Meshing;

public:
   Mesh(const std::vector<double>& coords);

   template<uint TopDim>
   Mesh(Mesh<Dim, TopDim>* mesh) : points(mesh->points)
   {
      static_assert(TopDim >= 0, "Dimension mismatch");
      vertices = mesh->vertices;
   }

   virtual std::size_t getNumBodies() const override;

   virtual MeshElement* getBody(std::size_t idx) const override;

   virtual MeshElement* getBodyByID(ID id) const override;

   virtual std::size_t getNumFacets() const override;

   virtual MeshElement* getFacet(std::size_t idx) const override;

   virtual MeshElement* getFacetByID(ID id) const override;

   virtual std::size_t getNumRidges() const override;

   virtual MeshElement* getRidge(std::size_t idx) const override;

   virtual MeshElement* getRidgeByID(ID id) const override;

   virtual std::size_t getNumPeaks() const override;

   virtual MeshElement* getPeak(std::size_t idx) const override;

   virtual MeshElement* getPeakByID(ID id) const override;

   virtual BaseAttribute* getAttributeOnBody(const std::string &name) const override;

   virtual BaseAttribute* getAttributeOnFacet(const std::string &name) const override;

   virtual BaseAttribute* getAttributeOnRidge(const std::string &name) const override;

   virtual BaseAttribute* getAttributeOnPeak(const std::string &name) const override;

   const Eigen::Matrix<double, Dim, 1>& getPoint(std::size_t idx) const;

   //std::vector<Eigen::Matrix<double, Dim, 1>>* getPoints() const;

   virtual const std::vector<double>& getCoords() const override;

   virtual std::vector<double>& getCoords() override;

#ifndef SWIG
   Eigen::Ref<const Eigen::VectorXd> getPoint(std::size_t idx) override;
#endif

   //Vertex<Dim>* getVertex(std::size_t idx);

   MeshElement* getVertex(std::size_t idx) const override;

   ID getVertexID(std::size_t idx) const override;

   std::size_t getVertexIdx(ID id) const override;

   //Vertex<Dim>* getVertexByID(ID id);

   MeshElement* getVertexByID(ID id) const override;

   Vertex<Dim>* getOrCreateVertexByID(ID id);

   std::size_t getNumVertices() const override;

   virtual IMesh* getHull() const override;

   template<class T>
   Attribute<T>* getAttributeOnVertex(const std::string &name)
   {
      return Mesh<Dim, 0>::template getAttribute<T>(name, vertexAttrs, refvertices.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnVertex(const std::string &name)
   {
      return Mesh<Dim, 0>::template getOrCreateAttribute<T>(name, vertexAttrs, refvertices.size());
   }

   template<class T>
   ListAttribute<T>* getListAttributeOnVertex(const std::string &name)
   {
      return Mesh<Dim, 0>::template getAttribute<T, ListAttribute>(name, vertexAttrs, refvertices.size());
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
   static inline BaseAttribute *getAttribute(const std::string &name,
                                             const std::unordered_map<std::string,
                                                     std::unique_ptr<BaseAttribute>> &attrlst,
                                             std::size_t size)
   {
      const auto it = attrlst.find(name);
      if (it == attrlst.end())
         return nullptr;
      return it->second.get();
   }

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

   template<typename T, template<typename> typename AttibuteType = Attribute>
   static inline AttibuteType<T> *getAttribute(const std::string &name,
                                               std::unordered_map<std::string,
                                                       std::unique_ptr<BaseAttribute>> &attrlst,
                                               std::size_t size)
   {
      const auto it = attrlst.find(name);
      if (it == attrlst.end())
         return nullptr;
      return static_cast<AttibuteType<T> *>(it->second.get());
   }

protected:
   // Coordinates
   std::vector<double> coords;
   std::vector<Eigen::Matrix<double, Dim, 1>> points;

   //Referenced elements
   std::vector<ID> refvertices;

   //Elements
   std::unique_ptr<std::vector<std::unique_ptr<Vertex<Dim>>>> vertices_owner;
   std::vector<std::unique_ptr<Vertex<Dim>>>* vertices;

   //Atributes
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> vertexAttrs;

};

template<uint Dim>
class Mesh<Dim, 1> : public Mesh<Dim, 0>
{
   static_assert(Dim >= 1, "Topological dimension not supported");
   friend SimplexTraits<Dim, 1>;
   friend Meshing;

public:
   /**
    *
    * @param points
    */
   Mesh(const std::vector<double>& coords);

   /**
    *
    * @tparam TopDim The topological dimension of the parent mesh
    * @param mesh The parent mesh
    */
   template<uint TopDim>
   Mesh(Mesh<Dim, TopDim>* mesh) : Mesh<Dim, 0>(mesh)
   {
      static_assert(TopDim >= 1, "Dimension mismatch");
      edges_lst = mesh->edges_lst;
      edges = mesh->edges;
      vhash2id = mesh->vhash2id;
      vertex2edge = mesh->vertex2edge;
      hull = std::make_unique<Mesh<Dim, 0>>(this);
   }

   virtual std::size_t getNumBodies() const override;

   virtual MeshElement* getBody(std::size_t idx) const override;

   virtual MeshElement* getBodyByID(ID id) const override;

   virtual std::size_t getNumFacets() const override;

   virtual MeshElement* getFacet(std::size_t idx) const override;

   virtual MeshElement* getFacetByID(ID id) const override;

   virtual std::size_t getNumRidges() const override;

   virtual MeshElement* getRidge(std::size_t idx) const override;

   virtual MeshElement* getRidgeByID(ID id) const override;

   virtual std::size_t getNumPeaks() const override;

   virtual MeshElement* getPeak(std::size_t idx) const override;

   virtual MeshElement* getPeakByID(ID id) const override;

   virtual BaseAttribute* getAttributeOnBody(const std::string &name) const override;

   virtual BaseAttribute* getAttributeOnFacet(const std::string &name) const override;

   virtual BaseAttribute* getAttributeOnRidge(const std::string &name) const override;

   virtual BaseAttribute* getAttributeOnPeak(const std::string &name) const override;

   /**
    * Returns the i-th edge of this mesh.
    *
    * @param idx The index
    * @return The pointer to the edge if found, otherwise
    * nullptr
    */
   Simplex<Dim, 1>* getEdge(std::size_t idx) const;

   /**
    * Returns the edge with the ID @param id
    *
    * @param id The ID of the edge
    * @return The pointer to the edge if found, otherwise
    * nullptr
    */
   Simplex<Dim, 1>* getEdgeByID(ID id) const;

   /**
    *
    * @param id
    * @return
    */
   Simplex<Dim, 1>* getOrCreateEdgeByID(ID id);

   Simplex<Dim, 1>* getEdge(ID vid1, ID vid2) const;

   Simplex<Dim, 1>* getOrCreateEdge(ID vid1, ID vid2);

   std::size_t getNumEdges() const;

   virtual IMesh* getHull() const override;

   std::pair<typename std::unordered_multimap<ID, Edge<Dim>*>::const_iterator,
           typename std::unordered_multimap<ID, Edge<Dim>*>::const_iterator> getEdgesOfVertex(std::size_t idx) const;

   template<class T>
   Attribute<T>* getAttributeOnEdge(const std::string &name)
   {
      return Mesh<Dim, 0>::template getAttribute<T>(name, edgeAttrs, refedges.size());
   }

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
   Edge<Dim>* insertEdge(ID id);

   //Edge lists
   std::unique_ptr<std::vector<ID>> edges_lst_owner;
   std::vector<ID>* edges_lst;

   //Referenced elements
   std::vector<ID> refedges;

   // Elements
   std::unique_ptr<std::vector<std::unique_ptr<Edge<Dim>>>> edges_owner;
   std::vector<std::unique_ptr<Edge<Dim>>>* edges;

   // Neighbor relations
   LineHash linhash;
   std::unique_ptr<std::unordered_map<ullong, ID>> vhash2id_owner;
   std::unordered_map<ullong, ID>* vhash2id;
   std::unique_ptr<std::unordered_multimap<ID, Edge<Dim>*>> vertex2edge_owner;
   std::unordered_multimap<ID, Edge<Dim>*>* vertex2edge;

   // Attributes
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> edgeAttrs;

   // Sub meshes
   std::unique_ptr<Mesh<Dim, 0>> hull;

};

template<uint Dim>
class Mesh<Dim, 2> : public Mesh<Dim, 1>
{
   static_assert(Dim >= 2, "Topological dimension not supported");
   friend SimplexTraits<Dim, 2>;
   friend Meshing;

public:
   Mesh(const std::vector<double>& coords);

   template<uint TopDim>
   Mesh(Mesh<Dim, TopDim>* mesh) : Mesh<Dim, 1>(mesh)
   {
      static_assert(TopDim >= 2, "Dimension mismatch");
      faces_lst = mesh->faces_lst;
      faces = mesh->faces;
      ehash2id = mesh->ehash2id;
      vertex2face = mesh->vertex2face;
      edge2face = mesh->edge2face;
      hull = std::make_unique<Mesh<Dim, 1>>(this);
   }

   virtual std::size_t getNumBodies() const override;

   virtual MeshElement* getBody(std::size_t idx) const override;

   virtual MeshElement* getBodyByID(ID id) const override;

   virtual std::size_t getNumFacets() const override;

   virtual MeshElement* getFacet(std::size_t idx) const override;

   virtual MeshElement* getFacetByID(ID id) const override;

   virtual std::size_t getNumRidges() const override;

   virtual MeshElement* getRidge(std::size_t idx) const override;

   virtual MeshElement* getRidgeByID(ID id) const override;

   virtual std::size_t getNumPeaks() const override;

   virtual MeshElement* getPeak(std::size_t idx) const override;

   virtual MeshElement* getPeakByID(ID id) const override;

   virtual BaseAttribute* getAttributeOnBody(const std::string &name) const override;

   virtual BaseAttribute* getAttributeOnFacet(const std::string &name) const override;

   virtual BaseAttribute* getAttributeOnRidge(const std::string &name) const override;

   virtual BaseAttribute* getAttributeOnPeak(const std::string &name) const override;

   Simplex<Dim, 2>* getFace(std::size_t idx) const;

   Simplex<Dim, 2>* getFaceByID(ID id) const;

   Simplex<Dim, 2>* getFace(ID eid1, ID eid2, ID eid3) const;

   Simplex<Dim, 2>* getOrCreateFace(ID eid1, ID eid2, ID eid3);

   std::size_t getNumFaces() const;

   std::pair<typename std::unordered_multimap<ID, Face<Dim>*>::const_iterator,
   typename std::unordered_multimap<ID, Face<Dim>*>::const_iterator> getFacesOfEdge(Edge<Dim>* edge) const;

   virtual IMesh* getHull() const override;

   template<class T>
   Attribute<T>* getAttributeOnFace(const std::string &name)
   {
      return Mesh<Dim, 0>::template getAttribute<T>(name, facesAttrs, reffaces.size());
   }

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
   Face<Dim>* insertFace(ID id);

   //Face list
   std::unique_ptr<std::vector<ID>> faces_lst_owner;
   std::vector<ID>* faces_lst;

   // Referenced elements
   std::vector<ID> reffaces;

   // Elements
   std::unique_ptr<std::vector<std::unique_ptr<Face<Dim>>>> faces_owner;
   std::vector<std::unique_ptr<Face<Dim>>>* faces;

   // Neighbor relations
   TriangleHash trihash;
   std::unique_ptr<std::unordered_map<ullong, ID>> ehash2id_owner;
   std::unordered_map<ullong, ID>* ehash2id;
   std::unique_ptr<std::unordered_multimap<ID, Face<Dim>*>> vertex2face_owner;
   std::unordered_multimap<ID, Face<Dim>*>* vertex2face;
   std::unique_ptr<std::unordered_multimap<ID, Face<Dim>*>> edge2face_owner;
   std::unordered_multimap<ID, Face<Dim>*>* edge2face;

   // Attributes
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> facesAttrs;

   // Sub meshes
   std::unique_ptr<Mesh<Dim, 1>> hull;

};

template<>
class Mesh<3, 3> : public Mesh<3, 2>
{
   friend SimplexTraits<3, 3>;
   friend Meshing;

public:
   //Mesh(vector<Matrix<double, 3, 1>>* points);
   Mesh(const std::vector<double>& coords);

   std::size_t getNumBodies() const override;

   MeshElement* getBody(std::size_t idx) const override;

   MeshElement* getBodyByID(ID id) const override;

   std::size_t getNumFacets() const override;

   MeshElement* getFacet(std::size_t idx) const override;

   MeshElement* getFacetByID(ID id) const override;

   std::size_t getNumRidges() const override;

   MeshElement* getRidge(std::size_t idx) const override;

   MeshElement* getRidgeByID(ID id) const override;

   std::size_t getNumPeaks() const override;

   MeshElement* getPeak(std::size_t idx) const override;

   MeshElement* getPeakByID(ID id) const override;

   BaseAttribute* getAttributeOnBody(const std::string &name) const override;

   BaseAttribute* getAttributeOnFacet(const std::string &name) const override;

   BaseAttribute* getAttributeOnRidge(const std::string &name) const override;

   BaseAttribute* getAttributeOnPeak(const std::string &name) const override;

   Simplex<3, 3>* getCell(std::size_t idx) const;

   Simplex<3, 3>* getCellByID(ID id) const;

   std::size_t getNumCells() const;

   template<class T>
   Attribute<T>* getAttributeOnCell(const std::string &name)
   {
      return Mesh<3, 0>::template getAttribute<T>(name, cellsAttrs, cells.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnCell(const std::string &name)
   {
      return Mesh<3, 0>::template getOrCreateAttribute<T>(name, cellsAttrs, cells.size());
   }

   template<class T>
   Attribute<T>* getOrCreateAttributeOnBody(const std::string &name, const T& def)
   {
      return Mesh<3, 0>::template getOrCreateAttribute<T>(name, cellsAttrs, cells.size(), def);
   }

protected:
   std::unique_ptr<ID[]> cells_lst_owner;
   std::vector<ID>* cells_lst;

   std::vector<std::unique_ptr<Cell>> cells;
   std::unordered_multimap<ullong, Cell*> point2body;
   std::unordered_multimap<ullong, Cell*> edge2body;
   std::unordered_multimap<ullong, Cell*> face2body;
   std::unordered_map<std::string, std::unique_ptr<BaseAttribute>> cellsAttrs;
   Mesh<3, 2>* hull;

};

#endif //LBM_POINT_H
