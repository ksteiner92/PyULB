//
// Created by klaus on 06.01.19.
//

#include <sstream>

#include "mesh.h"

using namespace std;
using namespace Eigen;

template<uint Dim>
Mesh<Dim, 0>::Mesh(vector<Matrix<double, Dim, 1>>* points) : points(points)
{
   vertices_owner = make_unique<vector<unique_ptr<Vertex<Dim>>>>();
   vertices = vertices_owner.get();
   const size_t npts = points->size();
   vertices->reserve(npts);
   /*vertices->resize(npts);
   refvertices.resize(npts);
   for (size_t iv = 0; iv < npts; iv++) {
      (*vertices)[iv] = make_unique<Vertex<Dim>>(Vertex<Dim>(iv, iv));
      refvertices[iv] = iv;
   }*/
}

template<uint Dim>
size_t Mesh<Dim, 0>::getNumBodies() const
{
   return getNumVertices();
}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::getBody(size_t idx) const
{
   return getVertex(idx);
}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::getBodyByID(size_t id) const
{
   return getVertexByID(id);
}

template<uint Dim>
size_t Mesh<Dim, 0>::getNumFacets() const
{
   throw logic_error("Mesh with 0 topology does not have facets");
}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::getFacet(size_t idx) const
{
   throw logic_error("Mesh with 0 topology does not have facets");
}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::getFacetByID(size_t id) const
{
   throw logic_error("Mesh with 0 topology does not have facets");
}

template<uint Dim>
size_t Mesh<Dim, 0>::getNumRidges() const
{
   throw logic_error("Mesh with 0 topology does not have ridges");
}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::getRidge(size_t idx) const
{
   throw logic_error("Mesh with 0 topology does not have ridges");
}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::getRidgeByID(size_t id) const
{
   throw logic_error("Mesh with 0 topology does not have ridges");
}

template<uint Dim>
size_t Mesh<Dim, 0>::getNumPeaks() const
{
   throw logic_error("Mesh with 0 topology does not have peaks");
}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::getPeak(size_t idx) const
{
   throw logic_error("Mesh with 0 topology does not have peaks");
}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::getPeakByID(size_t id) const
{
   throw logic_error("Mesh with 0 topology does not have peaks");
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 0>::getAttributeOnBody(const string &name) const
{
   return Mesh<Dim, 0>::getAttribute(name, vertexAttrs, refvertices.size());
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 0>::getAttributeOnFacet(const string &name) const
{
   throw logic_error("Mesh with topological dimension 0 does not have facets");
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 0>::getAttributeOnRidge(const string &name) const
{
   throw logic_error("Mesh with topological dimension 0 does not have ridges");
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 0>::getAttributeOnPeak(const string &name) const
{
   throw logic_error("Mesh with topological dimension 0 does not have peaks");
}

template<uint Dim>
const Matrix<double, Dim, 1>& Mesh<Dim, 0>::getPoint(size_t idx) const
{
   if (idx >= points->size())
      throw out_of_range("Index out of range");
   return (*points)[idx];
}

template<uint Dim>
vector<Matrix<double, Dim, 1>>* Mesh<Dim, 0>::getPoints() const
{
   return points;
}

template<uint Dim>
size_t Mesh<Dim, 0>::getNumPoints() const
{
   return points->size();
}

template<uint Dim>
Vertex<Dim>* Mesh<Dim, 0>::getVertex(size_t idx) const
{
   if (idx < refvertices.size())
      return (*vertices)[refvertices[idx]].get();
   return nullptr;
}

template<uint Dim>
size_t Mesh<Dim, 0>::getVertexID(size_t idx) const
{
   return refvertices[idx];
}

template<uint Dim>
size_t Mesh<Dim, 0>::getVertexIdx(size_t id) const
{
   const auto it = binary_find(refvertices.begin(), refvertices.end(), id);
   if (it == refvertices.end())
      throw logic_error("Vertex not found");
   return distance(refvertices.begin(), it);
}

template<uint Dim>
Vertex<Dim>* Mesh<Dim, 0>::getVertexByID(size_t id) const
{
   if (!binary_search(refvertices.begin(), refvertices.end(), id))
      return nullptr;
   if (id < vertices->size())
      return (*vertices)[id].get();
   return nullptr;
}

template<uint Dim>
Vertex<Dim>* Mesh<Dim, 0>::getOrCreateVertexByID(size_t id)
{
   if (!binary_search(refvertices.begin(), refvertices.end(), id)) {
      if (id < vertices->size()) {
         insertSorted(refvertices, id);
         return (*vertices)[id].get();
      }
      if (id >= points->size())
         throw out_of_range("Vertex index out of range");
      vertices->push_back(make_unique<Vertex<Dim>>(Vertex<Dim>(id, id)));
      Vertex<Dim>* vertex = vertices->back().get();
      insertSorted(refvertices, id);
      return vertex;
   }
   return (*vertices)[id].get();
}

template<uint Dim>
size_t Mesh<Dim, 0>::getNumVertices() const
{
   return refvertices.size();
}

template<uint Dim>
IMesh* Mesh<Dim, 0>::getHull() const
{
   throw logic_error("Mesh of 0 topological dimension has no hull");
   return nullptr;
}

template<uint Dim>
Mesh<Dim, 1>::Mesh(vector<Matrix<double, Dim, 1>>* points) : Mesh<Dim, 0>(points)
{
   edges_owner = make_unique<vector<unique_ptr<Edge<Dim>>>>();
   edges = edges_owner.get();
   vhash2id_owner = make_unique<unordered_map<ullong, size_t>>();
   vhash2id = vhash2id_owner.get();
   vertex2edge_owner = make_unique<unordered_multimap<size_t, Edge<Dim>*>>();
   vertex2edge = vertex2edge_owner.get();
   hull = make_unique<Mesh<Dim, 0>>(this);
}

template<uint Dim>
size_t Mesh<Dim, 1>::getNumBodies() const
{
   return getNumEdges();
}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::getBody(size_t idx) const
{
   return getEdge(idx);
}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::getBodyByID(size_t id) const
{
   return getEdgeByID(id);
}

template<uint Dim>
size_t Mesh<Dim, 1>::getNumFacets() const
{
   return Mesh<Dim, 0>::getNumVertices();
}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::getFacet(size_t idx) const
{
   return Mesh<Dim, 0>::getVertex(idx);
}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::getFacetByID(size_t id) const
{
   return Mesh<Dim, 0>::getVertexByID(id);
}

template<uint Dim>
size_t Mesh<Dim, 1>::getNumRidges() const
{
   throw logic_error("Mesh with topological dimension does not have ridges");
}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::getRidge(size_t idx) const
{
   throw logic_error("Mesh with topological dimension 1 does not have ridges");
}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::getRidgeByID(size_t id) const
{
   throw logic_error("Mesh with topological dimension 1 does not have ridges");
}

template<uint Dim>
size_t Mesh<Dim, 1>::getNumPeaks() const
{
   throw logic_error("Mesh with topological dimension 1 does not have peaks");
}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::getPeak(size_t idx) const
{
   throw logic_error("Mesh with topological dimension 1 does not have peaks");
}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::getPeakByID(size_t id) const
{
   throw logic_error("Mesh with topological dimension 1 does not have peaks");
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 1>::getAttributeOnBody(const string &name) const
{
   return Mesh<Dim, 0>::getAttribute(name, edgeAttrs, refedges.size());
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 1>::getAttributeOnFacet(const string &name) const
{
   return Mesh<Dim, 0>::getAttributeOnBody(name);
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 1>::getAttributeOnRidge(const string &name) const
{
   throw logic_error("Mesh with topological dimension 1 does not have ridges");
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 1>::getAttributeOnPeak(const string &name) const
{
   throw logic_error("Mesh with topological dimension 1 does not have peaks");
}

template<uint Dim>
Edge<Dim>* Mesh<Dim, 1>::getEdge(size_t idx) const
{
   if (idx < refedges.size())
      return (*edges)[refedges[idx]].get();
   return nullptr;
}

template<uint Dim>
Edge<Dim>* Mesh<Dim, 1>::getEdgeByID(size_t id) const
{
   if (!binary_search(refedges.begin(), refedges.end(), id))
      return nullptr;
   if (id < edges->size())
      return (*edges)[id].get();
   return nullptr;
}

template<uint Dim>
Edge<Dim>* Mesh<Dim, 1>::getEdge(size_t vid1, size_t vid2) const
{
   auto it = vhash2id->find(linhash(make_tuple(vid1, vid2)));
   if (it != vhash2id->end())
      return getEdgeByID(it->second);
   else {
      it = vhash2id->find(linhash(make_tuple(vid2, vid1)));
      if (it != vhash2id->end())
         return getEdgeByID(it->second);
   }
   return nullptr;
}

template<uint Dim>
Edge<Dim> *Mesh<Dim, 1>::getOrCreateEdgeByID(size_t id)
{
   if (!binary_search(refedges.begin(), refedges.end(), id)) {
      if (id < edges->size()) {
         insertSorted(refedges, id);
         Edge<Dim>* edge = (*edges)[id].get();
         Mesh<Dim, 0>::getOrCreateVertexByID((*edge)[0]);
         Mesh<Dim, 0>::getOrCreateVertexByID((*edge)[1]);
         return edge;
      } else {
         stringstream ss;
         ss << "No edge with id " << id << " known";
         throw runtime_error(ss.str());
      }
   }
   return (*edges)[id].get();
}

template<uint Dim>
Edge<Dim> *Mesh<Dim, 1>::getOrCreateEdge(size_t vid1, size_t vid2)
{
   const ullong key = linhash(make_tuple(vid1, vid2));
   auto it = vhash2id->find(key);
   if (it != vhash2id->end()) {
      if (!binary_search(refedges.begin(), refedges.end(), it->second)) {
         Mesh<Dim, 0>::getOrCreateVertexByID(vid1);
         Mesh<Dim, 0>::getOrCreateVertexByID(vid2);
         insertSorted(refedges, it->second);
      }
      return (*edges)[it->second].get();
   } else {
      it = vhash2id->find(linhash(make_tuple(vid2, vid1)));
      if (it != vhash2id->end()) {
         if (!binary_search(refedges.begin(), refedges.end(), it->second)) {
            Mesh<Dim, 0>::getOrCreateVertexByID(vid1);
            Mesh<Dim, 0>::getOrCreateVertexByID(vid2);
            insertSorted(refedges, it->second);
         }
         return (*edges)[it->second].get();
      }
   }
   Vertex<Dim>* v1 = Mesh<Dim, 0>::getOrCreateVertexByID(vid1);
   Vertex<Dim>* v2 = Mesh<Dim, 0>::getOrCreateVertexByID(vid2);
   const size_t id = edges->size();
   edges->push_back(make_unique<Edge<Dim>>(Edge<Dim>({v1, v2}, id)));
   (*vhash2id)[key] = id;
   Edge<Dim>* edge = edges->back().get();
   vertex2edge->insert(make_pair(vid1, edge));
   vertex2edge->insert(make_pair(vid2, edge));
   insertSorted(refedges, id);
   return edge;
}

template<uint Dim>
size_t Mesh<Dim, 1>::getNumEdges() const
{
   return refedges.size();
}

template<uint Dim>
IMesh* Mesh<Dim, 1>::getHull() const
{
   return hull.get();
}

template<uint Dim>
pair<typename unordered_multimap<size_t, Edge<Dim>*>::const_iterator,
     typename unordered_multimap<size_t, Edge<Dim>*>::const_iterator> Mesh<Dim, 1>::getEdgesOfVertex(size_t idx) const
{
   if (idx >= Mesh<Dim, 0>::refvertices.size())
      throw out_of_range("Vertex index out of range");
   return vertex2edge->equal_range(Mesh<Dim, 0>::refvertices[idx]);
}

template<uint Dim>
Mesh<Dim, 2>::Mesh(vector<Matrix<double, Dim, 1>>* points) : Mesh<Dim, 1>(points)
{
   faces_owner = make_unique<vector<unique_ptr<Face<Dim>>>>();
   faces = faces_owner.get();
   ehash2id_owner = make_unique<unordered_map<ullong, size_t>>();
   ehash2id = ehash2id_owner.get();
   vertex2face_owner = make_unique<unordered_multimap<size_t, Face<Dim>*>>();
   vertex2face = vertex2face_owner.get();
   edge2face_owner = make_unique<unordered_multimap<size_t, Face<Dim>*>>();
   edge2face = edge2face_owner.get();
   hull = make_unique<Mesh<Dim, 1>>(this);
}

template<uint Dim>
size_t Mesh<Dim, 2>::getNumBodies() const
{
   return getNumFaces();
}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::getBody(size_t idx) const
{
   return getFace(idx);
}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::getBodyByID(size_t id) const
{
   return getFaceByID(id);
}

template<uint Dim>
size_t Mesh<Dim, 2>::getNumFacets() const
{
   return Mesh<Dim, 1>::getNumEdges();
}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::getFacet(size_t idx) const
{
   return Mesh<Dim, 1>::getEdge(idx);
}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::getFacetByID(size_t id) const
{
   return Mesh<Dim, 1>::getEdgeByID(id);
}

template<uint Dim>
size_t Mesh<Dim, 2>::getNumRidges() const
{
   return Mesh<Dim, 0>::getNumVertices();
}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::getRidge(size_t idx) const
{
   return Mesh<Dim, 0>::getVertex(idx);
}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::getRidgeByID(size_t id) const
{
   return Mesh<Dim, 0>::getVertexByID(id);
}

template<uint Dim>
size_t Mesh<Dim, 2>::getNumPeaks() const
{
   throw logic_error("Mesh with topological dimension 2 does not have peaks");
}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::getPeak(size_t idx) const
{
   throw logic_error("Mesh with topological dimension 2 does not have peaks");
}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::getPeakByID(size_t id) const
{
   throw logic_error("Mesh with topological dimension 2 does not have peaks");
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 2>::getAttributeOnBody(const string &name) const
{
   return Mesh<Dim, 0>::getAttribute(name, facesAttrs, facesAttrs.size());
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 2>::getAttributeOnFacet(const string &name) const
{
   return Mesh<Dim, 1>::getAttributeOnBody(name);
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 2>::getAttributeOnRidge(const string &name) const
{
   return Mesh<Dim, 0>::getAttributeOnBody(name);
}

template<uint Dim>
BaseAttribute* Mesh<Dim, 2>::getAttributeOnPeak(const string &name) const
{
   throw logic_error("Mesh with topological dimension 2 does not have peaks");
}

template<uint Dim>
Face<Dim>* Mesh<Dim, 2>::getFace(size_t idx) const
{
   if (idx < reffaces.size())
      return (*faces)[reffaces[idx]].get();
   return nullptr;
}

template<uint Dim>
Face<Dim>* Mesh<Dim, 2>::getFaceByID(size_t id) const
{
   if (!binary_search(reffaces.begin(), reffaces.end(), id))
      return nullptr;
   if (id < faces->size())
      return (*faces)[id].get();
   return nullptr;
}

template<uint Dim>
Face<Dim>* Mesh<Dim, 2>::getFace(size_t eid1, size_t eid2, size_t eid3) const
{
   auto it = ehash2id->find(trihash(make_tuple(eid1, eid2, eid3)));
   if (it != ehash2id->end())
      return getFace(it->second);
   return nullptr;
}

template<uint Dim>
Face<Dim> *Mesh<Dim, 2>::getOrCreateFace(size_t eid1, size_t eid2, size_t eid3)
{
   const ullong key = trihash(make_tuple(eid1, eid2, eid3));
   auto it = ehash2id->find(key);
   if (it == ehash2id->end()) {
      Edge<Dim>* e1 = Mesh<Dim, 1>::getOrCreateEdgeByID(eid1);
      Edge<Dim>* e2 = Mesh<Dim, 1>::getOrCreateEdgeByID(eid2);
      Edge<Dim>* e3 = Mesh<Dim, 1>::getOrCreateEdgeByID(eid3);
      const size_t id = faces->size();
      faces->push_back(make_unique<Face<Dim>>(Face<Dim>({e1, e2, e3}, id)));
      (*ehash2id)[key] = id;
      if (!binary_search(reffaces.begin(), reffaces.end(), id))
         insertSorted(reffaces, id);
      Face<Dim>* face = faces->back().get();
      for (size_t iv = 0; iv < face->getNumVertices(); iv++) {
         const size_t vid = (*face)[iv];
         vertex2face->insert(make_pair(vid, face));
         Mesh<Dim, 0>::getOrCreateVertexByID(vid);
      }
      edge2face->insert(make_pair(e1->getID(), face));
      edge2face->insert(make_pair(e2->getID(), face));
      edge2face->insert(make_pair(e3->getID(), face));
      return face;
   }
   Face<Dim>* face = (*faces)[it->second].get();
   const size_t id = face->getID();
   if (!binary_search(reffaces.begin(), reffaces.end(), id)) {
      for (size_t iv = 0; iv < face->getNumVertices(); iv++)
         Mesh<Dim, 0>::getOrCreateVertexByID((*face)[iv]);
      Mesh<Dim, 1>::getOrCreateEdgeByID(eid1);
      Mesh<Dim, 1>::getOrCreateEdgeByID(eid2);
      Mesh<Dim, 1>::getOrCreateEdgeByID(eid3);
      insertSorted(reffaces, id);
   }
   return face;
}

template<uint Dim>
Face<Dim> *Mesh<Dim, 2>::getOrCreateFace(Edge<Dim>* e1, Edge<Dim>* e2, Edge<Dim>* e3)
{
   return getOrCreateFace(e1->getID(), e2->getID(), e3->getID());
}

template<uint Dim>
pair<typename unordered_multimap<size_t, Face<Dim> *>::const_iterator,
     typename unordered_multimap<size_t, Face<Dim> *>::const_iterator>
Mesh<Dim, 2>::getFacesOfEdge(Edge<Dim> *edge) const
{
   return edge2face->equal_range(edge->getID());
}

template<uint Dim>
size_t Mesh<Dim, 2>::getNumFaces() const
{
   return reffaces.size();
}

template<uint Dim>
IMesh* Mesh<Dim, 2>::getHull() const
{
   return hull.get();
}

Mesh<3, 3>::Mesh(vector<Matrix<double, 3, 1>>* points) : Mesh<3, 2>(points)
{
}

size_t Mesh<3, 3>::getNumBodies() const
{
   return getNumCells();
}

MeshElement* Mesh<3, 3>::getBody(size_t idx) const
{
   return getCell(idx);
}

MeshElement* Mesh<3, 3>::getBodyByID(size_t id) const
{
   return getCellByID(id);
}

size_t Mesh<3, 3>::getNumFacets() const
{
   return Mesh<3, 2>::getNumFaces();
}

MeshElement* Mesh<3, 3>::getFacet(size_t idx) const
{
   return Mesh<3, 2>::getFace(idx);
}

MeshElement* Mesh<3, 3>::getFacetByID(size_t id) const
{
   return Mesh<3, 2>::getFaceByID(id);
}

size_t Mesh<3, 3>::getNumRidges() const
{
   return Mesh<3, 1>::getNumEdges();
}

MeshElement* Mesh<3, 3>::getRidge(size_t idx) const
{
   return Mesh<3, 1>::getEdge(idx);
}

MeshElement* Mesh<3, 3>::getRidgeByID(size_t id) const
{
   return Mesh<3, 1>::getEdgeByID(id);
}

size_t Mesh<3, 3>::getNumPeaks() const
{
   return Mesh<3, 0>::getNumVertices();
}

MeshElement* Mesh<3, 3>::getPeak(size_t idx) const
{
   return Mesh<3, 0>::getVertex(idx);
}

MeshElement* Mesh<3, 3>::getPeakByID(size_t id) const
{
   return Mesh<3, 0>::getVertexByID(id);
}

BaseAttribute* Mesh<3, 3>::getAttributeOnBody(const string &name) const
{
   return Mesh<3, 0>::getAttribute(name, cellsAttrs, cellsAttrs.size());
}

BaseAttribute* Mesh<3, 3>::getAttributeOnFacet(const string &name) const
{
   return Mesh<3, 2>::getAttributeOnBody(name);
}

BaseAttribute* Mesh<3, 3>::getAttributeOnRidge(const string &name) const
{
   return Mesh<3, 1>::getAttributeOnBody(name);
}

BaseAttribute* Mesh<3, 3>::getAttributeOnPeak(const string &name) const
{
   return Mesh<3, 0>::getAttributeOnBody(name);
}

Cell* Mesh<3, 3>::getCell(size_t idx) const
{
   if (idx >= cells.size()) {
      cout << "Mesh3D getCell: Index out of bounds" << endl;
      cout.flush();
      throw out_of_range("Cell index out of range");
   }
   return cells[idx].get();
}

Cell* Mesh<3, 3>::getCellByID(size_t id) const
{
   return nullptr;
}

size_t Mesh<3, 3>::getNumCells() const
{
   return cells.size();
}

template class Mesh<0, 0>;
template class Mesh<1, 0>;
template class Mesh<2, 0>;
template class Mesh<3, 0>;
template class Mesh<1, 1>;
template class Mesh<2, 1>;
template class Mesh<2, 2>;
template class Mesh<3, 1>;
template class Mesh<3, 2>;
template class Mesh<3, 3>;