//#pragma SWIG nowarn=315
%{
#include <type_traits>
#include "mesh.h"
%}

%include <std_array.i>
%include <std_vector.i>
%include "attributes.i"
%include "mesh.h"

%extend Simplex
{
   int __getitem__(std::size_t i)
   {
      return (*$self)[i];
   }
};
%ignore Simplex<0, 1, 1>::getFacet;
%ignore Simplex<0, 2, 2>::getFacet;
%ignore Simplex<0, 3, 3>::getFacet;
%template(Simplex0D1D) Simplex<0, 1, 1>;
%template(Simplex0D2D) Simplex<0, 2, 2>;
%template(Simplex0D3D) Simplex<0, 3, 3>;
%template(Simplex1D1D) Simplex<1, 1, 1>;
%template(Simplex1D2D) Simplex<1, 2, 2>;
%template(Simplex1D3D) Simplex<1, 3, 3>;
%template(Simplex2D2D) Simplex<2, 2, 2>;
%template(Simplex2D3D) Simplex<2, 3, 3>;
%template(Simplex3D3D) Simplex<3, 3, 3>;

%template(Vertex1D) Vertex<1, 1>;
%template(Vertex1DList) std::vector<Vertex<1, 1>>;
%template(Vertex2D) Vertex<2, 2>;
%template(Vertex2DList) std::vector<Vertex<2, 2>>;
%template(Vertex3D) Vertex<3, 3>;
%template(Vertex3DList) std::vector<Vertex<3, 3>>;
%template(Edge1D) Edge<1, 1>;
%template(Edge1DList) std::vector<Edge<1, 1>>;
%template(Edge2D) Edge<2, 2>;
%template(Edge2DList) std::vector<Edge<2, 2>>;
%template(Edge3D) Edge<3, 3>;
%template(Edge3DList) std::vector<Edge<3, 3>>;
%template(Face2D) Face<2, 2>;
%template(Face2DList) std::vector<Face<2, 2>>;
%template(Face3D) Face<3, 3>;
%template(Face3DList) std::vector<Face<3, 3>>;
%template(Cell3D) Cell<3>;
%template(Cell3DList) std::vector<Cell<3>>;

//%template(FacetsVertex) std::array<int, 1>;
%template(FacetsEdge1D) std::array<Vertex<1, 1>, 2>;
%template(FacetsEdge2D) std::array<Vertex<2, 2>, 2>;
%template(FacetsEdge3D) std::array<Vertex<3, 3>, 2>;
%template(FacetsFace2D) std::array<Edge<2, 2>, 3>;
%template(FacetsFace3D) std::array<Edge<3, 3>, 3>;
%template(FacetsCell3D) std::array<Face<3, 3>, 4>;

%ignore Mesh<1, 1>::vertices;
%ignore Mesh<1, 1>::edges;
%ignore Mesh<1, 1>::edge2edge;
%ignore Mesh<1, 1>::point2edge;
%ignore Mesh<2, 1>::vertices;
%ignore Mesh<2, 1>::edges;
%ignore Mesh<2, 1>::edge2edge;
%ignore Mesh<2, 1>::point2edge;

%ignore Mesh<2, 2>::vertices;
%ignore Mesh<2, 2>::edges;
%ignore Mesh<2, 2>::edge2edge;
%ignore Mesh<2, 2>::point2edge;
%ignore Mesh<2, 2>::faces;
%ignore Mesh<2, 2>::face2face;
%ignore Mesh<2, 2>::edge2face;
%ignore Mesh<2, 2>::point2face;

%ignore Mesh<3, 3>::vertices;
%ignore Mesh<3, 3>::edges;
%ignore Mesh<3, 3>::edge2edge;
%ignore Mesh<3, 3>::point2edge;
%ignore Mesh<3, 3>::faces;
%ignore Mesh<3, 3>::face2face;
%ignore Mesh<2, 3>::edge2face;
%ignore Mesh<3, 3>::point2face;
%ignore Mesh<3, 3>::bodies;
%ignore Mesh<3, 3>::point2Cell;
%ignore Mesh<3, 3>::edge2Cell;
%ignore Mesh<3, 3>::face2Cell;
%ignore Mesh<3, 2>::vertices;
%ignore Mesh<3, 2>::edges;
%ignore Mesh<3, 2>::edge2edge;
%ignore Mesh<3, 2>::point2edge;
%ignore Mesh<3, 2>::faces;
%ignore Mesh<3, 2>::face2face;
%ignore Mesh<3, 2>::edge2face;
%ignore Mesh<3, 2>::point2face;
%ignore Mesh<3, 1>::vertices;
%ignore Mesh<3, 1>::edges;
%ignore Mesh<3, 1>::edge2edge;
%ignore Mesh<3, 1>::point2edge;
%template(Mesh1D) Mesh<1, 1>;
%template(Mesh2D1D) Mesh<2, 1>;
%template(Mesh2D) Mesh<2, 2>;
%template(Mesh3D1D) Mesh<3, 1>;
%template(Mesh3D2D) Mesh<3, 2>;
%template(Mesh3D) Mesh<3, 3>;
%extend Mesh<2, 2> {
   %template(getOrCreateDoubleAttributeOnVertex) getOrCreateAttributeOnVertex<double>;
   %template(getOrCreateIntListAttributeOnVertex) getOrCreateAttributeOnVertex<std::vector<int>>;
   %template(getOrCreatePoint2DAttributeOnFace) getOrCreateAttributeOnFace<Eigen::Matrix<double, 2, 1>>;
   %template(getOrCreatePoint2DAttributeOnEdge) getOrCreateAttributeOnEdge<Eigen::Matrix<double, 2, 1>>;
   %template(getOrCreatePoint2DListAttributeOnVertex) getOrCreateAttributeOnVertex<std::vector<Eigen::Matrix<double, 2, 1>>>;
};
