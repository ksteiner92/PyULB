//#pragma SWIG nowarn=315
#define decltype(auto) // \
void ignore_me();
%ignore ignore_me;
%{
#include <type_traits>
#include "mesh.h"

void ignore_me() {}
%}


%include <std_array.i>
%include <std_vector.i>
%include <stdint.i>
%include "attributes.i"
%include "mesh.h"

%extend Simplex
{
   int __getitem__(std::size_t i)
   {
      return (*$self)[i];
   }
};
%ignore Simplex<1, 0>::getFacet;
%ignore Simplex<2, 0>::getFacet;
%ignore Simplex<3, 0>::getFacet;
%template(Simplex0D1D) Simplex<1, 0>;
%template(Simplex0D2D) Simplex<2, 0>;
%template(Simplex0D3D) Simplex<3, 0>;
%template(Simplex1D1D) Simplex<1, 1>;
%template(Simplex1D2D) Simplex<2, 1>;
%template(Simplex1D3D) Simplex<3, 1>;
%template(Simplex2D2D) Simplex<2, 2>;
%template(Simplex2D3D) Simplex<3, 2>;
%template(Simplex3D3D) Simplex<3, 3>;

%template(Vertex1D) Vertex<1>;
%template(Vertex1DList) std::vector<Vertex<1>>;
%template(Vertex2D) Vertex<2>;
%template(Vertex2DList) std::vector<Vertex<2>>;
%template(Vertex3D) Vertex<3>;
%template(Vertex3DList) std::vector<Vertex<3>>;
%template(Edge1D) Edge<1>;
%template(Edge1DList) std::vector<Edge<1>>;
%template(Edge2D) Edge<2>;
%template(Edge2DList) std::vector<Edge<2>>;
%template(Edge3D) Edge<3>;
%template(Edge3DList) std::vector<Edge<3>>;
%template(Face2D) Face<2>;
%template(Face2DList) std::vector<Face<2>>;
%template(Face3D) Face<3>;
%template(Face3DList) std::vector<Face<3>>;
%template(Cell3DList) std::vector<Cell>;

%template(FacetsEdge1D) std::array<Vertex<1>, 2>;
%template(FacetsEdge2D) std::array<Vertex<2>, 2>;
%template(FacetsEdge3D) std::array<Vertex<3>, 2>;
%template(FacetsFace2D) std::array<Edge<2>, 3>;
%template(FacetsFace3D) std::array<Edge<3>, 3>;
%template(FacetsCell3D) std::array<Face<3>, 4>;

%ignore getFacesOfEdge;

%template(Mesh0D) Mesh<0, 0>;
%template(Mesh10D) Mesh<1, 0>;
%template(Mesh20D) Mesh<2, 0>;
%template(Mesh30D) Mesh<3, 0>;

%template(Mesh11D) Mesh<1, 1>;
%template(Mesh21D) Mesh<2, 1>;
%template(Mesh21DList) std::vector<Mesh<2, 1>*>;
%template(Mesh2D) Mesh<2, 2>;

%template(Mesh31D) Mesh<3, 1>;
%template(Mesh32D) Mesh<3, 2>;
%extend Mesh<2, 1> {
   %template(getOrCreatePoint2DAttributeOnEdge) getOrCreateAttributeOnEdge<Eigen::Matrix<double, 2, 1>>;
   %template(getOrCreatePoint2DAttributeOnVertex) getOrCreateAttributeOnVertex<Eigen::Matrix<double, 2, 1>>;
   %template(getOrCreateIntAttributeOnVertex) getOrCreateAttributeOnVertex<int>;
   %template(getOrCreateUInt8AttributeOnVertex) getOrCreateAttributeOnVertex<uint8_t>;
}
%extend Mesh<2, 2> {
   %template(getOrCreateDoubleAttributeOnVertex) getOrCreateAttributeOnVertex<double>;
   %template(getOrCreateIntListAttributeOnVertex) getOrCreateAttributeOnVertex<std::vector<int>>;
   %template(getOrCreatePoint2DListAttributeOnFace) getOrCreateAttributeOnFace<std::vector<Eigen::Matrix<double, 2, 1>>>;
   %template(getOrCreatePoint2DAttributeOnFace) getOrCreateAttributeOnFace<Eigen::Matrix<double, 2, 1>>;
   %template(getOrCreatePoint2DAttributeOnEdge) getOrCreateAttributeOnEdge<Eigen::Matrix<double, 2, 1>>;
   %template(getOrCreatePoint2DListAttributeOnVertex) getOrCreateAttributeOnVertex<std::vector<Eigen::Matrix<double, 2, 1>>>;
};