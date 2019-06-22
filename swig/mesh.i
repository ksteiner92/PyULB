//#pragma SWIG nowarn=315
#define decltype(auto) // \
void ignore_me();
%ignore ignore_me;
%{
#include "mesh.h"
#include "meshing.h"

void ignore_me() {}
%}


%include <std_array.i>
%include <std_vector.i>
%include <stdint.i>
%include "attributes.i"
%include "mesh.h"
%include "meshing.h"

%extend Simplex
{
   int __getitem__(std::size_t i)
   {
      return (*$self)[i];
   }
};

%extend MeshElement
{
    int __getitem__(std::size_t i)
    {
       return (*$self)[i];
    }
};

%ignore Simplex<1, 0>::getFacet;
%ignore Simplex<2, 0>::getFacet;
%ignore Simplex<3, 0>::getFacet;

%define DEFINE_SIMPLEX(Dim, SimplexDim, ElementName)
%{
   namespace swig
   {
      template <> struct traits<Simplex<Dim, SimplexDim>>
      {
         typedef pointer_category category;
         static const char* type_name()
         {
            return  #ElementName;
         }
      };
   }
%}
%template(ElementName ##Dim ##D) Simplex<Dim, SimplexDim>;
%enddef


DEFINE_SIMPLEX(1, 0, Vertex)
DEFINE_SIMPLEX(2, 0, Vertex)
DEFINE_SIMPLEX(3, 0, Vertex)
DEFINE_SIMPLEX(1, 1, Edge)
DEFINE_SIMPLEX(2, 1, Edge)
DEFINE_SIMPLEX(3, 1, Edge)
DEFINE_SIMPLEX(2, 2, Face)
DEFINE_SIMPLEX(3, 2, Face)
DEFINE_SIMPLEX(3, 3, Cell)

%template(FacetsEdge1D) std::array<Simplex<1, 0>*, 2>;
%template(FacetsEdge2D) std::array<Simplex<2, 0>*, 2>;
%template(FacetsEdge3D) std::array<Simplex<3, 0>*, 2>;
%template(FacetsFace2D) std::array<Simplex<2, 1>*, 3>;
%template(FacetsFace3D) std::array<Simplex<3, 1>*, 3>;
%template(FacetsCell3D) std::array<Simplex<3, 2>*, 4>;

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
%template(Mesh3D) Mesh<3, 3>;

%include "meshattr.i"

%template(Delaunay2D) Meshing::generate<2, 2>;