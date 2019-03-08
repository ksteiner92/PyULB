%{
#include "poisson.h"
%}

%include <std_vector.i>
%include <stdint.i>
%include "mesh.i"
%include "poisson.h"

typedef long unsigned int size_t;
namespace std {
using ::size_t;
}
using std::size_t;

%template(VectorSizeT) std::vector<std::size_t>;

%template(Poisson1D) Poisson<1, 1>;
%template(Poisson2D1D) Poisson<2, 1>;
%template(Poisson2D) Poisson<2, 2>;
%template(Poisson3D1D) Poisson<3, 1>;
%template(Poisson3D2D) Poisson<3, 2>;
%template(Poisson3D) Poisson<3, 3>;