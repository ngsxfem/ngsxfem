/// from ngxfem
#include "../xfem/cutinfo.hpp"
using namespace ngsolve;
namespace ngcomp
{
  CutInformation::CutInformation (shared_ptr<MeshAccess> ama)
    : ma(ama), cut_ratio_of_node(4), elems_of_domain_type(3), facets_of_domain_type(3)
  {
    for (auto dt : {NEG,POS,IF})
      elems_of_domain_type[dt] = make_shared<BitArray>(ma->GetNE(VOL));
    elems_of_domain_type[NEG]->Set();
    elems_of_domain_type[POS]->Clear();
    elems_of_domain_type[IF]->Clear();
    for (auto dt : {NEG,POS,IF})
      facets_of_domain_type[dt] = make_shared<BitArray>(ma->GetNFacets());
    facets_of_domain_type[NEG]->Set();
    facets_of_domain_type[POS]->Clear();
    facets_of_domain_type[IF]->Clear();
  }

  void CutInformation::Update(shared_ptr<CoefficientFunction> lset)
  {
    // NODE_TYPE: NT_VERTEX
    int nv = ma->GetNV();
    cut_ratio_of_node[NT_VERTEX] = make_shared<VVector<double>>(nv);
    // ...

    // NODE_TYPE: NT_EDGE
    int nedges = ma->GetNEdges();
    cut_ratio_of_node[NT_EDGE] = make_shared<VVector<double>>(nedges);
    // ...

    // NODE_TYPE: NT_FACE
    int nfaces = ma->GetNFaces();
    cut_ratio_of_node[NT_FACE] = make_shared<VVector<double>>(nfaces);
    // ...

    // NODE_TYPE: NT_CELL
    int ncells = ma->GetDimension() == 3 ? ma->GetNE() : 0;
    cut_ratio_of_node[NT_CELL] = make_shared<VVector<double>>(ncells);
    // ...
  }
}
