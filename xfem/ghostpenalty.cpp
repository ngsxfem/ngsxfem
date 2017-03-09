#define FILE_GHOSTPENALTY_CPP
#include "ghostpenalty.hpp"
#include <diffop_impl.hpp>

namespace ngfem
{

  void SetRefBaryCenter(ELEMENT_TYPE eltype, FlatVector<> & point)
  {
    switch (eltype)
    {
    case ET_POINT : point = 0.0; break;
    case ET_SEGM : point = 0.5; break;
    case ET_TRIG : point = 1.0/3.0; break;
    case ET_QUAD : point = 0.5; break;
    case ET_TET : point = 1.0/4.0; break;
    case ET_PYRAMID : point(0) = 1.0/2.0; point(1) = 1.0/2.0; point(2) = 1.0/4.0; break;
    case ET_PRISM : point(0) = 1.0/3.0; point(1) = 1.0/3.0; point(2) = 1.0/2.0; break;
    case ET_HEX : point = 1.0/2.0; break;
    };
  };

  template <int D, int difforder>
  void GhostPenaltyIntegrator<D,difforder>::CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                                                             const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                                                             const FiniteElement & volumefel2, int LocalFacetNr2,
                                                             const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                                                             FlatMatrix<double> & elmat,
                                                             LocalHeap & lh// ,
                                                             // BitArray* twice
                                                             ) const
  {
    static int timer = NgProfiler::CreateTimer ("GhostPenaltyIntegrator");

    if (LocalFacetNr2==-1) throw Exception("GhostPenaltyIntegrator: LocalFacetNr2==-1");
    if (D==3) throw Exception("scaling is not correct!");
    NgProfiler::RegionTimer reg (timer);

    const CompoundFiniteElement * dcfel1 =
      dynamic_cast<const CompoundFiniteElement*> (&volumefel1);
    const CompoundFiniteElement * dcfel2 =
      dynamic_cast<const CompoundFiniteElement*> (&volumefel2);
    if (dcfel1==NULL || dcfel2==NULL) {
      cout << "not compound!" << endl;
      cout << " leaving " << endl;
      return;
    }

    const CompoundFiniteElement & cfel1 = *dcfel1;
    const CompoundFiniteElement & cfel2 = *dcfel2;

    const XFiniteElement * xfe1 = NULL;
    const XDummyFE * dummfe1 = NULL;
    const ScalarFiniteElement<D> * scafe1 = NULL;
    // const ScalarSpaceTimeFiniteElement<D> * stscafe1 = NULL;
    const void * stscafe1 = NULL;
    const XFiniteElement * xfe2 = NULL;
    const XDummyFE * dummfe2 = NULL;
    const ScalarFiniteElement<D> * scafe2 = NULL;
    // const ScalarSpaceTimeFiniteElement<D> * stscafe2 = NULL;
    const void * stscafe2 = NULL;

    for (int i = 0; i < cfel1.GetNComponents(); ++i)
    {
      if (xfe1==NULL)
        xfe1 = dynamic_cast<const XFiniteElement* >(&cfel1[i]);
      if (dummfe1==NULL)
        dummfe1 = dynamic_cast<const XDummyFE* >(&cfel1[i]);
      if (scafe1==NULL)
        scafe1 = dynamic_cast<const ScalarFiniteElement<D>* >(&cfel1[i]);
      // if (stscafe1==NULL)
      //   stscafe1 = dynamic_cast<const ScalarSpaceTimeFiniteElement<D>* >(&cfel1[i]);
    }

    for (int i = 0; i < cfel2.GetNComponents(); ++i)
    {
      if (xfe2==NULL)
        xfe2 = dynamic_cast<const XFiniteElement* >(&cfel2[i]);
      if (dummfe2==NULL)
        dummfe2 = dynamic_cast<const XDummyFE* >(&cfel2[i]);
      if (scafe2==NULL)
        scafe2 = dynamic_cast<const ScalarFiniteElement<D>* >(&cfel2[i]);
      // if (stscafe2==NULL)
      //   stscafe2 = dynamic_cast<const ScalarSpaceTimeFiniteElement<D>* >(&cfel2[i]);
    }

    const bool spacetime = (stscafe1 != NULL);
    elmat = 0.0;

    DOMAIN_TYPE facetdt = IF;

    if (!xfe1 || !xfe2) { // not a ghost edge
      // if (true)
      //   return; // only ghost penalty at the interface
      if (!(xfe1 || xfe2))
        return;
      else
        facetdt = xfe1 ? dummfe2->GetDomainType() : dummfe1->GetDomainType();
    }

    ELEMENT_TYPE eltype1 = volumefel1.ElementType();
    int nd1 = volumefel1.GetNDof();
    int ndof_x1 = xfe1 ? xfe1->GetNDof() : 0;
    // int ndof_sca1 = spacetime ? stscafe1->GetNDof() : scafe1->GetNDof();
    int ndof_sca1 = scafe1->GetNDof();
    FlatVector<> mat1_dudn(nd1, lh);
    FlatVector<> mat1_dudn_sca(ndof_sca1, &mat1_dudn(0));
    FlatVector<> mat1_dudn_x(ndof_x1, &mat1_dudn(ndof_sca1));
    FlatMatrixFixWidth<D> mat1_du_sca(ndof_sca1, lh);


    ELEMENT_TYPE eltype2 = volumefel2.ElementType();
    int nd2 = volumefel2.GetNDof();
    int ndof_x2 = xfe2 ? xfe2->GetNDof() : 0;
    // int ndof_sca2 = spacetime ? stscafe2->GetNDof() : scafe2->GetNDof();
    int ndof_sca2 = scafe2->GetNDof();
    FlatVector<> mat2_dudn(nd2, lh);
    FlatVector<> mat2_dudn_sca(ndof_sca2, &mat2_dudn(0));
    FlatVector<> mat2_dudn_x(ndof_x2, &mat2_dudn(ndof_sca2));
    FlatMatrixFixWidth<D> mat2_du_sca(ndof_sca2, lh);

    // int maxorder = spacetime ?
    //   max(stscafe1->OrderSpace(),stscafe2->OrderSpace())
    //   : max(scafe1->Order(),scafe2->Order());
    int maxorder = max(scafe1->Order(),scafe2->Order());

    if (maxorder==0) maxorder=1;

    int maxordertime = 0; //spacetime ? max(stscafe1->OrderTime(),stscafe2->OrderTime()) : 0;

    FlatMatrixFixWidth<1> bmat(nd1+nd2, lh);

    // const MasterElement & masterel1 = xfe1->GetMasterElement();
    // const MasterElement & masterel2 = xfe2->GetMasterElement();

    // Edge * commonedge = NULL;
    // for (int i = 0; i < masterel1.NBaseEdges(); ++i)
    //   for (int j = 0; j < masterel2.NBaseEdges(); ++j)
    //     if (masterel1.GetBaseEdge(i) == masterel2.GetBaseEdge(j)){
    //       commonedge = masterel1.GetBaseEdge(i);
    //     }

    Facet2ElementTrafo transform1(eltype1,ElVertices1);
    Facet2ElementTrafo transform2(eltype2,ElVertices2);

    const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);
    const NORMAL * normals2 = ElementTopology::GetNormals(eltype2);

    // HeapReset hr(lh);

    ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    FlatVector<double> eltype1_ref_barcenter(3,lh);
    SetRefBaryCenter(eltype1,eltype1_ref_barcenter);
    IntegrationPoint refbarycenterv(eltype1_ref_barcenter,0.0);
    MappedIntegrationPoint<D,D> barycenterv (refbarycenterv, eltrans1);
    Vec<D> pointv = barycenterv.GetPoint();
    FlatVector<double> eltype2_ref_barcenter(3,lh);
    SetRefBaryCenter(eltype2,eltype2_ref_barcenter);
    IntegrationPoint refbarycenterv2(eltype2_ref_barcenter,0.0);
    MappedIntegrationPoint<D,D> barycenterv2 (refbarycenterv2, eltrans2);
    Vec<D> pointv2 = barycenterv2.GetPoint();
    Vec<D> pointsdiff = pointv - pointv2;

    Vec<D> normal_ref1, normal_ref2;
    for (int i=0; i<D; i++) {
      normal_ref1(i) = normals1[LocalFacetNr1][i];
      normal_ref2(i) = normals2[LocalFacetNr2][i];
    }

    const int p = maxorder;

    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      if (facetdt==POS && dt==NEG) continue;
      if (facetdt==NEG && dt==POS) continue;

      const IntegrationRule & ir_facet =
        SelectIntegrationRule (etfacet, 2*p);

      const IntegrationRule & ir_time =
        SelectIntegrationRule (ET_SEGM, 2*maxordertime);

      // IntegrationRule ir_facet;
      // commonedge->FillSurfaceIntegrationRule(2*p,sign,ir_facet);
      bmat = 0.0;

      for (int k = 0; k < ir_time.GetNIP(); k++)
      {
        for (int l = 0; l < ir_facet.GetNIP(); l++)
        {
          IntegrationPoint ip1 = transform1(LocalFacetNr1, ir_facet[l]);
          MappedIntegrationPoint<D,D> sip1 (ip1, eltrans1);
          double lam = (dt == POS) ? coef_lam_pos->Evaluate(sip1) : coef_lam_neg->Evaluate(sip1);
          const double delta = coef_delta->Evaluate(sip1);
          Mat<D> inv_jac1 = sip1.GetJacobianInverse();
          double det1 = sip1.GetJacobiDet();

          Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;
          double len1 = L2Norm (normal1);
          normal1 /= len1;
          Vec<D> invjac_normal1 = inv_jac1 * normal1;
          if (difforder==1)
          {
            if (spacetime)
            {
              throw Exception("no space time supported anymore");
              // stscafe1->CalcDxShapeSpaceTime (sip1.IP(), ir_time[k](0), mat1_du_sca, lh);
              // mat1_dudn_sca = mat1_du_sca * invjac_normal1;
            }
            else
              mat1_dudn_sca = scafe1->GetDShape (sip1.IP(), lh) * invjac_normal1;
          }
          else
          {
            if (spacetime)
              throw Exception("no higher order ghost penalty for spacetime yet");
            // DO num diff in normal direction
            double eps = 1e-7;

            IntegrationPoint ipl(ip1);
            IntegrationPoint ipc(ip1);
            IntegrationPoint ipr(ip1);

            for (int d = 0; d < D; ++d)
              ipl(d) += eps * invjac_normal1(d);
            for (int d = 0; d < D; ++d)
              ipr(d) -= eps * invjac_normal1(d);

            // double len = L2Norm(invjac_normal1);

            FlatVector<> shapeleft = scafe1->GetShape (ipl, lh);
            FlatVector<> shapecenter = scafe1->GetShape (ipc, lh);
            FlatVector<> shaperight = scafe1->GetShape (ipr, lh);

            mat1_dudn_sca = shaperight - 2 * shapecenter + shapeleft;
            mat1_dudn_sca *= 1.0/(eps*eps);

          }
          if (xfe1)
          {
            mat1_dudn_x = mat1_dudn_sca;
            for (int i = 0; i < ndof_sca1; ++i)
              if (xfe1->GetSignsOfDof()[i]!=dt)
                mat1_dudn_x(i) = 0;
          }

          const double orthdist = abs(InnerProduct(pointsdiff,normal1));

          IntegrationPoint ip2 = transform2(LocalFacetNr2, ir_facet[l]);
          MappedIntegrationPoint<D,D> sip2 (ip2, eltrans2);
          Mat<D> inv_jac2 = sip2.GetJacobianInverse();
          double det2 = sip2.GetJacobiDet();
          Vec<D> normal2 = det2 * Trans (inv_jac2) * normal_ref2;
          double len2 = L2Norm (normal2);
          if(abs(len1-len2)>1e-6) {
            std::cout << "len :\t" << len1 << "\t=?=\t" << len2 << std::endl;
            throw Exception ("GhostPenaltyIntegrator: len1!=len2");
          }
          normal2 /= len2;
          Vec<D> invjac_normal2;;
          invjac_normal2 = inv_jac2 * normal2;

          if (difforder == 1)
          {
            if (spacetime)
            {
              throw Exception("no space time supported anymore");
              // stscafe2->CalcDxShapeSpaceTime (sip2.IP(), ir_time[k](0), mat2_du_sca, lh);
              // mat2_dudn_sca = mat2_du_sca * invjac_normal2;
            }
            else
              mat2_dudn_sca = scafe2->GetDShape (sip2.IP(), lh) * invjac_normal2;
          }
          else
          {
            if (spacetime)
              throw Exception("no higher order ghost penalty for spacetime yet");
            // DO num diff in normal direction
            double eps = 1e-7;

            IntegrationPoint ipl(ip2);
            IntegrationPoint ipc(ip2);
            IntegrationPoint ipr(ip2);

            for (int d = 0; d < D; ++d)
              ipl(d) -= eps * invjac_normal2(d);
            for (int d = 0; d < D; ++d)
              ipr(d) += eps * invjac_normal2(d);

            // double len = L2Norm(invjac_normal2);

            FlatVector<> shapeleft = scafe2->GetShape (ipl, lh);
            FlatVector<> shapecenter = scafe2->GetShape (ipc, lh);
            FlatVector<> shaperight = scafe2->GetShape (ipr, lh);

            mat2_dudn_sca = shaperight - 2 * shapecenter + shapeleft;
            mat2_dudn_sca *= 1.0/(eps*eps);

          }

          if (xfe2)
          {
            mat2_dudn_x = mat2_dudn_sca;
            for (int i = 0; i < ndof_sca2; ++i)
              if (xfe2->GetSignsOfDof()[i]!=dt)
                mat2_dudn_x(i) = 0;
          }

          bmat = 0.0;
          bmat.Col(0).Range(  0,    nd1) = mat1_dudn;
          bmat.Col(0).Range(nd1,nd1+nd2) = mat2_dudn;

          double pen = std::pow(orthdist, 2*difforder-1);
          elmat += pen * delta * lam * len1 * tau * ir_time[k].Weight() * ir_facet[l].Weight() * bmat * Trans (bmat);
        }
      }
    }
  }

  template class GhostPenaltyIntegrator<2,1>;
  template class GhostPenaltyIntegrator<3,1>;

  static RegisterBilinearFormIntegrator<GhostPenaltyIntegrator<2,1> > init_gp_2d ("lo_ghostpenalty", 2, 3);
  static RegisterBilinearFormIntegrator<GhostPenaltyIntegrator<3,1> > init_gp_3d ("lo_ghostpenalty", 3, 3);
  static RegisterBilinearFormIntegrator<GhostPenaltyIntegrator<2,2> > init_gp2_2d ("sec_ghostpenalty", 2, 3);
  static RegisterBilinearFormIntegrator<GhostPenaltyIntegrator<3,2> > init_gp2_3d ("sec_ghostpenalty", 3, 3);
  static RegisterBilinearFormIntegrator<GhostPenaltyIntegrator<2,1> > init_gp_2d_st ("stx_lo_ghostpenalty", 2, 5);
  static RegisterBilinearFormIntegrator<GhostPenaltyIntegrator<3,1> > init_gp_3d_st ("stx_lo_ghostpenalty", 3, 5);



  class CentralFDStencils
  {
    static const int max_accuracy_half = 8; // order 2,4,6,8,..,16
    static const int max_order = 10; // order 1,2,3,4,..,10
    Table<double> * stencils;
  public:
    static CentralFDStencils & Instance()
    {
      static CentralFDStencils myInstance;
      return myInstance;
    }
    static const FlatVector<> Get(int order, int accuracy)
    {
      CentralFDStencils & instance = Instance();
      int row = order == 0 ? 0 : 1+(order-1)*max_accuracy_half+(accuracy-1)/2;
      const FlatArray<double> fa ((*(instance.stencils))[row]);
      FlatVector<> coeffs(fa.Size(),&fa[0]);
      return coeffs;
    }
    static double GetOptimalEps(int order, int accuracy)
    {
      const double EPS = std::numeric_limits<double>::epsilon();
      const int stencil_width = order == 0 ? 0 : (accuracy+1)/2 + int(order-0.5)/2;
      const double eps = order == 0 ? 1.0 : std::pow((2 * stencil_width + 1)*EPS,1.0/(accuracy + order)); // + 0.25;
      static bool first0 = true;
      if (first0 && order==0)
      {

        cout << " order, eps = " << order << ", " << eps << endl;
        first0 = false;
        // getchar();
      }
      static bool first1 = true;
      if (first1 && order==1)
      {

        cout << " order, eps = " << order << ", " << eps << endl;
        first1 = false;
        // getchar();
      }
      static bool first2 = true;
      if (first2 && order==2)
      {
        cout << " order, eps = " << order << ", " << eps << endl;
        first2 = false;
        // getchar();
      }
      static bool first3 = true;
      if (first3 && order==3)
      {
        cout << " order, eps = " << order << ", " << eps << endl;
        first3 = false;
        // getchar();
      }
      static bool first4 = true;
      if (first4 && order==4)
      {
        cout << " order, eps = " << order << ", " << eps << endl;
        first4 = false;
        // getchar();
      }
      static bool first5 = true;
      if (first5 && order==5)
      {
        cout << " order, eps = " << order << ", " << eps << endl;
        first5 = false;
        // getchar();
      }
      return eps;

    }

  protected:
    CentralFDStencils()
    {
      const int rows = max_accuracy_half*max_order+1;

      Array<int> cnt(rows);

      cnt[0] = 1;
      for (int i = 0; i < max_order; ++i)
        for (int j = 0; j < max_accuracy_half; ++j)
        {
          const int row = i * max_accuracy_half + j + 1;
          const int order = i+1;
          // const int accuracy = 2*(j+1);
          const int stencil_width = j+1 + int(order-0.5)/2;
          cnt[row] = 2 * stencil_width +1;
        }

      stencils = new Table<double>(cnt);

      (*stencils)[0][0] = 1.0;
      for (int i = 0; i < max_order; ++i)
        for (int j = 0; j < max_accuracy_half; ++j)
        {
          const int row = i * max_accuracy_half + j + 1;
          const int order = i+1;
          // const int accuracy = 2*(j+1);
          const int stencil_width = j+1 + int(order-0.5)/2;
          const int size = 2 * stencil_width + 1;

          Vector<> stencil(size);
          {
            Matrix<> A(size);
            Vector<> f(size);
            double inv_factorial=1.0;
            for (int k = 0; k < size; k++)
            {
              if (k>0)
                inv_factorial/=k;
              for (int l = 0; l < size; l++)
                A(k,l) = std::pow(l-stencil_width,k) * inv_factorial;
              f(k) = k==order ? 1.0 : 0;
            }
            stencil = Inv(A) * f;
          }

          for (int col = 0; col < size; ++col)
            (*stencils)[row][col] = stencil(col);
        }
    }

    ~CentralFDStencils()
    {
      delete stencils;
    }
  };


  template <int D, int ORDER>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpDuDnk<D,ORDER>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                             MAT & mat, LocalHeap & lh)
  {
    const int FD_ACCURACY = 4;
    int version = 2;
    // mat.Row(0) = 0.0;
    // return;
    // FlatVector<> center(3,lh);
    // SetRefBaryCenter(bfel.ElementType(),center);
    // IntegrationPoint centerip (center(0),center(1),center(2),0.0);
    // MappedIntegrationPoint<D,D> mipc(centerip,mip.GetTransformation());

    // cout << mipc.GetPoint() << endl;

    if (version == 1)
    // not higher order accurate on curved meshes (!),
    // but more stable and efficient (one derivate less to evaluate by FD)
    {
      const ScalarFiniteElement<D> & scafe =
        dynamic_cast<const ScalarFiniteElement<D> & > (bfel);
      const int ndof = scafe.GetNDof();

      Vec<D> normal = static_cast<const DimMappedIntegrationPoint<D>&>(mip).GetNV();

      // cout << "normal: " << normal << endl;
      // Vec<D> normal; normal(0) = -1.0; normal(1) = 1.0;
      // normal /= L2Norm(normal);

      Vec<D> invjac_normal = mip.GetJacobianInverse() * normal;

      FlatVector<> fdstencil (CentralFDStencils::Get(ORDER-1,FD_ACCURACY));
      const double eps = CentralFDStencils::GetOptimalEps(ORDER-1,FD_ACCURACY);
      const int stencilpoints = fdstencil.Size();
      const int stencilwidth = (stencilpoints-1)/2;
      FlatMatrix<> dshapedn  (ndof, stencilpoints, lh);
      FlatVector<> dshapednk (ndof, lh);

      const double norminvjacn = L2Norm(invjac_normal);
      Vec<D> normal_dir_wrt_ref = (1.0/norminvjacn)  * invjac_normal;


      for (int i = 0; i < stencilpoints; ++i)
      {
        IntegrationPoint ip(mip.IP());
        for (int d = 0; d < D; ++d)
          ip(d) += (i-stencilwidth) * eps * normal_dir_wrt_ref(d);
        dshapedn.Col(i) = scafe.GetDShape (ip, lh) * invjac_normal;
      }
      dshapednk = dshapedn * fdstencil;
      const double eps_fac = std::pow(norminvjacn/eps,ORDER-1);
      mat.Row(0) = eps_fac * dshapednk;

      const int feorder = scafe.Order();
      if (false && feorder > 1)
      {
        cout << "ORDER = " << ORDER << endl;
        cout << mip.GetPoint() << endl;
        cout << "dshapedn" << endl;
        cout << dshapedn << endl;
        cout << "coeffs" << endl;
        cout << fdstencil << endl;
        cout << "dshapednk" << endl;
        cout << dshapednk << endl;
        getchar();
      }
    }
    else
    {
      // VERSION 2 (fix points on phys. domain + transform points(!) back)
      // higher order accurate also on curved meshes (!),
      // but takes all derivatives by FD
      const ScalarFiniteElement<D> & scafe =
        dynamic_cast<const ScalarFiniteElement<D> & > (bfel);
      const int ndof = scafe.GetNDof();

      Vec<D> normal = static_cast<const DimMappedIntegrationPoint<D>&>(mip).GetNV();

      // cout << "normal: " << normal << endl;
      // Vec<D> normal; normal(0) = -1.0; normal(1) = 1.0;
      // normal /= L2Norm(normal);

      Vec<D> invjac_normal = mip.GetJacobianInverse() * normal;

      const double h = D==2 ? sqrt(mip.GetJacobiDet()) : cbrt(mip.GetJacobiDet());

      FlatVector<> fdstencil (CentralFDStencils::Get(ORDER,FD_ACCURACY));
      const double eps = h * CentralFDStencils::GetOptimalEps(ORDER,FD_ACCURACY);

      const int stencilpoints = fdstencil.Size();
      const int stencilwidth = (stencilpoints-1)/2;

      FlatMatrix<> shape  (ndof, stencilpoints, lh);
      FlatVector<> dshapednk (ndof, lh);

      Vec<D> normal_dir_wrt_ref = invjac_normal;

      for (int i = 0; i < stencilpoints; ++i)
      {
        Vec<D> vec = mip.GetPoint();
        vec += (i-stencilwidth) * eps * normal;

        IntegrationPoint ip_x0(mip.IP());
        for (int d = 0; d < D; ++d)
          ip_x0(d) += (i-stencilwidth) * eps * normal_dir_wrt_ref(d);
        MappedIntegrationPoint<D,D> mip_x0(ip_x0,mip.GetTransformation());

        Vec<D> diff = vec - mip_x0.GetPoint();
        Vec<D> update = 0;
        int its = 0;
        while (L2Norm(diff) > 1e-8*h && its < 20)
        {
          MappedIntegrationPoint<D,D> mip_x0(ip_x0,mip.GetTransformation());
          diff = vec - mip_x0.GetPoint();
          // cout << "mip_x0.GetPoint()" << mip_x0.GetPoint() << endl;
          // cout << "diff" << diff << endl;
          update = mip_x0.GetJacobianInverse() * diff;
          for (int d = 0; d < D; ++d)
            ip_x0(d) += update(d);
          // getchar();
          its++;
        }
        if (its >= 50)
          cerr << "its >= 50 " << endl;
        shape.Col(i) = scafe.GetShape (ip_x0, lh);
      }

      dshapednk = shape * fdstencil;
      const double eps_fac = std::pow(1.0/eps,ORDER);
      mat.Row(0) = eps_fac * dshapednk;

      const int feorder = scafe.Order();
      if (false && feorder > 1)
      {
        cout << "ORDER = " << ORDER << endl;
        cout << mip.GetPoint() << endl;
        cout << "dshapedn" << endl;
        cout << shape << endl;
        cout << "coeffs" << endl;
        cout << fdstencil << endl;
        cout << "dshapednk" << endl;
        cout << dshapednk << endl;
        getchar();
      }
      // TODO: Add version 3:
      // VERSION 3 (fix points on phys. domain + transform points(!) back)
      // take du/dn (n the path-adjusted) normal direction and take FD from this.

    }


  }// generate matrix

  template class T_DifferentialOperator<DiffOpDuDnk<2,1>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,2>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,3>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,4>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,5>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,6>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,7>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,8>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,1>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,2>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,3>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,4>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,5>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,6>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,7>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,8>>;

}
