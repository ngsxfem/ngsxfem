/*********************************************************************/
/* File:   SpaceTimeFE.cpp (based on mylittlengsolve)                */
/* Author: Janosch Preuss                                            */
/* Date:   June 2017                                                 */
/*********************************************************************/
#include <fem.hpp>
#include "SpaceTimeFE.hpp"
#include "../utils/ngsxstd.hpp"


namespace ngfem
{

   template <int D>
   SpaceTimeFE<D> :: SpaceTimeFE (ScalarFiniteElement<D>* s_FE, ScalarFiniteElement<1>* t_FE, bool aoverride_time, double atime)
    /*
      Call constructor for base class:
      number of dofs is (dofs in space) * (Dofs in time), maximal order is order
     */
      : ScalarFiniteElement<D> ((s_FE->GetNDof())*(t_FE->GetNDof()), s_FE->Order())
    {

        sFE = s_FE;
        tFE = t_FE;
        time = atime;
        override_time = aoverride_time;
    }




    template <int D>
    void SpaceTimeFE<D> :: CalcTimeShape (const IntegrationPoint & ip,
                                          BareSliceVector<> time_shape,
                                          int derivorder) const
    {
      IntegrationPoint z(override_time ? time : ip.Weight());

      if(!IsSpaceTimeIntegrationPoint(ip))
        throw Exception("SpaceTimeFE :: CalcShape called with a mere space IR");

          //Vector<> time_shape(tFE->GetNDof());
      FlatMatrix<> t_shape(tFE->GetNDof(),1,&time_shape(0));
      if (derivorder == 0)
        tFE->CalcShape(z,time_shape);
      else if (derivorder == 1)
        tFE->CalcDShape(z,t_shape);
      else if (derivorder == 2)
        tFE->CalcDDShape(z,t_shape);
      else 
        throw Exception("SpaceTimeFE :: CalcTimeShape called with an invalid time derivative order");
    }


   template <int D>
    void SpaceTimeFE<D> :: CalcShape (const IntegrationPoint & ip,
                                    BareSliceVector<> shape) const
    {
      Vector<> time_shape(tFE->GetNDof());
      Vector<> space_shape(sFE->GetNDof());

      GenericCalcShape(ip, shape,
        [&](const IntegrationPoint& ipx, BareSliceVector<> x_shape) { sFE->CalcShape(ipx,x_shape); },
        space_shape, 1,
        [&](const IntegrationPoint& ipt, BareSliceVector<> t_shape) { tFE->CalcShape(ipt,t_shape); },
        time_shape
      );
     }



    template <int D>
    void SpaceTimeFE<D> :: CalcDShape (const IntegrationPoint & ip,
                                    BareSliceMatrix<> dshape) const

    {

      Vector<> time_shape(tFE->GetNDof());
      Matrix<> space_shape(sFE->GetNDof(),D);

      GenericCalcShape(ip, dshape,
        [&](const IntegrationPoint& ipx, BareSliceMatrix<> x_shape) { sFE->CalcDShape(ipx,x_shape); },
        space_shape, D,
        [&](const IntegrationPoint& ipt, BareSliceVector<> t_shape) { tFE->CalcShape(ipt,t_shape); },
        time_shape
      );  
      
    }


    template <int D>
    void SpaceTimeFE<D> :: CalcMappedDDShape (const BaseMappedIntegrationPoint & mip, 
                                      BareSliceMatrix<> ddshape) const

    {

      Vector<> time_shape(tFE->GetNDof());
      Matrix<> space_shape(sFE->GetNDof(),D*D);

      GenericCalcMappedShape(mip, ddshape,
        [&](const BaseMappedIntegrationPoint& mipx, BareSliceMatrix<> x_shape) { sFE->CalcMappedDDShape(mipx,x_shape); },
        space_shape, D*D,
        [&](const IntegrationPoint& ipt, BareSliceVector<> t_shape) { tFE->CalcShape(ipt,t_shape); },
        time_shape
      );  

    }


    // for time derivatives

    template <int D>
    void SpaceTimeFE<D> :: CalcDtShape (const IntegrationPoint & ip,
                                     BareSliceVector<> dshape) const

    {
      Matrix<> time_dshape(tFE->GetNDof(),1);
      Vector<> space_shape(sFE->GetNDof());
      GenericCalcShape(ip, dshape,
        [&](const IntegrationPoint& ipx, BareSliceVector<> x_shape) { sFE->CalcShape(ipx,x_shape); },
        space_shape, 1,
        [&](const IntegrationPoint& ipt, BareSliceMatrix<> t_shape) { tFE->CalcDShape(ipt,t_shape); },
        time_dshape
      );
      return;
    }

    template <int D>
    void SpaceTimeFE<D> :: CalcDDtShape (const IntegrationPoint & ip,
                                     BareSliceVector<> ddshape) const

    {
      Matrix<> time_dshape(tFE->GetNDof(),1);
      Vector<> space_shape(sFE->GetNDof());
      GenericCalcShape(ip, ddshape,
        [&](const IntegrationPoint& ipx, BareSliceVector<> x_shape) { sFE->CalcShape(ipx,x_shape); },
        space_shape, 1,
        [&](const IntegrationPoint& ipt, BareSliceMatrix<> t_shape) { tFE->CalcDDShape(ipt,t_shape); },
        time_dshape
      );
      return;
    }

   void LagrangePolyHornerCalc::CalcNewtonBasisCoeffs() {
        NewtonBasisCoeffs.SetSize(nodes.Size(), nodes.Size());
        for(int i=0; i<nodes.Size(); i++){
            Matrix<double> p(nodes.Size(),nodes.Size());
            for(int j=0; j<nodes.Size(); j++) p(j,0) = 0.;
            p(i,0) = 1;
            for(int k=1; k<nodes.Size(); k++){
                for(int ii=0; ii<nodes.Size()-k; ii++){
                    p(k+ii,k) = ((p(k+ii,k-1) - p(k+ii-1,k-1))/(nodes[k+ii] - nodes[ii]));
                }
            }
            for(int j=0; j<nodes.Size(); j++) NewtonBasisCoeffs(j,i) = p(j,j);
        }
   }

   void LagrangePolyHornerCalc::SetUpChilds(){
        for(int k=0; k<nodes.Size(); k++){
            Array<double> nodes_minus(nodes);
            nodes_minus.RemoveElement(k);
            LagrangePolyHornerCalc child(nodes_minus, false);
            my_childs.Append(child);
        }
    }

    double LagrangePolyHornerCalc::Lagrange_Pol_Horner (double x, int i) const {
        Array<double> b (nodes.Size());
        b[nodes.Size()-1] = NewtonBasisCoeffs(nodes.Size()-1,i);
        for(int j=nodes.Size()-2; j>=0; j--) b[j] = b[j+1]*(x - nodes[j]) + NewtonBasisCoeffs(j,i);

        return b[0];
    }
    double LagrangePolyHornerCalc::Lagrange_Pol_D_Horner(double x, int i) const {
        if(my_childs.Size() == 0) throw Exception("LagrangePolyHornerCalc::Lagrange_Pol_D_Horner was called although instance was created in non-deriv mode");
        double sum = 0;
        for(int k=0; k<nodes.Size(); k++){
            if(k != i){
                sum += 1./(nodes[i] - nodes[k])*my_childs[k].Lagrange_Pol_Horner(x, (i > k ? i-1 : i));
            }
        }
        return sum;
    }

    NodalTimeFE :: NodalTimeFE (int order, bool askip_first_nodes, bool aonly_first_nodes, int ndof_first_node)
        : ScalarFiniteElement<1> (askip_first_nodes ? order + 1 - ndof_first_node          // skip_first_nodes
                                                    : (aonly_first_nodes ? ndof_first_node // only_first_nodes
                                                                         : order + 1),     // neither
                                  order), 
        skip_first_nodes(askip_first_nodes), only_first_nodes(aonly_first_nodes)
      {
         k_t = order;
         CalcInterpolationPoints ();

         if(order >= 5) do_horner_eval = true;

         if(do_horner_eval){
             LagrangePolyHornerCalc HornerLP2(nodes, true);
             HornerLP = HornerLP2;
         }
      }

      void NodalTimeFE :: CalcShape (const IntegrationPoint & ip,
                                     BareSliceVector<> shape) const
      {
         AutoDiff<1> adx (ip(0), 0);
         int begin = skip_first_nodes ? 1 : 0;
         int end = only_first_nodes ? 1 : ndof+begin;
         int cnt = 0;
         for(int i = begin; i < end; i++) {
             shape(cnt++) = (do_horner_eval ? HornerLP.Lagrange_Pol_Horner(ip(0),i) : Lagrange_Pol (adx, i).Value()) ;
         }
      }


      void NodalTimeFE :: CalcDShape (const IntegrationPoint & ip,
                                      BareSliceMatrix<> dshape) const
      {
         AutoDiff<1> adx (ip(0), 0);
         int begin = skip_first_nodes ? 1 : 0;
         int end = only_first_nodes ? 1 : ndof+begin;
         int cnt = 0;
         for(int i = begin; i < end; i++) {
             dshape(cnt++,0) = (do_horner_eval ? HornerLP.Lagrange_Pol_D_Horner(ip(0),i) : Lagrange_Pol(adx, i).DValue(0));
          }
      }

      Vector<double> CalcLobattoPointsRec(int order) {
          int N = order;
          Vector<double> x(N+1); Vector<double> x_old(N+1);
          for(int j=0; j<N+1; j++) {
              x[j] = cos(M_PI*((double)j)/N);
              x_old[j] = 0.0;
          }

          Matrix<double> P(N+1, N+1);
          int its = 0;
          while( (its < 100) && (L2Norm (x -x_old) > 1e-20) ) {
              x_old = x;
              for(int j=0; j<N+1; j++){
                  P(j, 0) = 1; P(j,1) = x[j];
              }
              for(int k=1; k<N; k++){
                  for(int j=0; j<N+1; j++) P(j,k+1)=( (2.*k+1)*(x[j] * P(j,k)) -((double)k)*P(j,k-1) )/((double)k+1.);
              }
              for(int j=0;j<N+1; j++) x[j]=x_old[j] -( x[j] * P(j,N) - P(j,N-1) )/( N*P(j,N) );
              its ++;
          }
          for(int j=0; j<N+1; j++) x[j] = 0.5*(1-x[j]);

          //This can be used for comparison with the hard-coded values for orders up to 5
          //NodalTimeFE fe(order, false, false, 1);
          //for(int j=0; j<N+1; j++) cout << x[j] - fe.GetNodes()[j] << endl;
          return x;
      }

      void NodalTimeFE :: CalcInterpolationPoints ()
      {
         nodes.SetSize(order+1);
         switch (order)
         {
          // Gauss-Lobatto integration points (Spectral FE)
          // The maximum order implemented here is mentioned in python_spacetime.cpp
          // in a documentation string. Please update that after inserting higher orders here.
          case 0 : nodes[0] = 0.0;  break;
          case 1 : nodes[0] = 0.0; nodes[1] = 1.0;  break;
          case 2 : nodes[0] = 0.0; nodes[1] = 0.5; nodes[2] = 1.0;  break;
          case 3 : nodes[0] = 0.0; nodes[1] = 0.5*(1.0-1.0/sqrt(5.0));
                   nodes[2] = 0.5*(1.0+1.0/sqrt(5.0)); nodes[3] = 1.0;  break;
          case 4 : nodes[0] = 0.0; nodes[1] = 0.5*(1.0-sqrt(3.0/7.0)); nodes[2] = 0.5;
                   nodes[3] = 0.5*(1.0+sqrt(3.0/7.0)); nodes[4] = 1.0;  break;
          case 5 : nodes[0] = 0.0;
                   nodes[1] = 0.5*(1.0 - sqrt(1.0/3.0 + 2.0*sqrt(7.0)/21.0));
                   nodes[2] = 0.5*(1.0 - sqrt(1.0/3.0 - 2.0*sqrt(7.0)/21.0));
                   nodes[3] = 0.5*(1.0 + sqrt(1.0/3.0 - 2.0*sqrt(7.0)/21.0));
                   nodes[4] = 0.5*(1.0 + sqrt(1.0/3.0 + 2.0*sqrt(7.0)/21.0));
                   nodes[5] = 1.0;  break;
          default :
                  auto nds = CalcLobattoPointsRec(order);
                  for(int j=0;j<order+1;j++) nodes[j] = nds[j];
         }
      }

    GCC3FE :: GCC3FE (bool askip_first_nodes, bool aonly_first_nodes)
        : NodalTimeFE (3, askip_first_nodes, aonly_first_nodes, 2)
      {
         ;
      }


      void GCC3FE :: CalcShape (const IntegrationPoint & ip,
                                     BareSliceVector<> shape) const
      {
         AutoDiff<1> x (ip(0), 0);
         int begin = skip_first_nodes ? 2 : 0;
         int end = only_first_nodes ? 1 : ndof+begin;
         int cnt = 0;
         if (!skip_first_nodes)
         {
           shape(cnt++) = ((1-x)*(1-x)*(1+2*x)).Value();
           shape(cnt++) = ((1-x)*(1-x)*x).Value();
         }
         
         if (!only_first_nodes)
         {
           shape(cnt++) = (x*x*(3-2*x)).Value();
           shape(cnt++) = (x*x*(x-1)).Value();
         }
      }


      void GCC3FE :: CalcDShape (const IntegrationPoint & ip,
                                      BareSliceMatrix<> dshape) const
      {
         AutoDiff<1> x (ip(0), 0);
         int begin = skip_first_nodes ? 2 : 0;
         int end = only_first_nodes ? 1 : ndof+begin;
         int cnt = 0;
         if (!skip_first_nodes)
         {
           dshape(cnt++,0) = ((1-x)*(1-x)*(1+2*x)).DValue(0);
           dshape(cnt++,0) = ((1-x)*(1-x)*x).DValue(0);
         }
         
         if (!only_first_nodes)
         {
           dshape(cnt++,0) = (x*x*(3-2*x)).DValue(0);
           dshape(cnt++,0) = (x*x*(x-1)).DValue(0);
         }
      }
      
      void GCC3FE :: CalcInterpolationPoints ()
      {
         nodes.SetSize(4);
        nodes[0] = 0.0; nodes[1] = 0.0; nodes[2] = 1.0; nodes[3] = 1.0;
      }

  
      template class SpaceTimeFE<0>;
      template class SpaceTimeFE<1>;
      template class SpaceTimeFE<2>;
      template class SpaceTimeFE<3>;

}
