/*********************************************************************/
/* File:   SpaceTimeFE.cpp (based on mylittlengsolve)                */
/* Author: Janosch Preuss                                            */
/* Date:   June 2017                                                 */
/*********************************************************************/

#include <fem.hpp>
#include "SpaceTimeFE.hpp"


namespace ngfem
{


   SpaceTimeFE :: SpaceTimeFE (ScalarFiniteElement<2>* s_FE, ScalarFiniteElement<1>* t_FE, bool aoverride_time, double atime)
    /*
      Call constructor for base class:
      number of dofs is (dofs in space) * (Dofs in time), maximal order is order
     */
      : ScalarFiniteElement<2> ((s_FE->GetNDof())*(t_FE->GetNDof()), s_FE->Order())
    {

        sFE = s_FE;
        tFE = t_FE;
        time = atime;
        override_time = aoverride_time;
    }

    void SpaceTimeFE :: CalcShape (const IntegrationPoint & ip,
                                    BareSliceVector<> shape) const
    {

       if (tFE->GetNDof() == 1)
          sFE->CalcShape(ip,shape);
       else
       {

            Vector<> time_shape(tFE->GetNDof());
            IntegrationPoint z(override_time ? time : ip(2));
            tFE->CalcShape(z,time_shape);

            Vector<> space_shape(sFE->GetNDof());
            sFE->CalcShape(ip,space_shape);

            // define shape functions
            int ii = 0;
            for(int j=0; j<tFE->GetNDof(); j++) {
              for(int i=0; i<sFE->GetNDof(); i++) {
                shape(ii++) = space_shape(i)*time_shape(j);
              }
            }
       }
     }

    void SpaceTimeFE :: CalcDShape (const IntegrationPoint & ip,
                                    BareSliceMatrix<> dshape) const

    {
      // matrix of derivatives:

         if (tFE->GetNDof() == 1)
            sFE->CalcDShape(ip,dshape);
         else {

            Vector<> time_shape(tFE->GetNDof());
            IntegrationPoint z(override_time ? time : ip(2));
            tFE->CalcShape(z,time_shape);

            Matrix<double> space_dshape(sFE->GetNDof(),2);
            sFE->CalcDShape(ip,space_dshape);


            int ii = 0;
            for(int j = 0; j < tFE->GetNDof(); j++) {
                for(int i=0; i< sFE->GetNDof(); i++) {
                    dshape(ii,0) = space_dshape(i,0)*time_shape(j);
                    dshape(ii,1) = space_dshape(i,1)*time_shape(j);
                    ii++;
                }
            }
         }

    }

    // for time derivatives

    void SpaceTimeFE :: CalcDtShape (const IntegrationPoint & ip,
                                     BareSliceVector<> dshape) const

    {
        // matrix of derivatives:

           Matrix<double> time_dshape(tFE->GetNDof(),1);
           IntegrationPoint z(override_time ? time : ip(2));
           tFE->CalcDShape(z,time_dshape);

           Vector<> space_shape(sFE->GetNDof());
           sFE->CalcShape(ip,space_shape);

           int ii = 0;
           for(int j = 0; j < tFE->GetNDof(); j++) {
              for(int i=0; i< sFE->GetNDof(); i++) {
                 dshape(ii++) = space_shape(i)*time_dshape(j,0);
              }
           }

    }


    NodalTimeFE :: NodalTimeFE (int order)
        : ScalarFiniteElement<1> (order+1, order)
      {
       k_t = order;
      }


      void NodalTimeFE :: CalcShape (const IntegrationPoint & ip,
                                     BareSliceVector<> shape) const
      {
         Vector<double> intp_pts(k_t+1);
         GetIntpPts (intp_pts);
         AutoDiff<1> adx (ip(0), 0);
         for(int i = 0; i < k_t+1; i++) {
             shape(i) = Lagrange_Pol (adx, intp_pts , i).Value() ;
         }
      }


      void NodalTimeFE :: CalcDShape (const IntegrationPoint & ip,
                                      BareSliceMatrix<> dshape) const
      {
         Vector<double> intp_pts(k_t+1);
         GetIntpPts (intp_pts);
         AutoDiff<1> adx (ip(0), 0);
         for(int i = 0; i < k_t+1; i++) {
             dshape(i,0) = Lagrange_Pol(adx, intp_pts , i).DValue(0);
          }
      }

      void NodalTimeFE :: GetIntpPts (Vector<>& intp_pts) const
      {
         switch (intp_pts.Size())
         {
          // Gauss-Lobatto integration points (Spectral FE)
          case 1 : intp_pts(0) = 0.0;  break;
          case 2 : intp_pts(0) = 0.0; intp_pts(1) = 1.0;  break;
          case 3 : intp_pts(0) = 0.0; intp_pts(1) = 0.5; intp_pts(2) = 1.0;  break;
          case 4 : intp_pts(0) = 0.0; intp_pts(1) = 0.5*(1.0-1.0/sqrt(5.0));
                   intp_pts(2) = 0.5*(1.0+1.0/sqrt(5.0)); intp_pts(3) = 1.0;  break;
          case 5 : intp_pts(0) = 0.0; intp_pts(1) = 0.5*(1.0-sqrt(3.0/7.0)); intp_pts(2) = 0.5;
                   intp_pts(3) = 0.5*(1.0+sqrt(3.0/7.0)); intp_pts(4) = 1.0;  break;
          case 6 : intp_pts(0) = 0.0;
                   intp_pts(1) = 0.5*(1.0 - sqrt(1.0/3.0 + 2.0*sqrt(7.0)/21.0));
                   intp_pts(2) = 0.5*(1.0 - sqrt(1.0/3.0 - 2.0*sqrt(7.0)/21.0));
                   intp_pts(3) = 0.5*(1.0 + sqrt(1.0/3.0 - 2.0*sqrt(7.0)/21.0));
                   intp_pts(4) = 0.5*(1.0 + sqrt(1.0/3.0 + 2.0*sqrt(7.0)/21.0));
                   intp_pts(5) = 1.0;  break;
          default : throw Exception("Requested TimeFE not implemented yet.");
         }
      }


}
