/*********************************************************************/
/* File:   myElement.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*

SpaceTimeFE

*/


#include <fem.hpp>
#include "myElement.hpp"


namespace ngfem
{


   SpaceTimeFE :: SpaceTimeFE (int order, ScalarFiniteElement<2>* s_FE, ScalarFiniteElement<1>* t_FE, double atime)
    /*
      Call constructor for base class:
      number of dofs is (dofs in space) * (Dofs in time), maximal order is order
     */
      : ScalarFiniteElement<2> ((s_FE->GetNDof())*(t_FE->GetNDof()), order)
    {

        sFE = s_FE;
        tFE = t_FE;
        time = atime;

    }

    void SpaceTimeFE :: CalcShape (const IntegrationPoint & ip,
                                    BareSliceVector<> shape) const
    {

       if (tFE->GetNDof() == 1)
          sFE->CalcShape(ip,shape);
       else {

            Vector<> time_shape(tFE->GetNDof());
            IntegrationPoint z(time);//ip(2));
            tFE->CalcShape(z,time_shape);

            Vector<> space_shape(sFE->GetNDof());
            sFE->CalcShape(ip,space_shape);

            // define shape functions
            int ii = 0;
            for(int j = 0; j < tFE->GetNDof(); j++) {
               for(int i=0; i< sFE->GetNDof() ; i++) {
                   shape(ii++) = space_shape(i)*time_shape(j);
               }
             }
       }
     }

    void SpaceTimeFE :: CalcDShape (const IntegrationPoint & ip,
                                     SliceMatrix<> dshape) const

    {
      // matrix of derivatives:

         if (tFE->GetNDof() == 1)
            sFE->CalcDShape(ip,dshape);
         else {

            Vector<> time_shape(tFE->GetNDof());
            IntegrationPoint z(time);
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
        IntegrationPoint z(time);
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


}
