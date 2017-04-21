#ifndef FILE_MYELEMENT_HPP
#define FILE_MYELEMENT_HPP

// based on:

/*********************************************************************/
/* File:   myElement.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/




namespace ngfem
{




    class SpaceTimeFE : public ScalarFiniteElement<2>
   {
        ScalarFiniteElement<2>* sFE = nullptr;
        ScalarFiniteElement<1>* tFE = nullptr;



    public:
      // constructor
      SpaceTimeFE (int order,ScalarFiniteElement<2>* s_FE,ScalarFiniteElement<1>*t_FE );

      virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

      /*
        Calculate the vector of shape functions in the point ip.
        ip is given in the reference element.
       */
      virtual void CalcShape (const IntegrationPoint & ip,
                              BareSliceVector<> shape) const;

      // for time derivatives

      virtual void CalcDtShape (const IntegrationPoint & ip,
                               SliceMatrix<> dshape) const;
      /*
        Calculate the matrix of derivatives of the shape functions in the point ip.
        dshape is an 3 by 2 matrix in our case.
       */
      virtual void CalcDShape (const IntegrationPoint & ip,
                               SliceMatrix<> dshape) const;

      // there are some more functions to bring in ...
      using ScalarFiniteElement<2>::CalcShape;
      using ScalarFiniteElement<2>::CalcDShape;
    };

 }


#endif

