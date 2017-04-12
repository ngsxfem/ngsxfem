#ifndef FILE_MYELEMENT_HPP
#define FILE_MYELEMENT_HPP

/*********************************************************************/
/* File:   myElement.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*
  
My own simple first and second order triangular finite elements

*/


namespace ngfem
{

  /*
    A linear segment finite element is a scalar element in one dimension, 
    thus we derive it from the ScalarFiniteElement<1> base class:
   */

  class MyLinearSegm : public ScalarFiniteElement<1>
  {
  public:
    // constructor
    MyLinearSegm ();

    virtual ELEMENT_TYPE ElementType() const { return ET_SEGM; }
    
    /*
      Calculate the vector of shape functions in the point ip.
      ip is given in the reference element.
     */
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const;
  
    /*
      Calculate the matrix of derivatives of the shape functions in the point ip.
      dshape is an 2 by 1 matrix in our case.
     */
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             SliceMatrix<> dshape) const;

    // there are some more functions to bring in ...
    using ScalarFiniteElement<1>::CalcShape;    
    using ScalarFiniteElement<1>::CalcDShape;
  };





   class MyLinearQuad : public ScalarFiniteElement<2>
  {
    int vnums[3];
  public:
    // constructor
    MyLinearQuad (int order);

    virtual ELEMENT_TYPE ElementType() const { return ET_QUAD; }
    void SetVertexNumber (int i, int v) { vnums[i] = v; }
    
    /*
      Calculate the vector of shape functions in the point ip.
      ip is given in the reference element.
     */
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const;
  
    /*
      Calculate the matrix of derivatives of the shape functions in the point ip.
      dshape is an 4 by 2 matrix in our case.
     */
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             SliceMatrix<> dshape) const;

    // there are some more functions to bring in ...
    using ScalarFiniteElement<2>::CalcShape;    
    using ScalarFiniteElement<2>::CalcDShape;
  };



}

#endif

