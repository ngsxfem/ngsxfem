/*********************************************************************/
/* File:   myElement.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*
  
My own simple first and second order triangular finite elements

*/


#include <fem.hpp>
#include "myElement.hpp"
#include "myHOElement.hpp"


namespace ngfem
{
  

  MyLinearSegm :: MyLinearSegm ()
  /*
    Call constructor for base class: 
    geometry is ET_SEGM, number of dofs is 2, maximal order is 1
   */
    : ScalarFiniteElement<1> (2, 1)
  { ; }


  void MyLinearSegm :: CalcShape (const IntegrationPoint & ip, 
                                  BareSliceVector<> shape) const
  {
    // coordinates in reference elements
    double x = ip(0);

    /*
      Vertex coordinates have been defined to be (1), (0),
      see ElementTopology::GetVertices(ET_SEGM)
     */

    // define shape functions
    shape(0) = x;
    shape(1) = 1-x;
  }
 
  void MyLinearSegm :: CalcDShape (const IntegrationPoint & ip, 
                                   SliceMatrix<> dshape) const   
  {
    // matrix of derivatives:

    dshape(0,0) = 1;
    dshape(1,0) = -1;
  }

 

  

  MyLinearQuad :: MyLinearQuad (int order)
    /*
      Call constructor for base class: 
      geometry is ET_QUAD, number of dofs is 4, maximal order is 1
    */
    : ScalarFiniteElement<2> ((order+1)*(order+1), order)
  { ; }


  void MyLinearQuad :: CalcShape (const IntegrationPoint & ip, 
                                  BareSliceVector<> shape) const
  {
      // Import shape functions from MyLinearSegm
      MyHighOrderSegm space_x(order),space_y(order);
      IntegrationPoint ipx(ip(0));
      IntegrationPoint ipy(ip(1));
      Vector<double> shape_x(order+1),shape_y(order+1);
      space_x.CalcShape(ipx,shape_x);
      space_y.CalcShape(ipy,shape_y);
      
   /* shape functions
      shape(0) = (1-x)*(1-y);
      shape(1) = x*(1-y);
      shape(2) = x*y;
      shape(3) = (1-x)*y;  */
 
    const POINT3D * vertices = ElementTopology::GetVertices (ET_QUAD);
    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD); 

    // shape functions are products of the shape functions from MyLinearSegm

    for(int i = 0; i<4; i++) {
       int xi = vertices[i][0];
       int yi = vertices[i][1];
       shape(i) = shape_x(1-xi) * shape_y(1-yi);
    } 
       if ( order == 2) {
 
            int ii = 4;    
   	   for (int i = 0; i < 4; i++) { 
	       int es = edges[i][0], ee = edges[i][1];

               
	       if (vnums[es] > vnums[ee]) swap (es, ee);                
               // temporary solution
               if ( (es == 0 && ee == 1) || (es == 1 && ee == 0) )
                       shape[ii++] = shape_x(2)*shape_y(1);
               if ( (es == 2 && ee == 3) || (es == 3 && ee == 2) )
                       shape[ii++] = shape_x(2)*shape_y(0);
               if ( (es == 3 && ee == 0) || (es == 0 && ee == 3) )
                      shape[ii++] = shape_x(1)*shape_y(2);
                if ( (es == 1 && ee == 2) || (es == 2 && ee == 1) ) 
                      shape[ii++] = shape_x(0)*shape_y(2);
                                            	       
	    }   

            shape(8) = shape_x(2)*shape_y(2);
      }
       
  }
  
  void MyLinearQuad :: CalcDShape (const IntegrationPoint & ip, 
                                   SliceMatrix<> dshape) const
    
  {
    // import shape functions from MyLinearSegm
    MyHighOrderSegm space_x(2),space_y(2);
    IntegrationPoint ipx(ip(0));
    IntegrationPoint ipy(ip(1));
    Vector<double> shape_x(order+1),shape_y(order+1);
    Matrix<double> dshape_x(order+1,2),dshape_y(order+1,2);
    space_x.CalcShape(ipx,shape_x);
    space_x.CalcDShape(ipx,dshape_x);
    space_y.CalcShape(ipy,shape_y);
    space_y.CalcDShape(ipy,dshape_y);

    // matrix of derivatives:
    /* dshape(0,0) = -(1-y);
    dshape(0,1) = -(1-x);
    dshape(1,0) = 1-y;
    dshape(1,1) = -x;
    dshape(2,0) = y;
    dshape(2,1) = x;
    dshape(3,0) = -y;
    dshape(3,1) = 1-x; */
   
     const POINT3D * vertices = ElementTopology::GetVertices (ET_QUAD);
     const EDGE * edges = ElementTopology::GetEdges (ET_QUAD); 
    
    // Using the tensor product structure the derivatives act
    // only on one of the functions
    for(int i = 0; i<4; i++) {
       int xi = vertices[i][0];
       int yi = vertices[i][1];
       dshape(i,0) = dshape_x(1-xi,0) * shape_y(1-yi);
       dshape(i,1) = shape_x(1-xi) * dshape_y(1-yi,0);
    } 

   if ( order == 2) {
 
            int ii = 4;    
   	   for (int i = 0; i < 4; i++) { 
	       int es = edges[i][0], ee = edges[i][1];

               
	       if (vnums[es] > vnums[ee]) swap (es, ee);                
               // temporary solution
               if ( (es == 0 && ee == 1) || (es == 1 && ee == 0) ) {
                       dshape(ii,0) = dshape_x(2,0)*shape_y(1);
                       dshape(ii,1) = shape_x(2)*dshape_y(1,0);
                       ii++;
               }
               if ( (es == 2 && ee == 3) || (es == 3 && ee == 2) ) {
                       dshape(ii,0) = dshape_x(2,0)*shape_y(0);
                       dshape(ii,1) = shape_x(2)*dshape_y(0,0);
                       ii++;
               }
               if ( (es == 3 && ee == 0) || (es == 0 && ee == 3) ) {
                      dshape(ii,0) = dshape_x(1,0)*shape_y(2);
                      dshape(ii,1) = shape_x(1)*dshape_y(2,0);
                      ii++;
               }
                if ( (es == 1 && ee == 2) || (es == 2 && ee == 1) ) {
                      dshape(ii,0) = dshape_x(0,0)*shape_y(2);
                      dshape(ii,1) = shape_x(0)*dshape_y(2,0);
                      ii++;
                }
                                            	       
	    }   

            dshape(8,0) = dshape_x(2,0)*shape_y(2);
            dshape(8,1) = shape_x(2)*dshape_y(2,0);
            
      }

   
      
    
  }


}
