#ifndef FILE_SPACETIMEFE_HPP
#define FILE_SPACETIMEFE_HPP

// SpaceTimeFE based on:

/*********************************************************************/
/* File:   myElement.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/


#include <fem.hpp>

namespace ngfem
{

  template <int D>
    class SpaceTimeFE : public ScalarFiniteElement<D>
   {
        ScalarFiniteElement<D>* sFE = nullptr;
        ScalarFiniteElement<1>* tFE = nullptr;
        double time;
        bool override_time = false;

    public:
      // constructor
      SpaceTimeFE (ScalarFiniteElement<D>* s_FE,ScalarFiniteElement<1>*t_FE, bool override_time, double time );

      virtual ELEMENT_TYPE ElementType() const { return sFE->ElementType(); }


      virtual void CalcShape (const IntegrationPoint & ip,
                              BareSliceVector<> shape) const;

      // for time derivatives

      virtual void CalcDtShape (const IntegrationPoint & ip,
                               BareSliceVector<> dshape) const;

      virtual void CalcDShape (const IntegrationPoint & ip,
                               BareSliceMatrix<> dshape) const;

      // there are some more functions to bring in ...
      //using ScalarFiniteElement<2>::CalcShape;
      //using ScalarFiniteElement<2>::CalcDShape;
    };


    class NodalTimeFE : public ScalarFiniteElement<1>
      {
        int vnums[2];
        int k_t;
        bool skip_first_node = false;
        bool only_first_node = false;
        Array<double> nodes;

      public:
        NodalTimeFE (int order, bool askip_first_node, bool aonly_first_node);
        virtual ELEMENT_TYPE ElementType() const { return ET_SEGM; }
        void SetVertexNumber (int i, int v) { vnums[i] = v; }

        virtual void CalcShape (const IntegrationPoint & ip,
                                BareSliceVector<> shape) const;

        virtual void CalcDShape (const IntegrationPoint & ip,
                                 BareSliceMatrix<> dshape) const;

        bool IsNodeActive(int i) const
        {
          if (i<0 || i > k_t+1)
            throw Exception("node outside node range");
          if (i==0 && skip_first_node) 
            return false;
          else if (i!=0 && only_first_node)
            return false;
          else
            return true;
        }

        void CalcInterpolationPoints ();
        Array<double> & GetNodes() { return nodes; }
        int order_time() const { return k_t; }

        template <class T>
        T Lagrange_Pol (T x, int i) const
        {
           T result  = 1;
           for (int j = 0; j < nodes.Size(); j++) {
               if ( j != i)
                    result *= ( x-nodes[j] ) / ( nodes[i] - nodes[j] );
           }

           return result;
        }

      };

 }


#endif

