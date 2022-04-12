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

      virtual void CalcMappedDDShape (const BaseMappedIntegrationPoint & mip, 
                                      BareSliceMatrix<> ddshape) const;

      // there are some more functions to bring in ...
      //using ScalarFiniteElement<2>::CalcShape;
      //using ScalarFiniteElement<2>::CalcDShape;
    };

    class LagrangePolyHornerCalc {
     protected:
        Array<double> nodes;

        Matrix<double> NewtonBasisCoeffs;
        Array<LagrangePolyHornerCalc> my_childs;
        void SetUpChilds();

     public:
        LagrangePolyHornerCalc() { }

        LagrangePolyHornerCalc(Array<double> & nodes2, bool deriv_calc_needed) : nodes(nodes2) { CalcNewtonBasisCoeffs(); if(deriv_calc_needed) SetUpChilds(); }

        void CalcNewtonBasisCoeffs();
        double Lagrange_Pol_Horner (double x, int i) const;
        double Lagrange_Pol_D_Horner (double x, int i) const;
    };

    class NodalTimeFE : public ScalarFiniteElement<1>
      {
      protected:
        int vnums[2];
        int k_t;
        bool skip_first_nodes = false;
        bool only_first_nodes = false;
        Array<double> nodes;

        bool do_horner_eval = false;
        LagrangePolyHornerCalc HornerLP;


      public:
        NodalTimeFE (int order, bool askip_first_nodes, bool aonly_first_nodes, int ndof_first_node = 1);
        virtual ELEMENT_TYPE ElementType() const { return ET_SEGM; }
        void SetVertexNumber (int i, int v) { vnums[i] = v; }

        virtual void CalcShape (const IntegrationPoint & ip,
                                BareSliceVector<> shape) const;

        virtual void CalcDShape (const IntegrationPoint & ip,
                                 BareSliceMatrix<> dshape) const;

        virtual bool IsNodeActive(int i) const
        {
          if (i<0 || i > k_t+1)
            throw Exception("node outside node range");
          if (i==0 && skip_first_nodes) 
            return false;
          else if (i!=0 && only_first_nodes)
            return false;
          else
            return true;
        }

        virtual void CalcInterpolationPoints ();
        virtual Array<double> & GetNodes() { return nodes; }
        virtual int order_time() const { return k_t; }

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

    class GCC3FE : public NodalTimeFE
      {
      public:
        GCC3FE (bool askip_first_nodes, bool aonly_first_nodes);

        virtual void CalcShape (const IntegrationPoint & ip,
                                BareSliceVector<> shape) const;

        virtual void CalcDShape (const IntegrationPoint & ip,
                                 BareSliceMatrix<> dshape) const;

        virtual bool IsNodeActive(int i) const
        {
          if (i<0 || i > 3+1 /*???*/)
            throw Exception("node outside node range");
          if (i<=1 && skip_first_nodes) 
            return false;
          else if (i>1 && only_first_nodes)
            return false;
          else
            return true;
        }

        virtual int order_time() const { return 3; }
        virtual void CalcInterpolationPoints ();

      };
  
 }


#endif

