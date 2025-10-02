#ifndef FILE_SPACETIMEFE_HPP
#define FILE_SPACETIMEFE_HPP

// SpaceTimeFE:

#include <fem.hpp>
#include "../utils/spacetimechecks.hpp"

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

      ScalarFiniteElement<1> * GetTimeFE() const { return tFE; }
      ScalarFiniteElement<D> * GetSpaceFE() const { return sFE; }

      virtual void CalcShape (const IntegrationPoint & ip,
                              BareSliceVector<> shape) const;

      virtual void CalcTimeShape (const IntegrationPoint & ip,
                                  BareSliceVector<> shape,
                                  int derivorder) const;    


      template <typename SpaceTimeShape, 
                typename SpaceCalcShape, typename SpaceShape, 
                typename TimeCalcShape, typename TimeShape>
      void GenericCalcShape(const IntegrationPoint& ip,
                            SpaceTimeShape shape,
                            SpaceCalcShape&& space_calcshape,
                            SpaceShape&& space_shape,
                            int MD,
                            TimeCalcShape&& time_calcshape,
                            TimeShape&& time_shape) const;


      template <typename SpaceTimeShape, 
                typename SpaceCalcShape, typename SpaceShape, 
                typename TimeCalcShape, typename TimeShape>
      void GenericCalcMappedShape(const BaseMappedIntegrationPoint& mip,
                            SpaceTimeShape shape,
                            SpaceCalcShape&& space_calcshape,
                            SpaceShape&& space_shape,
                            int MD,
                            TimeCalcShape&& time_calcshape,
                            TimeShape&& time_shape) const;

      // for time derivatives

      virtual void CalcDtShape (const IntegrationPoint & ip,
                               BareSliceVector<> dshape) const;

      virtual void CalcDDtShape (const IntegrationPoint & ip,
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
  


  template <int D>
  template <typename SpaceTimeShape, 
            typename SpaceCalcShape, typename SpaceShape, 
            typename TimeCalcShape, typename TimeShape>
  void SpaceTimeFE<D>::GenericCalcShape(
      const IntegrationPoint& ip,
      SpaceTimeShape shape,
      SpaceCalcShape&& space_calcshape,
      SpaceShape&& space_shape,
      int MD,
      TimeCalcShape&& time_calcshape,
      TimeShape&& time_shape
  ) const
  {
      if(!IsSpaceTimeIntegrationPoint(ip))
          throw Exception("SpaceTimeFE :: CalcShape called with a mere space IR");
      if (tFE->Order() == 0) {
          space_calcshape(ip, shape);
      } else {
          IntegrationPoint z(override_time ? time : ip.Weight());


          //Vector<> time_shape(tFE->GetNDof());
          time_calcshape(z, time_shape);
          FlatMatrix<> t_shape(tFE->GetNDof(), 1, &time_shape(0));
          //Vector<> space_shape(sFE->GetNDof());
          space_calcshape(ip, space_shape);
          FlatMatrix<> x_shape(sFE->GetNDof(), MD, &space_shape(0));

          int ii = 0;
          for (int j = 0; j < tFE->GetNDof(); j++) {
              for (int i = 0; i < sFE->GetNDof(); i++) {
                  for (int dimi = 0; dimi < MD; dimi++)
                    shape(ii,dimi) = x_shape(i,dimi)*t_shape(j,0);
                  ii++; 
              }
          }
      }
  }     

  template <int D>
  template <typename SpaceTimeShape, 
            typename SpaceCalcShape, typename SpaceShape, 
            typename TimeCalcShape, typename TimeShape>
  void SpaceTimeFE<D>::GenericCalcMappedShape(
      const BaseMappedIntegrationPoint& mip,
      SpaceTimeShape shape,
      SpaceCalcShape&& space_calcshape,
      SpaceShape&& space_shape,
      int MD,
      TimeCalcShape&& time_calcshape,
      TimeShape&& time_shape
  ) const
  {
      if (tFE->Order() == 0) {
          space_calcshape(mip, shape);
      } else {
          IntegrationPoint z(override_time ? time : mip.IP().Weight());

          //if(!IsSpaceTimeIntegrationPoint(ip)) /// <- re-introduce this
          //    throw Exception("SpaceTimeFE :: CalcShape called with a mere space IR");

          //Vector<> time_shape(tFE->GetNDof());
          time_calcshape(z, time_shape);
          FlatMatrix<> t_shape(tFE->GetNDof(), 1, &time_shape(0));
          //Vector<> space_shape(sFE->GetNDof());
          space_calcshape(mip, space_shape);
          FlatMatrix<> x_shape(sFE->GetNDof(), MD, &space_shape(0));

          int ii = 0;
          for (int j = 0; j < tFE->GetNDof(); j++) {
              for (int i = 0; i < sFE->GetNDof(); i++) {
                  for (int dimi = 0; dimi < MD; dimi++)
                    shape(ii,dimi) = x_shape(i,dimi)*t_shape(j,0);
                  ii++; 
              }
          }
      }
  }     


 }

#endif

