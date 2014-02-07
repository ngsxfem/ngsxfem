#ifndef FILE_SPACETIMEFE_HPP
#define FILE_SPACETIMEFE_HPP

/// from ngsolve
#include <comp.hpp>
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

using namespace ngsolve;

namespace ngfem
{

typedef std::pair<double,double> TimeInterval;

/**
     
 */
  class SpaceTimeFiniteElement : public FiniteElement
  {
  protected:
    const FiniteElement & base_space;
    const FiniteElement & base_time;
    int order_space;
    int order_time;
    int ndof_space;
    int ndof_time;
  public:
    SpaceTimeFiniteElement(const FiniteElement & base_space,
                           const FiniteElement & base_time);
    virtual ~SpaceTimeFiniteElement();

    /// the name
    virtual string ClassName(void) const;

    int GetNDofSpace() const { return ndof_space;}
    int GetNDofTime() const { return ndof_time;}
    int OrderSpace() const { return order_space;}
    int OrderTime() const { return order_time;}
    virtual int GetNDof() const { return ndof;}

  };

  template <int D>  
  class ScalarSpaceTimeFiniteElement : public SpaceTimeFiniteElement
  {
  protected:
    const ScalarFiniteElement<D> & scalar_space;
    const DGFiniteElement<1> & scalar_time;
  public:
    ScalarSpaceTimeFiniteElement(const ScalarFiniteElement<D> & base_space,
                                 const DGFiniteElement<1> & base_time);
    virtual ~ScalarSpaceTimeFiniteElement();

    virtual int GetNDof() const { return ndof;}

    /// the name
    virtual string ClassName(void) const;

    virtual bool IsDGFiniteElement() const;
    /// mass diag
    virtual void GetDiagMassMatrix(FlatVector<> diagmass,
                                   LocalHeap & lh) const;

    /// compute shape
    virtual void CalcShapeTime (double time,
                                FlatVector<> shape) const;
    /// compute shape
    virtual void CalcShapeSpace (const IntegrationPoint & ip,
                                 FlatVector<> shape) const;
    /// compute shape
    virtual void CalcShapeSpaceTime (const IntegrationPoint & ip, double time,
                                     FlatVector<> shape, LocalHeap & lh) const;

    /// compute dshape, vector: ndof 
    virtual void CalcDtShapeSpaceTime (const IntegrationPoint & ip, double time,
                                       FlatVector<> dshape, LocalHeap & lh) const;

    virtual ELEMENT_TYPE ElementType() const{ return base_space.ElementType(); }

    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcDxShapeSpaceTime (const IntegrationPoint & ip, double time,
                                       FlatMatrixFixWidth<D> dshape, LocalHeap & lh) const;

    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcMappedDxShapeSpaceTime (const SpecificIntegrationPoint<D,D> & sip, double time,
                                             FlatMatrixFixWidth<D> dshape, LocalHeap & lh) const;
  };
    
    // /// compute dshape, vector: ndof 
    // virtual void CalcMappedDtShapeSpaceTime (const SpecificIntegrationPoint<D,D> & sip, double time
    //                                          FlatVector<> dshape) const;

    // /// compute dshape, matrix: ndof x (spacedim spacedim)
    // // virtual void CalcDDShape (const IntegrationPoint & ip, 
    // //                           FlatMatrix<> ddshape) const;

  template <int D>
  class ScalarSpaceTimeFEEvaluator
  {
  protected:
    const ScalarSpaceTimeFiniteElement<D> & fe;
    FlatVector<> linvec;
    IntegrationPoint ip;
    FlatVector<> shape;
    LocalHeap & lh;
  public:
    ScalarSpaceTimeFEEvaluator(const ScalarSpaceTimeFiniteElement<D> & a_fe, FlatVector<> a_linvec, LocalHeap & a_lh)
      : fe(a_fe), linvec(a_linvec), lh(a_lh), shape(linvec.Size(), lh)
    {
    }
      
    double operator()(const Vec<D+1>& point);
  };



  /// Hat hier nichts zu suchen

  template <int D>
  class Simplex
  {
  public:
    Array< Vec<D> > p;
    Simplex(const Array< Vec<D> > & a_p): p(a_p)
    {
      ;
    }
  };

  template <int D>
  inline ostream & operator<< (ostream & ost, const Simplex<D> & s)
  {
    for (int i = 0; i < D+1; i++)
      ost << i << ":" << s.p[i] << "\t";
    ost << endl;
    return ost;
  }


  template<int D>
  void DecomposeIntoSimplices(ELEMENT_TYPE et, Array<Simplex<D+1> *>& ret, LocalHeap & lh)
  {
    ret.SetSize(0);
	switch (et)
	  {
	  case ET_TRIG:
	  case ET_TET:
      {
        const POINT3D * verts = ElementTopology::GetVertices(et);
        Array< Vec<D+1> > p(2*(D+1));
        for (int i = 0; i < D+1; ++i)
        {
          for (int d = 0; d < D; ++d)
          {
            p[i][d] = verts[i][d];
            p[i+D+1][d] = verts[i][d];
          }
          p[i][D] = 0.0;
          p[i+D+1][D] = 1.0;
        }

        Array< Vec<D+1> > tet(D+2);
        for (int i = 0; i < D+1; ++i)
        {
          for (int j = 0; j < D+2; ++j)
            tet[j] = p[i+j];
          ret.Append(new (lh) Simplex<D+1> (tet));
        }

        cout << " report \n";
        for (int i = 0; i < ret.Size(); ++i)
          cout << " simplex " << i << ": " << *ret[i] << endl;
      }
      break;
      default:  // for the compiler
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
        break;
	  }
  }

  template<int D>
  void Decompose(ELEMENT_TYPE et, const ScalarSpaceTimeFEEvaluator<D> & stfeeval)
  {

    
	switch (et)
	  {
	  case ET_TRIG:
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
	    break;
	  case ET_QUAD:
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
	    break;
	  case ET_TET:
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
	    break;
      default:  // for the compiler
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
        break;
	  }


  }




} // end of namespace

#endif
