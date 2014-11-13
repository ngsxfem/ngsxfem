#ifndef FILE_SPACETIME_PRECOND
#define FILE_SPACETIME_PRECOND


/*********************************************************************/
/* File:   stprecond.hpp                                             */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   June 2014                                                 */
/*********************************************************************/


namespace ngcomp
{

  
class SpaceTimePreconditioner : public Preconditioner
{
  shared_ptr<S_BilinearForm<double>> bfa;

  // Array<int> global_nums;
  // int ilower, iupper;
  SparseMatrix<double> * AssBlock = NULL;
  shared_ptr<BaseMatrix> InvAssBlock = NULL;
  string inversetype = "pardiso";
  Table<int> * blockjacobixtable = NULL;
  BaseBlockJacobiPrecond * blockjacobix = NULL;
  //Vector xdiag;

  const BitArray * freedofs;
  // const ParallelDofs * pardofs;

public:
    
  SpaceTimePreconditioner (const PDE & pde, const Flags & flags, const string & name);
  SpaceTimePreconditioner (const BaseMatrix & matrix, const BitArray * afreedofs); 

  ~SpaceTimePreconditioner ();
	
  virtual void Update();
  virtual void Mult (const BaseVector & f, BaseVector & u) const;

  virtual const BaseMatrix & GetAMatrix() const
  {
    return bfa->GetMatrix(); 
  }

  virtual const char * ClassName() const
  { return "SPACETIME Preconditioner"; }

private:
  void Setup (const BaseMatrix & matrix);
};

}



#endif






