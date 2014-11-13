template<bool inv>
inline void ApplyA(BaseMatrix & A, const BaseVector & v, BaseVector & w, 
              const BitArray * freedofs, FlatVector<double> fvdiaga)
{
  FlatVector<double> fvw = w.FVDouble();
  fvw = 0.0;
  A.MultAdd(1.0,v,w);

  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fvw(i) = 0.0;
    else
      if (inv)
        fvw(i) = fvw(i) * fvdiaga(i);
      else
        fvw(i) = fvw(i) / fvdiaga(i);
}

inline void ApplyATA(BaseMatrix & A, BaseMatrix * AT, const BaseVector & v, BaseVector & w, 
              const BitArray * freedofs, FlatVector<double> fvdiaga)
{
  BaseVector & c = *A.CreateVector();
  FlatVector<double> fvc = c.FVDouble();
  FlatVector<double> fvw = w.FVDouble();
  c = 0.0;
  w = v;
  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fvw(i) = 0.0;
    else
      fvw(i) = fvw(i) /sqrt(fvdiaga(i));

  A.MultAdd(1.0,w,c);

  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fvc(i) = 0.0;
    else
      fvc(i) = fvc(i) /(fvdiaga(i));

  w = 0.0;
  if(AT!=NULL)
    AT->MultAdd(1.0,c,w);
  else
    A.MultTransAdd(1.0,c,w);

  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fvw(i) = 0.0;
    else
      fvw(i) = fvw(i) /sqrt(fvdiaga(i));

  // for (int i = 0; i < freedofs->Size(); ++i)
  //   if (!freedofs->Test(i))
  //     fvw(i) = 0.0;

  // delete &d;
  delete &c;
}

inline void ApplyAAT(BaseMatrix & A, BaseMatrix * AT, const BaseVector & v, BaseVector & w, 
              const BitArray * freedofs, FlatVector<double> fvdiaga)
{
  BaseVector & c = *A.CreateVector();
  FlatVector<double> fvc = c.FVDouble();
  FlatVector<double> fvw = w.FVDouble();
  c = 0.0;
  w = v;

  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fvw(i) = 0.0;
    else
      fvw(i) = fvw(i) * sqrt(fvdiaga(i));

  if(AT!=NULL)
    AT->MultAdd(1.0,w,c);
  else
    A.MultTransAdd(1.0,w,c);
  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fvc(i) = 0.0;
    else
      fvc(i) = fvc(i) * fvdiaga(i);

  w = 0.0;
  A.MultAdd(1.0,c,w);

  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fvw(i) = 0.0;
    else
      fvw(i) = fvw(i) * sqrt(fvdiaga(i));
  delete &c;
}

template <bool symm>
inline double PowerIteration(BaseMatrix & A, const BitArray * freedofs, FlatVector<double> fvdiaga, Array<BaseVector*> & evs, bool printmuch = false)
{
  const int max_its = 100000;
  const double rel_acc = 1e-6;

  BaseVector & a = *A.CreateVector();
  BaseVector & b = *A.CreateVector();
  BaseVector & c = *A.CreateVector();

  FlatVector<double> fva = a.FVDouble();
  // FlatVector<double> fvb = b.FVDouble();
 

  a.SetRandom();
  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fva(i) = 0.0;
  double an = L2Norm(a);
  a /= an;

  double bn; //lambda
  double bn_last = an; //lambda
  if (printmuch)
    cout << endl;
  for (int i = 0; i < max_its; ++i)
  {
    bn_last = bn;

    for (int j = 0; j < evs.Size(); ++j)
    {
      const double ab = InnerProduct(a,*(evs[j]) );
      a -= ab * *(evs[j]);
    }

    if (symm)
      ApplyA<false>(A,a,b,freedofs,fvdiaga);
    else
    {
      // ApplyA<false>(A,a,c,freedofs,fvdiaga);
      
      // ApplyA<false>(A,c,b,freedofs,fvdiaga);
      ApplyATA(A,NULL,a,b,freedofs, fvdiaga);
    }

    bn = L2Norm(b);
    b /= bn;
    a = b;
    if (printmuch)
      std::cout << ""
        // << "\r                         \r"
                << "it = " << i << ", bn = " << bn << endl;

      // std::cout << "\r                         \r"
      //           << "it = " << i << ", bn = " << bn;
    if (abs(bn-bn_last) < rel_acc * abs(bn_last)) break;
  }
  evs.Append(&a);
  delete &b;
  delete &c;
  // getchar();
  if (symm)
    return bn;
  else
    return sqrt(bn);  
}


template <bool symm>
inline double InversePowerIteration(BaseMatrix & A, BaseMatrix & invA, const BitArray * freedofs, FlatVector<double> fvdiaga, 
                             Array<BaseVector*> & invevs, bool printmuch = false)
{
  throw Exception("transpose inverse hack deactivated right now");
  BaseMatrix & invAT = * dynamic_cast<BaseSparseMatrix&> (A) . InverseMatrix(freedofs);
  // for transpose:
  // BaseMatrix & invAT = * dynamic_cast<BaseSparseMatrix&> (A) . InverseMatrix(freedofs,true);

  const int max_its = 100000;
  const double rel_acc = 1e-6;
  BaseVector & a = *A.CreateVector();
  BaseVector & b = *A.CreateVector();

  FlatVector<double> fva = a.FVDouble();
  // FlatVector<double> fvb = b.FVDouble();

  double an;
  double bn; //lambda

  a.SetRandom();
  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fva(i) = 0.0;
  an = L2Norm(a);
  a /= an;

  if (printmuch)
    cout << endl;

  double bn_last = an; //lambda

  for (int i = 0; i < max_its; ++i)
  {
    bn_last = bn;

    for (int j = 0; j < invevs.Size(); ++j)
    {
      const double ab = InnerProduct(a,*(invevs[j]) );
      a -= ab * *(invevs[j]);
    }

    if (symm)
      ApplyA<true>(invA,a,b,freedofs,fvdiaga);
    else
      ApplyAAT(invA,&invAT,a,b,freedofs, fvdiaga);

    bn = L2Norm(b);
    b /= bn;
    a = b;
    if (printmuch)
      std::cout << ""
        // << "\r                         \r"
                << "it = " << i << ", bn = " << bn << endl;

      // std::cout << "\r                         \r"
      //           << "it = " << i << ", bn = " << bn;
    if (abs(bn-bn_last) < rel_acc * abs(bn_last)) break;
  }

  invevs.Append(&b);

  delete &a;
  delete &invAT;
  // getchar();
  if (symm)
    return 1.0/bn;
  else
    return 1.0/sqrt(bn);
}






inline void CalcCond(BaseMatrix & A, BaseMatrix & invA, const BitArray * freedofs, bool printmuch = true, bool jacobiprec = false, ofstream * outs = NULL, bool symmetric = false, shared_ptr<GridFunction> gfu = NULL)
{
  Array<BaseVector*> evs(0);
  Array<BaseVector*> invevs(0);

  BaseVector & diaga = *A.CreateVector();
  FlatVector<double> fvdiaga = diaga.FVDouble();
  fvdiaga = 1.0;
  if (jacobiprec)
  {
    SparseMatrixTM<double> * AA = dynamic_cast<SparseMatrixTM<double> *>(&A);
    if (AA != NULL)
    {
      for (int i = 0; i < freedofs->Size(); ++i)
      {
        if (freedofs->Test(i))
          fvdiaga(i) = (*AA)(i,i);
      }
    }
    else
      throw Exception("no diag access....");
  }

  const double cup = symmetric ? 
    PowerIteration<true>(A, freedofs, fvdiaga, evs, false)
    : PowerIteration<false>(A, freedofs, fvdiaga, evs, false);

  // std::cout << "\n lambda max_0 : " << cup << std::endl;
  // for (int i = 1; i < 100; ++i)
  //   std::cout << "\n lambda max_" << i << " : " << PowerIteration( A, freedofs, fvdiaga, evs, true) << std::endl;

  const double clow = symmetric ? 
    InversePowerIteration<true>( A, invA, freedofs, fvdiaga, invevs, false)
    : InversePowerIteration<false>( A, invA, freedofs, fvdiaga, invevs, false);

  // std::cout << "\n lambda min_0 : " << cup << std::endl;
  // for (int i = 1; i < 100; ++i)
  //   std::cout << "\n lambda min_" << i << " : " << InversePowerIteration( A, invA, freedofs, fvdiaga, invevs, true) << std::endl;

  // const double clow = 1.0/sqrt(bn);

  cout << "condition number: " << std::setw(12) << cup/clow 
       << " (up: " << std::setw(12) << cup << ", low: " << std::setw(12) << clow << ")";
  if (printmuch)
    cout << endl;

  if (outs)
  {
    *outs << cup/clow << "\t";
  }

  if (gfu)
    gfu->GetVector() = * invevs[0];

  // if (gfu)
  //   gfu->GetVector() = * evs[0];
  
  for (int i = 0; i < evs.Size(); ++i)
    delete evs[i];

  for (int i = 0; i < invevs.Size(); ++i)
    delete invevs[i];

  delete &diaga;
}

