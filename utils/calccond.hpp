void ApplyATA(BaseMatrix & A, BaseMatrix * AT, const BaseVector & v, BaseVector & w, const BitArray * freedofs)
{
  BaseVector & c = *A.CreateVector();
  FlatVector<double> fvc = c.FVDouble();
  FlatVector<double> fvw = w.FVDouble();
  c = 0.0;
  A.MultAdd(1.0,v,c);

  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fvc(i) = 0.0;

  w = 0.0;
  if(AT!=NULL)
    AT->MultAdd(1.0,c,w);
  else
    A.MultTransAdd(1.0,c,w);

  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fvw(i) = 0.0;

  delete &c;
}

void ApplyAAT(BaseMatrix & A, BaseMatrix * AT, const BaseVector & v, BaseVector & w, const BitArray * freedofs)
{
  BaseVector & c = *A.CreateVector();
  // FlatVector<double> fvc = c.FVDouble();
  // FlatVector<double> fvw = w.FVDouble();
  c = 0.0;
  if(AT!=NULL)
    AT->MultAdd(1.0,v,c);
  else
    A.MultTransAdd(1.0,v,c);
  // for (int i = 0; i < freedofs->Size(); ++i)
  //   if (!freedofs->Test(i))
  //     fvc(i) = 0.0;
  w = 0.0;
  A.MultAdd(1.0,c,w);
  // for (int i = 0; i < freedofs->Size(); ++i)
  //   if (!freedofs->Test(i))
  //     fvw(i) = 0.0;
  delete &c;
}

void CalcCond(BaseMatrix & A, BaseMatrix & invA, const BitArray * freedofs)
{
  // BaseMatrix AT(A);
  // AT.AsVector() = 0.0;
  
  BaseMatrix & invAT = * dynamic_cast<BaseSparseMatrix&> (A) . InverseMatrix(freedofs,true);

  const int max_its = 100000;
  const double rel_acc = 1e-6;
  BaseVector & a = *A.CreateVector();
  BaseVector & b = *A.CreateVector();
  BaseVector & c = *A.CreateVector();

  FlatVector<double> fva = a.FVDouble();
  FlatVector<double> fvb = b.FVDouble();
  FlatVector<double> fvc = c.FVDouble();
  // a.SetRandom();
  // b = invAT * a;
  // c = invAT * a;
  // a = b - c;
  // std::cout << " L2Norm(a) = " << L2Norm(a) << std::endl;
 

  a.SetRandom();
  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fva(i) = 0.0;
  double an = L2Norm(a);
  a /= an;

  double bn; //lambda
  double bn_last = an; //lambda
  for (int i = 0; i < max_its; ++i)
  {
    bn_last = bn;
    ApplyATA(A,NULL,a,b,freedofs);
    bn = L2Norm(b);
    b /= bn;
    a = b;
    std::cout << "\r                         \r"
              << "it = " << i << ", bn = " << bn;
    if (abs(bn-bn_last) < rel_acc * abs(bn_last)) break;
  }
  // cout << endl;
  // cout << " sqrt(lambda_max(ATA)) = " << sqrt(bn) << endl;
  const double cup = sqrt(bn);

  a.SetRandom();
  for (int i = 0; i < freedofs->Size(); ++i)
    if (!freedofs->Test(i))
      fva(i) = 0.0;
  an = L2Norm(a);
  a /= an;

  
  bn_last = an; //lambda
  for (int i = 0; i < max_its; ++i)
  {
    bn_last = bn;
    ApplyATA(invA,&invAT,a,b,freedofs);
    // b = invA * a;
    // double sum = 0.0;
    // int cnt = 0;
    // for (int j = 0; j < freedofs->Size(); ++j)
    //   if (!freedofs->Test(j))
    //     fvb(j) = 0.0;
    //   else
    //   {
    //     sum += fvb(j);
    //     cnt ++;
    //   }
    // for (int j = 0; j < freedofs->Size(); ++j)
    //   if (freedofs->Test(j))
    //     fvb(j) -= sum/cnt;

    bn = L2Norm(b);
    b /= bn;
    a = b;
    // std::cout << " a = " << a << std::endl;
    // getchar();
    std::cout << "\r                         \r"
              << "it = " << i << ", bn = " << bn;
    if (abs(bn-bn_last) < rel_acc * abs(bn_last)) break;
  }
  // cout << endl;
  // cout << " sqrt(lambda_min(ATA)) = " << 1.0/sqrt(bn) << endl;
  const double clow = 1.0/sqrt(bn);

  // ApplyATA(A,a,b);
  // ApplyAAT(A,b,c);
  // while (

  cout << "\r condition number is = " << cup/clow 
       << " (upper bound: " << cup << ", lower bound: " << clow << ")" << endl;

  delete &a;
  delete &b;
  delete &c;
  delete &invAT;
}
