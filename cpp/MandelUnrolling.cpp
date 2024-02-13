//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
 //                 ``This    source  code is Copyrighted in
 //                 U.S.,  for  an  indefinite  period,  and anybody
 //                 caught  using it without our permission, will be
 //                 mighty good friends of ourn, cause we don't give
 //                 a  darn.  Hack it. Compile it. Debug it. Run it.
 //                 Yodel  it.  Enjoy it. We wrote it, that's all we
 //                 wanted to do.''
 //
 //
 // COPYRIGHT (C):     :-))
 // PROJECT:           Object Oriented Finite Element Program
 // FILE:
 // CLASS:
 // MEMBER FUNCTIONS:
 //
 // MEMBER VARIABLES
 //
 // PURPOSE:
 //
 // RETURN:
 // VERSION:
 // LANGUAGE:          C++
 // TARGET OS:
 // DESIGNER:          Zhao Cheng, Boris Jeremic
 // PROGRAMMER:        Zhao Cheng
 // DATE:              Fall 2005
 // UPDATE HISTORY:    06/2006, add functions for matrix based elements, CZ
 //                    10/2006, add various more algorithms, CZ
 //                    Guanzhou Jie updated for parallel Dec 2006
 //                    Mahdi Taiebat & Boris Jeremic debugged the code for the
 //                    implicit method, April2007
 //                    Feb2008 adding ScaledExplicit, Mahdi's idea
 //                    Nima Tafazzoli updated for API (Feb 2009)

int NewTemplate3Dep::Tensor2MatrixSysR4(const tensor& T, Matrix& M)
 {
   int rank = T.rank();
   if (rank != 4) {
     cout << "NewTemplate3Dep::Tensor2MatrixSysR4 - tensor must be of rank 4" << endln;
     return 1;
   }

   int nr = M.noRows();
   int nc = M.noCols();
   if (nr != 6 || nc != 6) {
     cout << "NewTemplate3Dep::Tensor2MatrixSysR4 - matrix must be of (6, 6)" << endln;;
     return 1;
   }

   double sqrt2 = sqrt(2.0);
   double two = 2.0;

   // Adopt method from Helnwein (2001):

   M(0,0) = T.cval(1,1,1,1);
   M(0,1) = T.cval(1,1,2,2);
   M(0,2) = T.cval(1,1,3,3);
   M(0,3) = T.cval(1,1,1,2) *sqrt2;
   M(0,4) = T.cval(1,1,2,3) *sqrt2;
   M(0,5) = T.cval(1,1,1,3) *sqrt2;

   M(1,0) = T.cval(2,2,1,1);
   M(1,1) = T.cval(2,2,2,2);
   M(1,2) = T.cval(2,2,3,3);
   M(1,3) = T.cval(2,2,1,2) *sqrt2;
   M(1,4) = T.cval(2,2,2,3) *sqrt2;
   M(1,5) = T.cval(2,2,1,3) *sqrt2;

   M(2,0) = T.cval(3,3,1,1);
   M(2,1) = T.cval(3,3,2,2);
   M(2,2) = T.cval(3,3,3,3);
   M(2,3) = T.cval(3,3,1,2) *sqrt2;
   M(2,4) = T.cval(3,3,2,3) *sqrt2;
   M(2,5) = T.cval(3,3,1,3) *sqrt2;

   M(3,0) = T.cval(1,2,1,1) *sqrt2;
   M(3,1) = T.cval(1,2,2,2) *sqrt2;
   M(3,2) = T.cval(1,2,3,3) *sqrt2;
   M(3,3) = T.cval(1,2,1,2) *two;
   M(3,4) = T.cval(1,2,2,3) *two;
   M(3,5) = T.cval(1,2,1,3) *two;

   M(4,0) = T.cval(2,3,1,1) *sqrt2;
   M(4,1) = T.cval(2,3,2,2) *sqrt2;
   M(4,2) = T.cval(2,3,3,3) *sqrt2;
   M(4,3) = T.cval(2,3,1,2) *two;
   M(4,4) = T.cval(2,3,2,3) *two;
   M(4,5) = T.cval(2,3,1,3) *two;

   M(5,0) = T.cval(1,3,1,1) *sqrt2;
   M(5,1) = T.cval(1,3,2,2) *sqrt2;
   M(5,2) = T.cval(1,3,3,3) *sqrt2;
   M(5,3) = T.cval(1,3,1,2) *two;
   M(5,4) = T.cval(1,3,2,3) *two;
   M(5,5) = T.cval(1,3,1,3) *two;

     return 0;
 }
