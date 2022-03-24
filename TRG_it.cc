// Last modified date : 2020-02-12

// How to improve code? (later?)
// 1. Make functions so that we downsize the length of hotrg function. 
// 2. Save some logs useful (Truncation information i.e. error, dimension change, S.value spectrum,...)


#include "itensor/all.h"
#include <iostream>
#include <cstring>
#include <vector>
#include <list>
#include <sys/time.h>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <cfenv>
#include <fstream>
#include <complex_bessel.h>


using namespace std;
using namespace itensor;

void initial_tensor_Ising(ITensor &T, complex<double> beta);
void initial_tensor_XY(ITensor &T, complex<double> beta, complex<double> h, int dcut);
void hotrg(ITensor &T, complex<double> &log_fact, int Niter, int dcut);


//void heev_f(const char jobz, const char uplo, int n, dcomplex* a, int lda, double* w);

int main(int argc, char* argv[])
{
        int L = 7;      //L : size 2^L x 2^L
        int Lp = pow(4,L);      // number of sites
        int dcut = 100;  //dcut : Bond dimension cutoff

        cout.precision(14);

        Cplx K (1.1199 + 0.000*atoi(argv[1]),0);

        double hr = 0;
        double hi = atof(argv[2]);

        Cplx h0 (hr, 0.0);
        Cplx h (hr, hi);

        auto T = ITensor(0);
	initial_tensor_XY(T,K,h,dcut);
	complex<double> log_fact = 0;
        hotrg(T,log_fact,2*L,dcut);

        Cplx Z = log_fact;

        if(imag(Z)>=0)
                printf("%.16g+%.16gj\n",real(Z),imag(Z));
        else
                printf("%.16g%.16gj\n",real(Z),imag(Z));


        return 0;

}

void initial_tensor_Ising(ITensor &T, complex<double> beta)
{
   auto s = Index(2);
   auto x1 = addTags(s,"x1");
   auto y1 = addTags(s,"y1");
   auto x2 = addTags(s,"x2");
   auto y2 = addTags(s,"y2");

   auto al = Index(2);
   auto W = ITensor(s,al);
   W.set(1,1,sqrt(cosh(beta)));
   W.set(1,2,sqrt(sinh(beta)));
   W.set(2,1,sqrt(cosh(beta)));
   W.set(2,2,-sqrt(sinh(beta)));
   T = ITensor(x1,y1,x2,y2);
   for(auto s1: range1(2))
   for(auto s2: range1(2))
   for(auto s3: range1(2))
   for(auto s4: range1(2))
   {
      complex<double> P = 0;
      for(auto s5: range1(2)){
         P += eltC(W,s5,s1)*eltC(W,s5,s2)*eltC(W,s5,s3)*eltC(W,s5,s4);
      }

      T.set(s1,s2,s3,s4,P);
   }
   return;
}

void initial_tensor_XY(ITensor &T, complex<double> beta,complex<double> h, int dcut)
{
   auto s = Index(dcut-1);
   auto x1 = addTags(s,"x1");
   auto y1 = addTags(s,"y1");
   auto x2 = addTags(s,"x2");
   auto y2 = addTags(s,"y2");
   T = ITensor(x1,y1,x2,y2);

   complex<double> *bes_b = new complex<double>[dcut-1];
   complex<double> *bes_bh = new complex<double>[2*(dcut-1)-1];

   for(auto s1: range(dcut-1))
      bes_b[s1] = (sp_bessel::besselI(0.5*(dcut-2)-(double)s1,beta));

   for(auto s1: range(2*(dcut-1)-1))
      bes_bh[s1] = sp_bessel::besselI(s1,beta*h);


   complex<double> P = 0;
   for(auto s1: range(dcut-1))
   for(auto s2: range(dcut-1))
   for(auto s3: range(dcut-1))
   for(auto s4: range(dcut-1))
   {
      P = sqrt(bes_b[s1]*bes_b[s2]*bes_b[s3]*bes_b[s4])*bes_bh[abs(s1+s2-s3-s4)];
      T.set(s1+1,s2+1,s3+1,s4+1,P);
   }
   return;
}

void hotrg(ITensor &T, complex<double> &log_fact, int Niter, int dcut)
{
   auto x1 = findIndex(T,"x1");
   auto y1 = findIndex(T,"y1");
   auto x2 = findIndex(T,"x2");
   auto y2 = findIndex(T,"y2");
   int dxn = 0;
   log_fact = 0;
   complex<double> P;

   for(auto n: range1(Niter))
   {
      x1 = findIndex(T,"x1");
      x2 = findIndex(T,"x2");
      y1 = findIndex(T,"y1");
      y2 = findIndex(T,"y2");
      T = permute(T,{x1,y1,x2,y2});
      auto x3 = replaceTags(x1,"x1","x3");
      auto x4 = replaceTags(x1,"x1","x4");
      auto dx = dim(x1);
      auto dy = dim(y1);

      if(dx*dx > dcut)
         dxn = dcut;
      else
         dxn = dx*dx;

      auto T_p = T;

      T_p = T_p.conj();

//     Obtaining R_L

      auto M1 = replaceTags(prime(T,x2),"y2","i")*replaceTags(T_p,"y2","j");
      auto temp = replaceTags(prime(replaceTags(T,"x2","x4"),x4),"y1","i")*
                  replaceTags(replaceTags(T_p,"x2","x4"),"y1","j");
      M1 = M1*temp;
      temp = ITensor(0);

      auto [U1, D1] = diagHermitian(M1,{"Tags","u1"});
      auto u1 = commonIndex(U1,D1);
      auto u1p = prime(u1);
      auto xn1 = Index(dx*dx,"xn1");
      auto ind_RL = replaceTags(xn1,"xn1","RL");

      auto sD1 = ITensor(u1p,ind_RL);

      for(auto s1: range1(dx*dx))
         sD1.set(s1,s1,1.0);

      sD1 = D1*sD1;

      for(auto s1: range1(dx*dx))
      {
         P = sqrt(fabs(eltC(sD1,u1=s1,ind_RL=s1)));
         sD1.set(u1=s1,ind_RL=s1,P);
      }

      auto R_L = U1*sD1;
      U1 = ITensor(0);
      D1 = ITensor(0);
      sD1 = ITensor(0);


//     Obtaining R_R

      auto M2 = replaceTags(prime(T,x1),"y2","i")*replaceTags(T_p,"y2","j");
      temp = replaceTags(prime(replaceTags(T,"x1","x3"),x3),"y1","i")*
             replaceTags(replaceTags(T_p,"x1","x3"),"y1","j");
      M2 = M2*temp;
      temp = ITensor(0);

      auto [U2, D2] = diagHermitian(M2, {"Tags","u2"});
      auto u2 = commonIndex(U2,D2);
      auto u2p = prime(u2);
      auto xn2 = replaceTags(xn1,"xn1","xn2");
      auto ind_RR = replaceTags(xn2,"xn2","RR");
     
      auto sD2 = ITensor(u2p,ind_RR);
      for(auto s1: range1(dx*dx))
         sD2.set(s1,s1,1.0);

      sD2 = D2*sD2;

      for(auto s1: range1(dx*dx))
      {
         P = sqrt(fabs(eltC(sD2,u2=s1,ind_RR=s1)));
         sD2.set(u2=s1,ind_RR=s1,P);
      }

      auto R_R = U2*sD2;
      U2 = ITensor(0);
      D2 = ITensor(0);
      sD1 = ITensor(0);


//      SVD PHASE : A (ind_RL, ind_RR)
      auto A = replaceTags(replaceTags(R_L,"x2","x1"),"x4","x3")*R_R;
      auto [U, S, V] = svd(A, {ind_RL},{"MaxDim",dxn,"SVDMethod","gesvd","LeftTags","ulink","RightTags","vlink"});
      auto ul = commonIndex(U,S);
      auto vl = commonIndex(S,V);
      auto vlp = prime(vl);

      auto S_copy = ITensor(vl,vlp);

      dxn = dim(vl);

      for(auto s1: range1(dxn))
         S_copy.set(s1,s1,1.0);


      S_copy = S*S_copy;

      S = ITensor(0);

      S_copy = noPrime(S_copy);


      for(auto s1: range1(dxn))
      {
         P = 1.0/sqrt(eltC(S_copy,ul=s1,vl=s1));
         S_copy.set(ul=s1,vl=s1,P);
      }


      auto F1 = U*(S_copy); auto F2 = (S_copy)*V;

      F1 = replaceTags(F1,"vlink","link");
      F2 = replaceTags(F2,"ulink","link");


      // To correct it and thus conserve the symmetry, we insert dagger operation to PR.

      auto PL = dag(F2)*(replaceTags(replaceTags(R_R,"x1","x2"),"x3","x4"));
      auto PR = dag(dag(F1)*(replaceTags(replaceTags(R_L,"x2","x1"),"x4","x3")));


      auto it1 = findIndex(PL,"link");
      auto it2 = findIndex(PR,"link");
      auto itt = Index(dim(it1),"link");
      PL = PL * delta(it1,itt);
      PR = PR * delta(it2,itt);   
      PL.replaceTags("link","xn2");
      PR.replaceTags("link","xn1");

      xn1 = findIndex(PR,"xn1");
      xn2 = findIndex(PL,"xn2");

      auto Uf = ITensor(0);
      auto Tn = ITensor(xn1,y1,xn2,y2);
      auto Tp = ITensor(0);


/*
      Tn = PR*(replaceTags(T,"y2","i"));
      Tn = Tn*(replaceTags(replaceTags(replaceTags(T,"y1","i"),"x1","x3"),"x2","x4"));
      Tn = Tn*PL;
*/


      for(auto s3: range1(dxn))
      {
         Uf = ITensor(x1,x3);
         for(auto s1: range1(dx))
         for(auto s2: range1(dx))
            Uf.set(x1=s1,x3=s2,eltC(PR,x1=s1,x3=s2,xn1=s3));


         Tp = Uf*(replaceTags(T,"y2","i"));
         Tp = Tp*(replaceTags(replaceTags(replaceTags(T,"y1","i"),"x1","x3"),"x2","x4"));
         Tp = Tp*PL;

         Uf = ITensor(0);

         for(auto i1: range1(dy))
         for(auto i2: range1(dxn))
         for(auto i3: range1(dy))
         {
            P = eltC(Tp,y1=i1,xn2=i2,y2=i3);
            Tn.set(xn1=s3,y1=i1,xn2=i2,y2=i3,P);
         }

         Tp = ITensor(0);
      }

      T = Tn;
      Tn = ITensor(0);
      PL = ITensor(0);
      PR = ITensor(0);
 
      T = replaceTags(T,"y2","x1");
      T = replaceTags(T,"y1","x2");
      T = replaceTags(T,"xn1","y1");
      T = replaceTags(T,"xn2","y2");
      
      complex<double> Nr = norm(T);

      T = T/Nr;

      log_fact = log_fact + log(Nr)/pow(2.0,n);
   }

   x1 = findIndex(T,"x1");
   y1 = findIndex(T,"y1");
   x2 = findIndex(T,"x2");
   y2 = findIndex(T,"y2");

   auto dx = dim(x1);
   auto dy = dim(y1);

   complex<double> Tra = 0.0;

   for(auto s1: range1(dx))
   for(auto s2: range1(dy))
      Tra = Tra + eltC(T,x1=s1,x2=s1,y1=s2,y2=s2);

   Tra = log(Tra)/pow(2.0,Niter);
   log_fact = log_fact + Tra;

   return;
}


