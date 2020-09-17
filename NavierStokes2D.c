#include <stdlib.h>
#include <stdio.h>

//#include "fields.h"
//#include "blas.h"
//#include "real.h"

#include <math.h>
#define MIN(a,b) (((a)<(b))?(a):(b))

double** generateArray2D(int n,int m)
{
   int i;
   double** A = calloc(n, sizeof(double*));
   for (i=0;i<n;i++)
   {
      A[i] = calloc(m, sizeof(double));
   }
   return A;
}

double FIND_ABSMAX(double** A,int imax,int jmax)
{
   int i,j;
   double absMaxA, tmp;
   for(i=0;i<=imax+1;i++)
   {
      for(j=0;j<=jmax+1;j++)
      {
         tmp = fabs(A[i][j]);
	 if (tmp > absMaxA)
	 {
             absMaxA = tmp;
	 }
      }
   }
   return absMaxA;
}


void INIT_UVP(double** U,double** V,double** P,int imax,int jmax,double UI,double VI,double PI)
{
   int i,j;
   for(i=0;i<=imax+1;i++)
   {
      for(j=0;j<=jmax+1;j++)
      {
         U[i][j]=UI;
         V[i][j]=VI;
         P[i][j]=PI;
      }
   }
  return;
}


void COMP_DELT(double** U,double** V,int imax,int jmax,double* delt,double delx,double dely,double Re, double tau)
{
   double threshold = 1.0e-10;
   double tmp1 = 0.5*Re/( 1.0/(delx*delx) + 1.0/(dely*dely) );
   double tmp2, tmp3;
   double uAbsMax = FIND_ABSMAX(U,imax,jmax);
   double vAbsMax = FIND_ABSMAX(V,imax,jmax);
   if (uAbsMax < threshold)
   {
      uAbsMax = threshold;
   }
   if (vAbsMax < threshold)
   {
      vAbsMax = threshold;
   }
   tmp2 = delx/uAbsMax;
   tmp3 = dely/vAbsMax;

   if (tau>0)
   {
      *delt = tau*MIN(MIN(tmp1,tmp2),tmp3);
   }
   return;
}

void SETBCOND(double** U,double** V,int imax,int jmax,int wW,int wE,int wN,int wS)
{
   int i,j;
   for (j=0;j<=jmax+1;j++)
   {  
      //west
      U[     0][j]=        0.0;
      V[     0][j]=-V[   1][j];
      //east
      U[imax  ][j]=        0.0;
      V[imax+1][j]=-V[imax][j];
   }

   for (i=0;i<=imax+1;i++)
   {
      //north
      V[i][jmax  ]=        0.0;
      U[i][jmax+1]=-U[i][jmax];
      //south
      V[i][     0]=        0.0;
      U[i][     0]=-U[i][   1];
   }

}

void SETSPECBCOND(double** U,double** V,int imax,int jmax)
{
   int i;
   for (i=0;i<=imax+1;i++)
   {
      U[i][jmax+1]=2.0 - U[i][jmax];
   }
}

void COMP_FG(double** U,double** V,double** F,double** G,int imax,int jmax,double delt,double delx,double dely,double GX,double GY,double gamma,double Re)
{
   int i,j;

   double UXX,UYY,U00pp0,Um0p00,U00mp0,Um0m00,U2X,U00p0p,U0mp00,U00m0p,U0mm00,V00pp0,V0mppm,UVY;
   double VXX,VYY,V00p0p,V0mp00,V00m0p,V0mm00,V2Y,Vm0p00,V00mp0,Vm0m00,Um0pmp,UVX;
   
   printf("DELT  %f\n",delt);
   //3.36 and 3.37
   for (i=1;i<=imax-1;i++)
   {
       for (j=1;j<=jmax;j++)
       {
	    UXX = (U[i+1][j  ]-2.0*U[i][j]+U[i-1][j  ])/delx/delx;
            UYY = (U[i  ][j+1]-2.0*U[i][j]+U[i  ][j-1])/dely/dely;

            U00pp0 = U[i  ][j]+U[i+1][j];
            Um0p00 = U[i-1][j]+U[i  ][j];
            U00mp0 = U[i  ][j]-U[i+1][j];
            Um0m00 = U[i-1][j]-U[i  ][j];
            U2X    = ((pow(U00pp0,2.0)-pow(Um0p00,2.0)) + gamma*(fabs(U00pp0)*U00mp0-fabs(Um0p00)*Um0m00))/(4.0*delx);

            U00p0p = U[i][  j  ]+U[i  ][j+1];
            U0mp00 = U[i][  j-1]+U[i  ][j  ];
            U00m0p = U[i][  j  ]-U[i  ][j+1];
            U0mm00 = U[i][  j-1]-U[i  ][j  ];
            V00pp0 = V[i][  j  ]+V[i+1][j  ];
            V0mppm = V[i][  j-1]+V[i+1][j-1];
            UVY    = ((V00pp0*U00p0p-V0mppm*U0mp00) + gamma*(fabs(V00pp0)*U00m0p-fabs(V0mppm)*U0mm00))/(4.0*dely);

            F[i][j] = U[i][j] + delt*((UXX+UYY)/Re - U2X-UVY+GX);
	}
   }

   for (i=1;i<=imax;i++)
   {
       for (j=1;j<=jmax-1;j++)
       {
            VXX = (V[i+1][j  ]-2.0*V[i][j]+V[i-1][j  ])/delx/delx;
            VYY = (V[i  ][j+1]-2.0*V[i][j]+V[i  ][j-1])/dely/dely;

            V00p0p = V[i  ][j  ]+V[i  ][j+1];
            V0mp00 = V[i  ][j-1]+V[i  ][j  ];
            V00m0p = V[i  ][j  ]-V[i  ][j+1];
            V0mm00 = V[i  ][j-1]-V[i  ][j  ];
            V2Y    = ((pow(V00p0p,2.0)-pow(V0mp00,2.0)) + gamma*(fabs(V00p0p)*V00m0p-fabs(V0mp00)*V0mm00))/(4.0*dely);

            Vm0p00 = V[i-1][j  ]+V[i  ][j  ];
            V00mp0 = V[i  ][j  ]-V[i+1][j  ];
            Vm0m00 = V[i-1][j  ]-V[i  ][j  ];
            V00pp0 = V[i  ][j  ]+V[i+1][j  ];
            U00p0p = U[i  ][j  ]+U[i  ][j+1];
            Um0pmp = U[i-1][j  ]+U[i-1][j+1];
            UVX    = ((U00p0p*V00pp0-Um0pmp*Vm0p00) + gamma*(fabs(U00p0p)*V00mp0-fabs(Um0pmp)*Vm0m00))/(4.0*delx);

            G[i][j] = V[i][j] + delt*((VXX+VYY)/Re - UVX-V2Y+GY);
       }
   }

   // (3.42)
   for (i=1;i<=imax;i++)
   {
       G[i][   0] = V[i][   0];
       G[i][jmax] = V[i][jmax];
   }
   for (j=1;j<=jmax;j++)
   {
       F[0   ][j] = U[0   ][j];
       F[imax][j] = U[imax][j];
   }
   return;
}

void COMP_RHS(double** F,double** G,double** RHS,int imax,int jmax,double delt,double delx,double dely)
{
   int i,j;
   for (i=1;i<=imax;i++)
   {
      for (j=1;j<=jmax;j++)
      {
	 //3.38
         RHS[i][j] = ((F[i][j]-F[i-1][j])/delx + (G[i][j]-G[i][j-1])/dely)/delt;
         //printf("%f ",RHS[i][j]);
      }
   }
   return;
}

void POISSON(double** P, double** RHS,int imax,int jmax,double delx,double dely,double eps,int itermax,double omg,double* res)
{
   int i,j,it;
   int epsN,epsW,epsS,epsE;
   double tmp1a,tmp1b,tmp2a,tmp2b,tmp2c;
   double rit;
   double delx2 = delx*delx;
   double dely2 = dely*dely;
   
   for (it=1;it<=itermax;it++)
   {
      for (i=1;i<=imax;i++)
      {
         for (j=1;j<=jmax;j++)
         {
            if (i==   1) {epsW=0;} else {epsW=1;}
            if (i==imax) {epsE=0;} else {epsE=1;}
            if (j==   1) {epsN=0;} else {epsN=1;}
            if (j==jmax) {epsS=0;} else {epsS=1;}
            //3.44
	    tmp1a = (epsE             + epsW            )/delx2;
	    tmp1b = (epsN             + epsS            )/dely2;
            tmp2a = (epsE*P[i+1][j  ] + epsW*P[i-1][j  ])/delx2;
	    tmp2b = (epsN*P[i  ][j+1] + epsS*P[i  ][j-1])/dely2;
	    tmp2c = RHS[i][j];

	    P[i][j] = (1.0-omg)*P[i][j] + omg/(tmp1a+tmp1b)*(tmp2a + tmp2b - tmp2c);
	    
	    printf("%f  ",P[i][j]);
	    
         }
	 printf("\n");
      }
      //3.48
      for (i=1;i<=imax;i++)
      {     
        P[i][     0] = P[i][   1];
        P[i][jmax+1] = P[i][jmax];
      }
      for (j=1;j<=jmax;j++)
      {     
        P[0][     j] = P[1][   j];
        P[imax+1][j] = P[imax][j];
      }
      P[     0][     0] = P[   1][   1];
      P[     0][jmax+1] = P[   0][jmax];
      P[imax+1][     0] = P[imax][   0];
      P[imax+1][jmax+1] = P[imax][jmax];

      //3.46
      *res = 0.0;
      for (i=1;i<=imax;i++)
      {
         for (j=1;j<=jmax;j++)
         {

            if (i==   1) {epsW=0;} else {epsW=1;}
            if (i==imax) {epsE=0;} else {epsE=1;}
            if (j==   1) {epsN=0;} else {epsN=1;}
            if (j==jmax) {epsS=0;} else {epsS=1;}
            tmp2a = (epsE*(P[i+1][j  ]-P[i][j])-epsW*(P[i][j]-P[i-1][j  ]))/delx2;
	    tmp2b = (epsN*(P[i  ][j+1]-P[i][j])-epsS*(P[i][j]-P[i  ][j-1]))/dely2;
	    tmp2c = RHS[i][j];
	    rit   = tmp2a + tmp2b - tmp2c;
	    *res += rit*rit;
	    //printf("   rit_%d_%d  %f %f %d %d %d %d\n",i,j,rit,RHS[i][j],epsN,epsS,epsW,epsE);
         }
      }
      *res = sqrt((*res)/(imax*jmax));
      printf("    %d %f\n",it,*res);
   }
   return;
}


void ADAP_UV(double** U,double** V,double** F,double** G,double** P,int imax,int jmax,double delt,double delx,double dely)
{
   int i,j;
   //3.34
   for (i=1;i<=imax-1;i++)
   {
      for (j=1;j<=jmax;j++)
      {
         U[i][j] = F[i][j] - (P[i+1][j]-P[i][j])*delt/delx;
      }
   }
   //3.35
   for (i=1;i<=imax;i++)
   {
      for (j=1;j<=jmax-1;j++)
      {
         V[i][j] = G[i][j] - (P[i][j+1]-P[i][j])*delt/dely;
      }
   }
}

void WRITE_UVP(double** U,double** V,double** P,int imax,int jmax)
{
   int i,j,k;
   FILE * fp;
   fp = fopen("/Users/christiant/MyDocuments/Programming/C/UVP.dat","w");
   k = 8;
   for (i=0;i<=imax+1;i++)
   {
      for (j=0;j<=jmax+1;j++)
      {
        if ( i%k == 0 && (j%k==0 || j>=110))
        {
           fprintf(fp, "%d %d %f %f %f\n",i,j,U[i][j],V[i][j],P[i][j]);
        }
      }
   }
   fclose(fp);
   return;
}

int main(int argc, char* argv[])
{
   int n;
   int i,j,k;
   double res;
   
   int    imax    = 20; //128;
   int    jmax    = 20; //128;
   double xlength = 1.0;
   double ylength = 1.0;
   double delt    = 0.01; //0.02;
   double t       = 0.0;
   double t_end   = 5000; //10.0;
   double tau     = 0.5;
   int    itermax = 5;  //10000;    //100;
   double eps     = 1e-9; //0.001;
   double omg     = 1.7;
   double gamma   = 0.9;
   double Re      = 985.7; //1000.0;
   double UI      = 0.0;
   double VI      = 0.0;
   double PI      = 0.0;
   double GX      = 0.0;
   double GY      = 0.0;
   int    wW      = 2;
   int    wE      = 2;
   int    wS      = 2;
   int    wN      = 2;


   double delx = xlength/imax;
   double dely = ylength/jmax;

   double** U   = generateArray2D(imax+2,jmax+2);
   double** V   = generateArray2D(imax+2,jmax+2);
   double** P   = generateArray2D(imax+2,jmax+2);
   double** RHS = generateArray2D(imax+2,jmax+2);
   double** F   = generateArray2D(imax+2,jmax+2);
   double** G   = generateArray2D(imax+2,jmax+2);
   FILE * fp;
   fp = fopen ("/Users/christiant/MyDocuments/Programming/C/UVP.dat","w");


   INIT_UVP(U,V,P,imax,jmax,UI,VI,PI);
   n=0;
   while ((t<t_end) && (n<1))
   {
      COMP_DELT(   U,V,          imax,jmax,&delt,delx,dely,            Re,tau);
      SETBCOND(    U,V,          imax,jmax,                                   wW,wE,wN,wS);
      SETSPECBCOND(U,V,          imax,jmax);
      COMP_FG(     U,V,F,G,      imax,jmax, delt,delx,dely,GX,GY,gamma,Re);
      COMP_RHS(        F,G,  RHS,imax,jmax, delt,delx,dely);
      POISSON(             P,RHS,imax,jmax,      delx,dely,eps,itermax,omg,&res);
      ADAP_UV(     U,V,F,G,P,    imax,jmax, delt,delx,dely);
      WRITE_UVP(   U,V,    P,    imax,jmax);
      t += delt;
      n += 1;
      printf("%d %f %f \n",n, t, delt);
   }

   printf("---Done---");
   return 0;
}
