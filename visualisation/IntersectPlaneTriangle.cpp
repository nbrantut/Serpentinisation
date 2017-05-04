//==============================================================
//
// N : Normal to plane , |N| = 1, N->(nx3)
// Q : Point into the plane, Q->(nx3)
//
// P1,P2,P3 : Triangle vertices, P1,P2,P3->(nx3)
//
// return : 0 if no intersection
// U : Intersection points, U->()
//
//==============================================================

int Intersection(double *Q,double *N,double *P1,double *P2,double *P3,double *U){

    double d1 = (P1[0]-Q[0])*N[0] + (P1[1]-Q[1])*N[1] + (P1[2]-Q[2])*N[2] ;// (P1-Q) dot N
    double d2 = (P2[0]-Q[0])*N[0] + (P2[1]-Q[1])*N[1] + (P2[2]-Q[2])*N[2] ;// (P2-Q) dot N
    double d3 = (P3[0]-Q[0])*N[0] + (P3[1]-Q[1])*N[1] + (P3[2]-Q[2])*N[2] ;// (P3-Q) dot N

    bool s1 = d1<0 ;
    bool s2 = d2<0 ;
    bool s3 = d3<0 ;

    if(s2!=s1){
        U[0] = (d2*P1[0] - d1*P2[0])/(d2-d1) ;
        U[1] = (d2*P1[1] - d1*P2[1])/(d2-d1) ;
        U[2] = (d2*P1[2] - d1*P2[2])/(d2-d1) ;
        if(s1!=s3){
            U[3] = (d1*P3[0] - d3*P1[0])/(d1-d3) ;
            U[4] = (d1*P3[1] - d3*P1[1])/(d1-d3) ;
            U[5] = (d1*P3[2] - d3*P1[2])/(d1-d3) ;
        } else {
            U[3] = (d2*P3[0] - d3*P2[0])/(d2-d3) ;
            U[4] = (d2*P3[1] - d3*P2[1])/(d2-d3) ;
            U[5] = (d2*P3[2] - d3*P2[2])/(d2-d3) ;
        }
    } else if(s3!=s1){
        U[0] = (d3*P1[0] - d1*P3[0])/(d3-d1) ;
        U[1] = (d3*P1[1] - d1*P3[1])/(d3-d1) ;
        U[2] = (d3*P1[2] - d1*P3[2])/(d3-d1) ;

        U[3] = (d2*P3[0] - d3*P2[0])/(d2-d3) ;
        U[4] = (d2*P3[1] - d3*P2[1])/(d2-d3) ;
        U[5] = (d2*P3[2] - d3*P2[2])/(d2-d3) ;
    } else
    return 0 ;

    return 1;
}
//==============================================================

#include "mex.h"
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *N, *Q, *P1, *P2, *P3, *U;
    double *flag;
    
    N = mxGetPr(prhs[0]);
    Q = mxGetPr(prhs[1]);
    P1 = mxGetPr(prhs[2]);
    P2 = mxGetPr(prhs[3]);
    P3 = mxGetPr(prhs[4]);
    
    plhs[1] = mxCreateDoubleMatrix(3, 2, mxREAL);
    U =  mxGetPr(plhs[1]);
    

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    flag =  mxGetPr(plhs[0]);

    int out=Intersection(Q, N, P1, P2, P3, U);
    if (out==0) { *flag=0; }  else { *flag=2;}
}