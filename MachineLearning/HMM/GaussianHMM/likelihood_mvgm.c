/* likelihood_mvgm

  Compute the likelihood of a Multivariate Gaussian Mixture L(Z|x) where x = (M , S)

  Usage             L = likelihood_mvgm(Z , M , S , [P]);
  ------

  Inputs 
  ------

     Z              Measurements (m x K  x [v1] x ... x [vp]) 
     M              Mean vector (m x 1 x d x [n1] x ... x [nl]) 
     S              Covariance  (m x m x d x [n1] x ... x [nl] )
	 [P]            Weights (1 x 1 x d x [n1] x ... x [nl])  (default P = ones(1 , 1 , d x [n1] x ... x [nl]) )

  Outputs
  -------

     L              Likelihood of Z (d x K x [v1] x ... x [vp] x [n1] x ... x [nl])


  To compile
  ----------
  
  mex likelihood_mvgm.c
  
  Myself, I use Intel CPP compiler as : 

  mex  -f mexopts_intel10.bat likelihood_mvgm.c

		
  Example1
  -------

  d                                   = 2;                % dimension
  m                                   = 2;                % number of conpounds
  N                                   = 1000;             % number of training data


  P                                   = cat(3 , [0.4] , [0.6]);
  M                                   = cat(3 , [-1 ; -1] , [1 ; 1]);
  S                                   = cat(3 , [1 0.3 ; 0.3 0.8] , [0.7 0.6; 0.6 1]);
  [Z ,  X]                            = sample_mvgm(N , M , S , P);
  L                                   = likelihood_mvgm(Z , M , S , P);


  Example2
  -------


  d                                   = 2;                % dimension
  N                                   = 1000;             % number of training data


  M                                   = [-1 ; 1];
  S                                   = [1 0.3 ; 0.3 0.8] ;
  [Z ,  X]                            = sample_mvgm(N , M , S);
  L                                   = likelihood_mvgm(Z , M , S);


  

  Author : Sébastien PARIS (sebastien.paris@lsis.org) 
  -------		  

  
  Ver 1.1 (02/02/09) Correct a bug when d = 1.
  ------- 
				  
*/

#include <math.h>
#include "mex.h"

#define NUMERICS_FLOAT_MIN 1.0E-37
#define M_PI 3.14159265358979323846
#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif


/*--------------------------------------------------------------- */

double gauss(double *, double * , int , int);
void lubksb(double *, int , int *, double *);
int ludcmp(double *, int , int *, double * , double *);
double inv(double * , double * , double * , double * , int * , int );
void likelihood_mvgm(double * , double * , double * , double * , double * , int  , int  , int  , int  , int , double * , double * , double *, double * , double * , double * , int * , double *);

/*--------------------------------------------------------------- */

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *Z , *M , *S , *P;
	double *L;
	const int  *dimsZ , *dimsM , *dimsS;
	int  *dimsL , *indx;
	double *vect , *vv , *invsigma  , *temp_sigma , *temp_invsigma , *det_sigma , *res;
	int  numdimsZ , numdimsM , numdimsS;
	int  numdimsL;
	int  i , d = 1 , m , K , R=1 , V=1;
	int m2;
		
	/* Check input */	
	if(nrhs < 3)		
	{     
		mexErrMsgTxt(" At least 3 inputs argument are required for likelihood_mvgm");	
	}
		
	/* Input 1 */	
	Z               = mxGetPr(prhs[0]);	
	numdimsZ        = mxGetNumberOfDimensions(prhs[0]);	
	dimsZ           = mxGetDimensions(prhs[0]);	
	K               = dimsZ[1];	
	for (i = 2 ; i<numdimsZ ; i++)		
	{
		V *= dimsZ[i];	
	}
	
	/* Input 2 */	
	
	M               = mxGetPr(prhs[1]);	
	numdimsM        = mxGetNumberOfDimensions(prhs[1]);	
	dimsM           = mxGetDimensions(prhs[1]);	
	if ( (dimsM[1] != 1))		
	{
		mexErrMsgTxt("M must be (m x 1 x d)");	
	}	
	m               = dimsM[0];
	m2              = m*m;

	if(numdimsM > 2)
	{		
		d           = dimsM[2];		
		for(i = 3 ; i < numdimsM ; i++)			
		{			
			R         *= dimsM[i];
		}	
	}
		
	/* Input 3 */
		
	S              = mxGetPr(prhs[2]);	
	numdimsS       = mxGetNumberOfDimensions(prhs[2]);	
	dimsS          = mxGetDimensions(prhs[2]);
	
	if ( (dimsS[0] != m) && (dimsS[1] !=m) && (dimsS[2] != d))		
	{
		mexErrMsgTxt("S must be (m x m x d)");	
	}
	if(nrhs < 4)
	{		
		P    = (double *)mxMalloc(d*R*sizeof(double));
		for(i = 0 ; i < d*R ; i++)
		{
			P[i] = 1.0;
		}
	}
	else
	{
		P              = mxGetPr(prhs[3]);		
	}

	/* Output 1 */
	
	numdimsL      = 2 + max((numdimsZ - 2) , 0) + max((numdimsM - 3) , 0);
	dimsL         = (int *)malloc(numdimsL*sizeof(int));
	dimsL[0]      = d;
	dimsL[1]      = K;
	for (i = 2 ; i < numdimsZ ; i++)
	{
        dimsL[i] = dimsZ[i] ;
	}
	
	for (i = 3 ; i < numdimsM ; i++)	
	{	
        dimsL[ i - 1 + numdimsZ - 2] = dimsM[i] ;	
	}
	
	plhs[0]            = mxCreateNumericArray(numdimsL , dimsL , mxDOUBLE_CLASS, mxREAL);
	L                  = mxGetPr(plhs[0]);
	
	vect               = (double *)malloc(m*sizeof(double));
	vv                 = (double *)malloc(m*sizeof(double));
	temp_sigma         = (double *)malloc(m2*sizeof(double));
	temp_invsigma      = (double *)malloc(m2*sizeof(double));
	det_sigma          = (double *)malloc(d*sizeof(double));
	invsigma           = (double *)malloc((m2*d)*sizeof(double));
	indx               = (int *)malloc(m*sizeof(int));
	res                = (double *)malloc(m*sizeof(double));

		 
	/*---------------------------------------------------------------*/
	/*------------------------ MAIN CALL ----------------------------*/
	/*---------------------------------------------------------------*/
	
	likelihood_mvgm(Z , M , S , P , L , d , m , K , V , R , invsigma , temp_sigma , temp_invsigma , det_sigma , vect , vv , indx , res);
	
	/*---------------------------------------------------------------*/
	/*------------------------ FREE MEMORY --------------------------*/
	/*---------------------------------------------------------------*/
	
	free(vect);
	free(vv);
	free(temp_sigma);
	free(temp_invsigma);
	free(det_sigma);
	free(invsigma);
	free(indx);
	free(res);
	free(dimsL);	
 }
 
 /* ----------------------------------------------------------------------- */
 
 void likelihood_mvgm(double *Z , double *M , double *S , double *P , double *L , int d , int m , int K , int V , int R , double *invsigma , double *temp_sigma , double *temp_invsigma , double *det_sigma , double *vect , double *vv , int *indx , double *res) 
 {
	 
	 int  j , l , r , h , k, m2 = m*m , km , kd , lm2  , md = m*d , m2d = m2*d, rmd , rm2d;
	 int v , mK = m*K , dK = d*K , vmK , vdK , VdK = V*dK , rVdK , jm , rd , ind;
	 double cte = 1.0/pow(2*M_PI , m/2.0);
	  
	 for(r = 0 ; r < R ; r++)	 
	 {
		 rd   = r*d;
		 rmd  = r*md;
		 rm2d = r*m2d;
		 rVdK = r*VdK;
		 for (l = 0 ; l < d ; l++)		 
		 {		 
			 lm2   = l*m2;
			 ind   = lm2 + rm2d;		 
			 for(j = 0 ; j < m2 ; j++)			 
			 {				 
				 temp_sigma[j] = S[j + ind];				 
			 }
			 			 
			 det_sigma[l]  = inv(temp_sigma , temp_invsigma , vect , vv , indx , m);
			 			 
			 for(j = 0 ; j < m2 ; j++)				 
			 {				 
				 invsigma[j + lm2] = temp_invsigma[j];				 
			 }			 
			 det_sigma[l] = (cte*sqrt(fabs(det_sigma[l])));			 
		 }	 
		 for (v = 0 ; v < V ; v++) 	 
		 {	 
	     	 vmK = v*mK;	 
			 vdK = v*dK + rVdK;	 
			 for (k = 0 ; k < K ; k++)			 
			 {			 
    			 km    = k*m + vmK;			 
				 kd    = k*d + vdK;			 
				 for (j = 0 ; j < d ; j++)					 
				 {				 
					 jm            = j*m + rmd;				 
					 for(h = 0 ; h < m ; h++)				 
					 {	
						 res[h] = (Z[h + km] - M[h + jm]);	 
					 }
					 L[j + kd]  = P[j + rd]*det_sigma[j]*exp(-0.5*gauss(res , invsigma , m , j*m2));
				 }
			 }			 
		 }
	 }
 }
 
/*----------------------------------------------------------------------------------------------*/
double gauss(double *y, double *R , int d , int offset)
{	
	int  i , j , id;	
	register double temp , Q = 0.0;				
	id      = offset;
	for (i = 0 ; i < d ; i++)
	{		
		temp        = 0.0;		
		for(j = 0 ; j < d ; j++)			
		{		
			temp   += y[j]*R[j + id];		
		}	
		Q   += temp*y[i];
		id  += d;		
	}	
	return Q;	
}
/*------------------------------------------------------------------*/
double inv(double *temp , double *invQ  , double *vect , double *vv , int *indx , int d)
{
	int i , j , jd , d1 = d+1;	
	double dd , det = 1.0;
	if(ludcmp(temp , d , indx , &dd , vv ))
	{
		jd     = 0;
		for(j = 0 ; j < d ; j++)		
		{			
			det *= temp[jd];
			jd  += d1;		
		}
		jd    = 0;	
		for(j = 0; j < d; j++)
		{            
			for(i = 0; i < d; i++) 				
			{
				vect[i] = 0.0;
			}						
			vect[j] = 1.0;			
			lubksb(temp , d , indx , vect);			
			for(i = 0 ; i < d ; i++) 				
			{	
				invQ[jd + i] = vect[i];				
			}
			jd    += d;
		}		
	}
	return (1.0/det);	
}
/*-------------------------------------------------------------------------------*/
void lubksb(double *m, int n, int *indx, double *b)
{
    int i, ii = -1, ip, j , nn = n*n, in;
    double sum;
    for(i = 0; i < n; i++)
	{
        ip        = *(indx + i);
        sum       = *(b + ip);
        *(b + ip) = *(b + i);
        if(ii > -1)
		{
            for(j = ii; j <= i - 1; j++)
			{
                sum -= m[i + j*n] * *(b + j);
            }
        } 
		else if(sum)
		{
            ii = i;
        }
        *(b + i) = sum;
    }
    
	for(i = n - 1; i >= 0; i--)
	{
        sum = *(b + i);
		in  = i*n;
        for(j = i + 1; j < n; j++)
		{
            sum -= m[i + j*n] * *(b + j);
        }
		*(b + i) = sum / m[i + in];
    }
}

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
int ludcmp(double *m, int n, int *indx, double *d , double *vv)
{
    int i, imax, j, k , jn , kn , n1 = n - 1;
    double big, dum, sum , temp;
  
    d[0] = 1.0;
    for(i = 0; i < n; i++)
	{
        big = 0.0;
		jn  = 0;	
        for(j = 0; j < n; j++)
		{
            if((temp = fabs(m[i + jn])) > big)		
			{
                big = temp;
            }
			jn   += n;		
		}
        if(big == 0.0)
		{		
            return 0;
        }	
        vv[i] = 1.0 / big;
    }

	jn   = 0;
    for(j = 0; j < n; j++)
	{	
        for(i = 0; i < j; i++)	
		{
            sum = m[i + jn];
			kn    = 0;			
            for(k = 0 ; k < i; k++)				
			{
                sum -= m[i + kn ] * m[k + jn];
            }          
			m[i + jn] = sum;
			kn       += n;
        }
		
        big = 0.0;		
        for(i = j ; i < n; i++)	
		{
			kn  = 0;
            sum = m[i + jn];			
            for(k = 0; k < j; k++)		
			{            
				sum -= m[i + kn] * m[k + jn];
				kn  += n;         
			}          
			m[i + jn] = sum;			
            if((dum = vv[i] * fabs(sum)) >= big)				
			{
                big  = dum;				
                imax = i;
            }
        }      
		if(j != imax)		
		{
			kn    = 0;
            for(k = 0; k < n; k++)				
			{								
                dum           = m[imax + kn];				
                m[imax + kn]  = m[j + kn];				
                m[j + kn]     = dum;
				kn           += n;				
            }            
			d[0]       = -d[0];			
            vv[imax]   = vv[j];
        }
        indx[j] = imax;
        if(m[j + jn] == 0.0)
		{
			m[j + jn] = NUMERICS_FLOAT_MIN;
        
		}
		if(j != n1)
		{
            dum = 1.0 / (m[j + jn]);
            for(i = j + 1; i < n; i++)
			{
				m[i + jn] *= dum;
			}
        }
		jn    += n;
    }
    return 1;
};
