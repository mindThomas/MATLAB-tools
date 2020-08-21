/* sample_ghmm.c 


  Returns samples of a  Multivariate Gaussian process driven by HMM.

 Usage  
 -------

  [Z , X]  = sample_ghmm(K , PI , A , M , S , [v1] , ... , [vp] );

 Inputs
 -------

  PI                         Initial state probabilities (d x 1 x [n1] , ... , [nl]) , sum(PI) = ones(1 , 1 , [n1] , ... , [nl] )
  A                          State transition probabilities (d x d x [n1] , ... , [nl]), sum(A) = ones(1 , d , [n1] , ... , [nl])
  M                          Mean vectors (m x 1 x d x [n1] , ... , [nl])
  S                          Covariance matrices (m x m x d x [n1] , ... , [nl])

 Ouputs
 -------

  Z                          Measurements (m x K x [n1] , ... , [nl] , [v1] , ... , [vp]) 
  X                          State/label sequence (1 x K x [n1] , ... , [nl] , [v1] , ... , [vp])
 		
		  
  Example
  -------

   d       = 3;
   m       = 2;
   N       = 1;
   K       = 3;
   L       = 4;
 
   PI      = rand(d , 1 , N);
   sumPI   = sum(PI);
   PI      = PI./sumPI(ones(d , 1) , : , :);


   A       = rand(d , d);
   sumA    = sum(A);
   A       = A./sumA(ones(d , 1) , : , :);


   M       = 3*randn(m , 1 , d);
   S       = cat(3 , [2 , 1.1 ; 1.1 1.5] , [1.3 , 0.3 ; 0.3 1.5] ,  [1 , 0 ; 0 1]);


   %[Z , X] = sample_ghmm(K , PI(: , : , ones(1 , 4) ,  ones(1 , 2)) , A(: , : , ones(1 , 4), ones(1 , 2)) , M(: , : , : , ones(1 , 4), ones(1 , 2)) , S(: , : , : , ones(1 , 4), ones(1 , 2)) );


   [Z , X] = sample_ghmm(K , PI , A , M , S , L);


   %[Z , X] = sample_ghmm(K , PI(: , : , ones(1,L)) , A(: , : , ones(1,L)) , M(: , : , : , ones(1,L)) ,  S(: , : , : , ones(1,L)));

  
 
				
  To compile
  ----------
  
  mex  -DranSHR3 sample_ghmm.c

  
  Myself, I use Intel CPP compiler as : 

  
  mex -DranSHR3 -f mexopts_intel10.bat sample_ghmm.c

  Ver 1.0 (03/16/05)


  Author : Sébastien PARIS  © (sebastien.paris@lsis.org)
  --------
  

*/


#include <math.h>
#include <time.h>
#include "mex.h"

/*---------------- Basic generators definition ------------------- */

#define mix(a , b , c) \
{ \
	a -= b; a -= c; a ^= (c>>13); \
	b -= c; b -= a; b ^= (a<<8); \
	c -= a; c -= b; c ^= (b>>13); \
	a -= b; a -= c; a ^= (c>>12);  \
	b -= c; b -= a; b ^= (a<<16); \
	c -= a; c -= b; c ^= (b>>5); \
	a -= b; a -= c; a ^= (c>>3);  \
	b -= c; b -= a; b ^= (a<<10); \
	c -= a; c -= b; c ^= (b>>15); \
}

#define zigstep 128 /*  Number of Ziggurat'Steps */
#define znew   (z = 36969*(z&65535) + (z>>16) )
#define wnew   (w = 18000*(w&65535) + (w>>16) )
#define MWC    ((znew<<16) + wnew )
#define SHR3   ( jsr ^= (jsr<<17), jsr ^= (jsr>>13), jsr ^= (jsr<<5) )
#define CONG   (jcong = 69069*jcong + 1234567)
#define KISS   ((MWC^CONG) + SHR3)

#ifdef ranKISS
#define randint KISS
#define rand() (randint*2.328306e-10)
#endif 

#ifdef ranSHR3
#define randint SHR3
#define rand() (0.5 + (signed)randint*2.328306e-10)
#endif 

/*--------------------------------------------------------------- */

#ifdef __x86_64__
    typedef int UL;
#else
    typedef unsigned long UL;
#endif

/*--------------------------------------------------------------- */
static UL jsrseed = 31340134 , jsr;
#ifdef ranKISS
static UL z=362436069, w=521288629, jcong=380116160;
#endif

static UL iz , kn[zigstep];
static long hz;
static double wn[zigstep] , fn[zigstep];
/*--------------------------------------------------------------- */

void randini(void);  
void randnini(void);
double nfix(void);
double randn(void); 
void  matvect(double * , double * , double *, int , int , int); 
void  chol(double * , double * , int  , int); 
void  sample_ghmm(int , double * , double * , double * , double * , int , int , double * , double * , int , int , double * , double * , double *); 

/*--------------------------------------------------------------- */

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *PI , *A , *M , *S ;
	double *Z , *X;
	double *choles , *v , *b;
	const int *dimsPI , *dimsA , *dimsM , *dimsS;
	int    *dimsZ;
	int    numdimsPI , numdimsA , numdimsM , numdimsS;
	int    numdimsZ;
	int    K , d  , m , i , N = 1 , V = 1;
		
	/* Check input */
	
	if(nrhs < 5)	
	{     	
		mexErrMsgTxt("At least 5 inputs argument are required for sample_ghmm");		
	}
	
	/* Input 1 */
	if(mxGetNumberOfElements(prhs[0]) != 1)	
	{
		mexErrMsgTxt("First Input must be a scalar");	
	}
	
	K             = (int) mxGetScalar(prhs[0]);
	
	/* Input 2 */
		
	PI            = mxGetPr(prhs[1]);	
	numdimsPI     = mxGetNumberOfDimensions(prhs[1]);	
	dimsPI        = mxGetDimensions(prhs[1]);
		
	d             = dimsPI[0];	
	if (dimsPI[1] != 1)		
	{		
		mexErrMsgTxt("Second input must be (d x 1 x [n1] x ... x [nl])");			
	}
			
	for (i = 2 ; i < numdimsPI ; i++)		
	{		
		N  *=dimsPI[i];		
	}
	
	
	/* Input 3 */
		
	A             = mxGetPr(prhs[2]);	
	numdimsA      = mxGetNumberOfDimensions(prhs[2]);	
	dimsA         = mxGetDimensions(prhs[2]);
		
	/* Input 4 */
	
	M             = mxGetPr(prhs[3]);	
	numdimsM      = mxGetNumberOfDimensions(prhs[3]);	
	dimsM         = mxGetDimensions(prhs[3]);

	if ((dimsM[1] != 1) || (dimsM[2] != d))		
	{		
		mexErrMsgTxt("Fourth input must be (m x 1 x d x [n1] x ... x [nl])");			
	}
	m             = dimsM[0];

	/* Input 5 */
		
	S             = mxGetPr(prhs[4]);	
	numdimsS      = mxGetNumberOfDimensions(prhs[4]);	
	dimsS         = mxGetDimensions(prhs[4]);	
	numdimsZ      = 2 + (numdimsM - 3) + (nrhs - 5);	
	dimsZ         = (int *)malloc(numdimsZ*sizeof(int));	
	dimsZ[0]      = m;	
	dimsZ[1]      = K;

	for(i = 3 ; i < numdimsM ; i++)	
	{	
		dimsZ[i - 1] = dimsM[i];	
	}
	for (i = 5 ; i < nrhs ; i++)		
	{		
        dimsZ[(numdimsM - 3) + i  - 3  ] = (int) mxGetScalar(prhs[i]) ;		
		V                               *= dimsZ[(numdimsM - 3) + i  - 3  ];		
	}

	/* Output 1 */
	
	plhs[0]        = mxCreateNumericArray(numdimsZ , dimsZ, mxDOUBLE_CLASS, mxREAL);
	Z              = mxGetPr(plhs[0]);
	
	dimsZ[0]       = 1;	
	plhs[1]        = mxCreateNumericArray(numdimsZ, dimsZ, mxDOUBLE_CLASS, mxREAL);
	X              = mxGetPr(plhs[1]);
	
    choles         = (double *)malloc((m*m*d*N)*sizeof(double)); 
    v              = (double *)malloc(m*sizeof(double)); 
    b              = (double *)malloc(m*sizeof(double)); 
	
	
	/*---------------------------------------------------------------*/
	/*------------------------ MAIN CALL ----------------------------*/
	/*---------------------------------------------------------------*/
	
	randini();
	randnini();

	sample_ghmm(K , PI , A , M , S , N , V , Z , X , d , m , choles , v , b);
	
	/*---------------------------------------------------------------*/
	/*------------------------ FREE MEMORY --------------------------*/
	/*---------------------------------------------------------------*/

	free(dimsZ);
	free(choles);
	free(v);
	free(b); 	 
 }
 
 /* ----------------------------------------------------------------------- */
 void  sample_ghmm(int K , double *PI , double *A , double *M, double *S , int N , int V , double *Z , double *X , int d , int m , double *choles , double *v , double *b) 
 {	 
	 int l , i , s ,  k;	 
	 int lK, lmK , km , ld , lmd , KN = K*N , mKN = m*KN , sKN , smKN , lmmd , ldd , Xprevious  , Xcurrent , mXprevious , dXprevious , m2 = m*m;	 
	 register double temp , cP;

	 chol(S , choles , m , d*N); 

	 for (s = 0 ; s < V ; s++)		 
	 {		 		 
		 sKN  = s*KN;
		 smKN = m*sKN;
		 for (l = 0 ; l < N ; l++)			 
		 {			 
			 lK         = l*K + sKN;			 
			 lmK        = m*l*K + smKN;	 			 
			 ld         = l*d;			 
			 lmd        = ld*m;			 
			 lmmd       = m*lmd;			 
			 ldd        = ld*d;

			 Xprevious  = 1;
			 cP         = PI[0 + ld];
			 temp       = rand();
			 while ( (temp > cP) && (Xprevious < d) ) 
			 {  	
				 cP        += PI[Xprevious + ld];
				 Xprevious++;		
			 }  
			 X[0 + lK]  = Xprevious;
			 mXprevious = (Xprevious - 1)*m;
			 for (i = 0; i < m ; i++)
			 {
				 v[i] = randn();
			 }
			 matvect(choles , v , b , m , m , (Xprevious - 1)*m2 + lmmd); 
			 for (i = 0; i < m ; i++)
			 {
				 Z[i + lmK] = b[i] + M[i + mXprevious + lmd];
			 }

			 for (k = 1 ; k < K ; k++)
			 {		 
				 dXprevious = d*Xprevious;
				 Xcurrent   = 1;
				 cP         = A[0 + dXprevious + ldd];
				 temp       = rand();
				 while ( (temp > cP) && (Xcurrent < d) ) 
				 {  	
					 cP        += A[Xcurrent + dXprevious + ldd];
					 Xcurrent++;		
				 }  
				 X[k + lK]  = Xcurrent;
				 Xprevious  = (Xcurrent - 1);
				 mXprevious = Xprevious*m;
				 for (i = 0; i < m ; i++)
				 {
					 v[i] = randn();
				 }
				 matvect(choles , v , b , m , m , Xprevious*m2 + lmmd); 
				 km       = k*m + lmK;	 
				 for (i = 0 ; i < m ; i++)			 
				 {			 
					 Z[i + km] = b[i] + M[i + mXprevious + lmd];				 
				 }			 
			 }
		 }	  
	 }
 }
 /*----------------------------------------------------------*/
 void matvect(double *A , double *v , double *w, int d , int n , int off) 
/*
   w = Av, A(d x n), v(n x 1)   
 */   	   
 { 
	 int t , i ;	 
	 register double temp;	 	 
	 for (t = 0 ; t < d ; t++)
	 {
		 temp   = 0.0;
		 for(i = 0 ; i < n ; i++)
		 {				
			 temp += A[t + i*d + off]*v[i];
		 }
		 w[t] = temp;
	 }
 }
 /*----------------------------------------------------------*/
 void chol(double *Q , double *D , int d , int M) 	 
 {	 
	 int i , j , r , d2=d*d; 	 
	 int id , d1 = d - 1 , i1 , l , knnn , jd , v;
	 double sum , p , inv_p;

	 for (r = 0 ; r  < M ; r++)
	 {
		 v = r*d2;
		 for (i = 0 ; i < d2 ; i++)
		 {
			 D[i + v]    = Q[i + v];
		 }
		 p           = sqrt(D[0 + v]);
		 inv_p       = 1.0/p;
		 D[0 + v]    = p;

		 for(i = 1 ; i < d; i++)
		 {
			 D[d*i + v]  *= inv_p;
		 }

		 for(i = 1 ; i < d; i++)
		 {
			 id   = i*d;
			 i1   = i - 1;
			 sum  = D[i + id + v];    /* sum = B[i][i] */
			 for(l = 0; l < i; ++l)
			 {			
				 knnn = id + l;
				 sum -= D[knnn + v]*D[knnn + v];
			 }

			 p     = sqrt(sum);
			 inv_p = 1.0/p;
			 for(j = d1; j > i ; --j)
			 {
				 jd   = j*d;				 
				 sum  = D[jd + i + v];
				 for(l = 0; l < i ; ++l)
				 {				
					 sum   -= D[jd + l + v]*D[id + l + v];
				 }
				 D[jd + i + v] = sum*inv_p;
			 }

			 D[i + id + v] = p;
			 for(l = d1  ; l>i1 ; l--)
			 {
				 D[l + i1*d + v] = 0.0;	 
			 }
		 }

		 /* D = D'; */

		 for (j = 0 ; j < d ; j++)
		 {			
			 jd = j*d;			 
			 for(i = j + 1 ; i < d ; i++)				 
			 {
				 D[i + jd + v]  = D[j + i*d + v];
				 D[j + i*d + v] = 0.0;
			 }
		 }
	 }
 }
 /* --------------------------------------------------------------------------- */
 void randini(void) 
 {
	 /* SHR3 Seed initialization */

	 jsrseed  = (UL) time( NULL );
	 jsr     ^= jsrseed;
	 /* KISS Seed initialization */
#ifdef ranKISS
	 z        = (UL) time( NULL );
	 w        = (UL) time( NULL ); 
	 jcong    = (UL) time( NULL );
	 mix(z , w , jcong);

#endif 
 }

 /* --------------------------------------------------------------------------- */
 void randnini(void) 
 {	  
	 register const double m1 = 2147483648.0, m2 = 4294967296.0 ;
	 register double  invm1;
	 register double dn = 3.442619855899 , tn = dn , vn = 9.91256303526217e-3 , q; 
	 int i;

	 /* Ziggurat tables for randn */	 

	 invm1             = 1.0/m1;
	 q                 = vn/exp(-0.5*dn*dn);  
	 kn[0]             = (dn/q)*m1;	  
	 kn[1]             = 0;
	 wn[0]             = q*invm1;	  
	 wn[zigstep - 1 ]  = dn*invm1;
	 fn[0]             = 1.0;	  
	 fn[zigstep - 1]   = exp(-0.5*dn*dn);		
	 for(i = (zigstep - 2) ; i >= 1 ; i--)      
	 {   
		 dn              = sqrt(-2.*log(vn/dn + exp(-0.5*dn*dn)));          
		 kn[i+1]         = (dn/tn)*m1;		  
		 tn              = dn;          
		 fn[i]           = exp(-0.5*dn*dn);          
		 wn[i]           = dn*invm1;      
	 }
 }
/* --------------------------------------------------------------------------- */
double nfix(void)
{
    const double r = 3.442620; 	/* The starting of the right tail */
    static double x, y;
    
    for(;;)
    {
        
        x = hz*wn[iz];
        if(iz == 0)
        {	/* iz==0, handle the base strip */
            do
            {
                x = -log(rand())*0.2904764;  /* .2904764 is 1/r */
                y = -log(rand());
            }
            while( (y + y) < (x*x));
            return (hz > 0) ? (r + x) : (-r - x);
        }
        if( (fn[iz] + rand()*(fn[iz-1] - fn[iz])) < ( exp(-0.5*x*x) ) )
        {
            return x;
        }
        hz = randint;
        iz = (hz & (zigstep - 1));
        if(abs(hz) < kn[iz])
        {
            return (hz*wn[iz]);
        }
    }
}

 /* --------------------------------------------------------------------------- */
 double randn(void) 
 { 
	 hz = randint;
	 iz = (hz & (zigstep - 1));
	 return (abs(hz) < kn[iz]) ? (hz*wn[iz]) : ( nfix() );

 }
/* --------------------------------------------------------------------------- */




