
/*

 Forward Backward HMM algorithm


 Usage  
 -------


 [X , gamma , logll , [sum_alpha] , [sum_beta]]  = forward_backward(PI , A , L , [filter]);


 Inputs
 -------

 PI            Initial proabilities (N x 1) : Pr(x_1 = i) , i=1,...,N

 A             State transition probabilities matrix Pr(x_{k} = i| x_{k - 1} = j) such
               sum_{x_k}(A) = 1 => sum(A , 1) = 1

 L             Time indexed Likelihood matrix Pr(z_k | x_k = i) (N x K), i=1,...,N, k=1,...,K. 
               Extracted from B matrix such that B = Pr(z | x) (M x N), sum(B , 1) = 1 and B(z_k , :)' = L(: , k).

 filter        Optional flag. If filter = 0 => Just the alpha probabilities are computed, else (default filter = 1)
               the two passes alpha & beta.


 Ouputs
 -------

 X             MAP sequence estimated from gamma probabilies Pr(x_k|Z_1:K)

 gamma         gamma probabilities (N x K) or alpha probabilities if filter = 0

 logll         Loglikelihood of (PI,A,L)

 If the -Dbetanormalize flag is used to compile, two more outputs are given :

 sum_alpha     Pr(z_k|Z_{1:k-1})  (1 x K)

 sum_beta      Pr(z_k|Z_{k+1:K})  (1 x K)


 To compile
 ----------


mex (-Dbetanormalize) -output forward_backward.dll forward_backward.c

The betanormalize flag forces to scale gamma's probabilities with the sum of beta instead of the alpha's.


mex -f mexopts_intel10.bat -output forward_backward.dll forward_backward.c

mex -Dbetanormalize -f mexopts_intel10.bat -output forward_backward.dll forward_backward.c


Example 1
----------


N                        = 256;
M                        = 256;
K                        = 300;

PI                       = rand(N , 1);
sumPI                    = sum(PI);
PI                       = PI./sumPI(ones(N , 1) , :);

A                        = rand(N , N);
sA                       = sum(A);
A                        = A./sA(ones(N , 1) , :);

L                        = rand(N , K);

[X , gamma , logll]      = forward_backward(PI , A , L);


Example 2
----------


N                        = 2;
M                        = 4;
K                        = 30;
PI                       = rand(N , 1);
sumPI                    = sum(PI);
PI                       = PI./sumPI(ones(N , 1) , :);

A                        = rand(N , N);
sumA                     = sum(A);
A                        = A./sumA(ones(N , 1) , :);

B                        = rand(M , N);
sumB                     = sum(B);
B                        = B./sumB(ones(M , 1) , :);

[Z , X]                  = sample_dhmm(K , PI , A , B);

L                        = densdis(Z , B); %Likelihood

[X_MAP , gamma , logll]  = forward_backward(PI , A , L);


Example3
----------


d                                      = 2;
M                                      = 256;
K                                      = 300;
W                                      = [1 ; M];
H                                      = 1;
R1                                     = 2;
R2                                     = 3;


PI                                     = [0.5 ; 0.5];
A                                      = [0.9 0.1 ; 0.1 0.9];
X                                      = sample_dmc(K , PI , A);


PI1                                    = iniPI(W);
A1                                     = matA_gauss(W , H , R1);
B1                                     = A1;


PI2                                    = iniPI(W);
A2                                     = matA_gauss(W , H , R2);
B2                                     = A2;


[Z1 , X1]                              = sample_dhmm(K , PI1 , A1 , B1);
[Z2 , X2]                              = sample_dhmm(K , PI2 , A2 , B2);

temp                                   = [Z1 ; Z2];
Z                                      = temp((0:K-1)*d + X);

L1                                     = densdis(Z , B1);
L2                                     = densdis(Z , B2);

[X1 , gamma1 , ll1 , calpha1 , cgamma1] = forward_backward(PI1 , A1 , L1);
[X2 , gamma2 , ll2 , calpha2 , cgamma2] = forward_backward(PI2 , A2 , L2);

L                                      = [calpha1.*cgamma1 ; calpha2.*cgamma2];
X_FB                                   = forward_backward(PI , A , L);


 Author : Sébastien PARIS  © (sebastien.paris@lsis.org)
 --------


*/

#include <math.h>
#include "mex.h"

#define NUMERICS_FLOAT_MIN 1.0E-37

/*--------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------------------*/

#ifdef betanormalize
void forward_backward(double * , double *  , double * , double * , double * , double *  , double * , int , int , int , double * , double * , double * , double *);
#else
void forward_backward(double * , double *  , double *  , double * , double * , double *  , double * , int , int , int , double * , double * , double *);
#endif
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
		
	double *PI , *A, *L ;
	double *X , *gamma , *logll; 
	double *alpha , *At , *temp_gamma;
	
#ifdef betanormalize
	double *sum_alpha , *sum_gamma;
#else
	double *invsum_alpha;
#endif
	
	const int *dimsPI , *dimsA , *dimsL ;
	int K , N , NK , numdimsPI , numdimsA  , numdimsL  , filter = 0;
	
	/*---------------------------------------------------------------*/
	/*---------------------- PARSE INPUT ----------------------------*/	
	/*---------------------------------------------------------------*/
	
	if( (nrhs < 3) | (nrhs >5))		
	{
		mexErrMsgTxt("3 or 4 input are requiered");	
	}
	PI         = mxGetPr(prhs[0]);
    numdimsPI  = mxGetNumberOfDimensions(prhs[0]);
	dimsPI     = mxGetDimensions(prhs[0]);
	if ( (numdimsPI>2) & (dimsPI[1] > dimsPI[0]) )
	{
		mexErrMsgTxt("PI must be (N x 1)");			 	
	}
	
	A         = mxGetPr(prhs[1]);
    numdimsA  = mxGetNumberOfDimensions(prhs[1]);
	dimsA     = mxGetDimensions(prhs[1]);
	if ( (numdimsA>2) & (dimsA[1] != dimsA[0]) )
	{
		mexErrMsgTxt("A must be (N x N)");			 	
	}
	
	L         = mxGetPr(prhs[2]);
    numdimsL  = mxGetNumberOfDimensions(prhs[2]);
	dimsL     = mxGetDimensions(prhs[2]);
	
	if ( (numdimsL>2) & (dimsL[0] != dimsA[0]) )
	{	
		mexErrMsgTxt("L must be (N x K)");			 	
	}
    N         = dimsL[0];
	K         = dimsL[1];
	NK        = N*K;
	if (nrhs == 4)
	{
		filter = (int)mxGetScalar(prhs[3]);
	}
	
	
	/*---------------------------------------------------------------*/
	/*---------------------- PARSE OUTPUT ---------------------------*/	
	/*---------------------------------------------------------------*/
	
	plhs[0]       = mxCreateDoubleMatrix(1 , K , mxREAL);
	X             = mxGetPr(plhs[0]);
	
	plhs[1]       = mxCreateDoubleMatrix(N , K , mxREAL);
	gamma         = mxGetPr(plhs[1]);
	
	plhs[2]       = mxCreateDoubleMatrix(1 , 1 , mxREAL);
	logll         = mxGetPr(plhs[2]);
	
#ifdef betanormalize
	plhs[3]       = mxCreateDoubleMatrix(1 , K , mxREAL);
	sum_alpha     = mxGetPr(plhs[3]);

	plhs[4]       = mxCreateDoubleMatrix(1 , K , mxREAL);
	sum_gamma     = mxGetPr(plhs[4]);
	
#else
	
#endif
	
	
	/*---------------------------------------------------------------*/
	/*--------- Internal Tempory vector & matrices ------------------*/
	/*---------------------------------------------------------------*/
	
	alpha         = (double *)malloc(NK*sizeof(double));
	At            = (double *)malloc(N*N*sizeof(double));
	temp_gamma    = (double *)malloc(N*sizeof(double));
	
#ifdef betanormalize
	
#else
	invsum_alpha  = (double *)malloc(K*sizeof(double));
	
#endif
	

	/*---------------------------------------------------------------*/
	/*------------------------ MAIN CALL ----------------------------*/
	/*---------------------------------------------------------------*/
	
#ifdef betanormalize
	forward_backward(PI , A  , L , X , gamma , logll , alpha , N , K , filter , At , temp_gamma , sum_alpha , sum_gamma);
#else
	forward_backward(PI , A  , L , X , gamma , logll , alpha , N , K , filter , At , temp_gamma , invsum_alpha);
#endif
	
	/*---------------------------------------------------------------*/
	/*------------------------ FREE MEMORY --------------------------*/
	/*---------------------------------------------------------------*/
	
    free(alpha);
	free(At);
	free(temp_gamma);
	
#ifdef betanormalize
#else	
	free(invsum_alpha);
#endif	
}

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

#ifdef betanormalize
void forward_backward(double *PI , double *A  , double *L , double *X  , double *gamma , double *logll , double *alpha   , int N , int K , int filter , double *At , double *temp_gamma , double *sum_alpha , double *sum_gamma )
#else
void forward_backward(double *PI , double *A  , double *L , double *X  , double *gamma , double *logll , double *alpha  , int N , int K , int filter , double *At , double *temp_gamma , double *invsum_alpha )
#endif
{
	int i , j , k , iN , kN1 , kN  , KN = (K - 1)*N;	
	double  cteN = 1.0/N, sumsum;
	
#ifdef betanormalize
	
#else
	double  sum_alpha=0.0 , sum_gamma = 0.0;
#endif
	
    register double sum , invsum;
	register int v;
	double maxi;
	
	/* At = A'; To speed up alpha computation */ 
	
	for(i = 0; i < N ; i++)	
	{
		iN = i*N;
		for(j = 0 ; j < N ; j++)
		{
			At[j + iN] = A[i + j*N];	
		}
	}
	
	logll[0]          = 0.0;
	
	/* Initialization of alpha */
	
#ifdef betanormalize
	sumsum            = 0.0;
	for (i=0 ; i<N ; i++)
	{
		alpha[i]      = PI[i]*L[i];
		sumsum       += alpha[i];
	}
    sum_alpha[0]     = sumsum;
	logll[0]         = log(sumsum);
	invsum           = 1.0/(sumsum + NUMERICS_FLOAT_MIN);
	for  (i=0 ; i<N ; i++)
	{ 
		alpha[i] *= invsum;
	}
	
	/* Loop of alpha : alpha(: , k +1) = L(: , k).*(A*alpha(: , k - 1)) */
	
	kN            = 0;
	for(k = 1 ; k < K ; k++)
	{
		
		/* 1a */
		
		kN           += N;
		kN1          = kN - N;
        sumsum       = 0.0;
		for(i = 0 ; i<N ; i++)
		{
			iN             = i*N;		
			sum            = 0.0;
			for(j = 0 ; j<N ; j++)
			{       
				sum   += At[j + iN]*alpha[j + kN1];
			}
			v              = i + kN;
			alpha[v]       = sum*L[v];
			sumsum        += alpha[v];
		}
		
		sum_alpha[k]       = sumsum;
		logll[0]          += log(sumsum);
		
		/* 1b normalization*/
		
		invsum           = 1.0/(sumsum + NUMERICS_FLOAT_MIN);
		for(i = 0 ; i<N ; i++)
		{
			alpha[i + kN] *= invsum;
		}
	}
	
	if (filter == 0)
	{	
		/* 2a */	
		sum_gamma[K-1] = cteN;
		for (i=0 ; i<N ; i++)
		{	
			gamma[i + KN] = cteN;
		}
		
        /* Loop of gamma : gamma(: , k) = At*(L(: , k+1).*gamma(: , k+1)) */
		
		for(k = K - 2 ; k>=0 ; k--)
		{
			/* 2b */
			
			kN         = (k+1)*N;
			for(i = 0 ; i < N ; i++)			
			{
				v              = i + kN;	
				temp_gamma[i]  = L[v]*gamma[v]; /* gamma(: , k+1).*L(: , k+1) */
			}
			sumsum           = 0.0;
			kN              -= N;
			for(i = 0 ; i<N ; i++)
			{
				sum       = 0.0;	
				iN        = i*N;
				for(j = 0 ; j<N ; j++)
				{       										
					sum   += A[j + iN]*temp_gamma[j];	 /* At*temp_gamma  */
				}
				gamma[i + kN] = sum;
				sumsum       += sum;
			}
			sum_gamma[k]     = sumsum;
			invsum           = 1.0/(sumsum + NUMERICS_FLOAT_MIN);
			for(i = 0 ; i<N ; i++)
			{
				gamma[i + kN] *= invsum;	
			}
		}
	}
#else	
	for (i=0 ; i<N ; i++)
	{
		alpha[i]     = PI[i]*L[i];	
		sum_alpha   += alpha[i];
	}
	
	logll[0]         = log(sum_alpha);
	invsum           = 1.0/(sum_alpha + NUMERICS_FLOAT_MIN);
	invsum_alpha[0]  = invsum;
	for  (i=0 ; i<N ; i++)
	{ 
		alpha[i] *= invsum;
	}

	/* Loop of alpha : alpha(: , k +1) = L(: , k).*(A*alpha(: , k - 1)) */
	
	kN            = 0;
	for(k = 1 ; k < K ; k++)
	{	
		/* 1a */	
		kN        += N;
		kN1       = kN - N;
		sum_alpha = 0.0;
		for(i = 0 ; i<N ; i++)
		{
			iN             = i*N;
			sum            = 0.0;
			for(j = 0 ; j<N ; j++)
			{       
				sum   += At[j + iN]*alpha[j + kN1];
			}
			
			v              = i + kN;
			alpha[v]       = sum*L[v];
			sum_alpha     += alpha[v];
		}
		logll[0]       += log(sum_alpha);
		
		/* 1b normalization*/
		invsum           = 1.0/(sum_alpha + NUMERICS_FLOAT_MIN);
		invsum_alpha[k]  = invsum;
		for(i = 0 ; i<N ; i++)
		{	
			alpha[i + kN] *= invsum;
		}
	}
	if (filter == 0)
	{		
		/* 2a */
			
		for (i=0 ; i<N ; i++)	
		{	
			gamma[i + KN] = cteN;
		}
		
        /* Loop of gamma : gamma(: , k) = At*(L(: , k+1).*gamma(: , k+1)) */
		
		for(k = K - 2 ; k>=0 ; k--)
		{	
			/* 2b */	
			kN         = (k+1)*N;
			for(i = 0 ; i<N ; i++)
			{
				v              = i + kN;	
				temp_gamma[i]  = L[v]*gamma[v]; /* gamma(: , k+1).*L(: , k+1) */
			}
			kN       -= N;
			for(i = 0 ; i<N ; i++)
			{
				sum       = 0.0;	
				iN        = i*N;
				for(j = 0 ; j<N ; j++)
				{       										
					sum   += A[j + iN]*temp_gamma[j];	 /* At*temp_gamma  */
				}
				gamma[i + kN] = sum*invsum_alpha[k];
			}
		}
	}
	
#endif
	
    /* X & Gamma outputs */ 	
   	
	if (filter == 0)
	{
		for (k=0 ; k<K ; k++)		
		{
			/* 3a */
			
			kN        = k*N;	
			sum       = 0.0;
			for (j = 0 ; j<N ; j++)
			{
				v              = j + kN;	
				gamma[v]      *= alpha[v];
				sum           += gamma[v];
			}
			invsum          = 1.0/(sum + NUMERICS_FLOAT_MIN);
			maxi            = 0.0;
			for (j=0 ; j<N ; j++)
			{
				v         = j + kN;
				gamma[v] *= invsum;
				if(gamma[v] > maxi)
				{
					X[k] = j;	
					maxi = gamma[v];
				}
			}
			X[k]++;
		}
	}
	if (filter == 1)
	{	
		for (i=0 ; i<K*N ; i++)		
		{
			gamma[i] = alpha[i];	
		}
		for (k=0 ; k<K ; k++)	
		{
			maxi = 0.0;	
			kN   = k*N;
			for (j=0 ; j<N ; j++)
			{
				v = j + kN;
				if(gamma[v] > maxi)
				{
					X[k] = j;	
					maxi = gamma[v];
				}
			}
			X[k]++;
		}
	}
}