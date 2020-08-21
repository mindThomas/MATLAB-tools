
/*

 em_ghmm : Expectation-Maximization algorithm for a HMM with Multivariate Gaussian measurements


 Usage  
 -------

 [logl , PI , A , M , S]  = em_ghmm(Z , PI0 , A0 , M0 , S0 , [options]);

 Inputs
 -------

 Z             Measurements (m x K x n1 x ... x nl)
 PI0           Initial proabilities (d x 1) : Pr(x_1 = i) , i=1,...,d. PI0 can be  (d x 1 x v1 x ... x vr)
 A0            Initial state transition probabilities matrix Pr(x_{k} = i| x_{k - 1} = j)  such
               sum_{x_k}(A0) = 1 => sum(A , 1) = 1. A0 can be (d x d x v1 x ... x vr).
 M0            Initial mean vector. M0 can be (m x 1 x d x v1 x ... x vr)
 S0            Initial covariance matrix. S0 can be (m x m x d x v1 x ... x vr)
 options       nb_ite         Number of iteration (default [30])
               update_PI      Update PI (0/1 = no/[yes])
               update_A       Update PI (0/1 = no/[yes])
               update_M       Update M  (0/1 = no/[yes])
               update_S       Update S  (0/1 = no/[yes])

 Ouputs
 -------
 logl          Final loglikelihood (n1 x ... x nl x v1 x ... x vr)
 PI            Estimated initial probabilities (d x 1 x n1 x ... x nl v1 x ... x vr)
 A             Estimated state transition probabilities matrix (d x d x n1 x ... x nl v1 x ... x vr)
 M             Estimated mean vector (m x 1 x d x n1 x ... x nl v1 x ... x vr)
 S             Estimated covariance vector (m x m x d x n1 x ... x nl v1 x ... x vr)


To compile
-----------

mex -output em_ghmm.dll em_ghmm.c


mex -f mexopts_intel10.bat em_ghmm.c


Example 1
----------

   d                                   = 2;
   m                                   = 2;
   L                                   = 1;
   R                                   = 1;
   Ntrain                              = 3000;
   Ntest                               = 10000;
   options.nb_ite                      = 30;

   PI                                  = [0.5 ; 0.5];
   A                                   = [0.95 0.05 ; 0.05 0.95];
   M                                   = cat(3 , [-1 ; -1] , [2 ; 2]);
   S                                   = cat(3 , [1 0.3 ; 0.3 0.8] , [0.7 0.6; 0.6 1]);

   [Ztrain , Xtrain]                   = sample_ghmm(Ntrain , PI , A , M , S , L);
   Xtrain                              = Xtrain - 1;

   %%%%% initial parameters %%%%

   PI0                                 = rand(d , 1 , R);
   sumPI                               = sum(PI0);
   PI0                                 = PI0./sumPI(ones(d , 1) , : , :);
   
   A0                                  = rand(d , d , R);
   sumA                                = sum(A0);
   A0                                  = A0./sumA(ones(d , 1) , : , :);

   M0                                  = randn(m , 1 , d , R);
   S0                                  = repmat(cat(3 , [2 0 ; 0 2] , [3 0; 0 2]) , [1 , 1 , 1, R]);

   %%%%% EM algorithm %%%%

   [logl , PIest , Aest , Mest , Sest] = em_ghmm(Ztrain , PI0 , A0 , M0 , S0 , options);


   [x , y]                             = ndellipse(M , S);
   [xest , yest]                       = ndellipse(Mest , Sest);

   Ltrain_est                          = likelihood_mvgm(Ztrain , Mest , Sest);
   Xtrain_est                          = forward_backward(PIest , Aest , Ltrain_est);
   Xtrain_est                          = Xtrain_est - 1;

   ind1                                = (Xtrain_est == 0);
   ind2                                = (Xtrain_est == 1);

   Err_train                           = min(sum(Xtrain ~= Xtrain_est , 2)/Ntrain , sum(Xtrain ~= ~Xtrain_est , 2)/Ntrain);

   figure(1) , 
   h                                   = plot(Ztrain(1 , ind1) , Ztrain(2 , ind1) , 'k+' , Ztrain(1 , ind2) , Ztrain(2 , ind2) , 'g+' , x , y , 'b' , xest  , yest ,'r', 'linewidth' , 2);
   legend([h(1) ; h(3:m:end)] , 'Train data' , 'True'  , 'Estimated' , 'location' , 'best')
   title(sprintf('Train data, Error rate = %4.2f%%' , Err_train*100))

   %%%%% Test data  %%%%


   [Ztest , Xtest]                     = sample_ghmm(Ntest , PI , A , M , S , L);
   Xtest                               = Xtest - 1;


   Ltest_est                           = likelihood_mvgm(Ztest , Mest , Sest);
   Xtest_est                           = forward_backward(PIest , Aest , Ltest_est);
   Xtest_est                           = Xtest_est - 1;


   ind1                                = (Xtest_est == 0);
   ind2                                = (Xtest_est == 1);

   Err_test                            = min(sum(Xtest ~= Xtest_est , 2)/Ntest , sum(Xtest ~= ~Xtest_est , 2)/Ntest);

   figure(2), 
   h                                   = plot(Ztest(1 , ind1) , Ztest(2 , ind1) , 'k+' , Ztest(1 , ind2) , Ztest(2 , ind2) , 'g+' , x , y , 'b' , xest  , yest ,'r', 'linewidth' , 2);
   legend([h(1) ; h(3:m:end)] , 'Test data' , 'True'  , 'Estimated' , 'location' , 'best')
   title(sprintf('Test data, Error rate = %4.2f%%' , Err_test*100))



   Author : Sébastien PARIS  © (sebastien.paris@lsis.org)
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


struct opts
{
  int    nb_ite ;
  int    update_PI;
  int    update_A;
  int    update_M;
  int    update_S;
};



/*--------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------------------*/

double ytRy(double *, double * , int , int);
void lubksb(double *, int , int *, double *);
int ludcmp(double *, int , int *, double * , double *);
double inv(double * , double * , double * , double * , int * , int );
void em_ghmm(double * , double * , double * , double * , double * , struct opts , 
			 double *  , double * , double * , double * , double * , 
			 int  , int  , int  , int  , int , 
			 double * , double * , double * , double * , double * , double * , double * , double * , double * , double * , double * , 
			 double * , double * , 
			 double * , double * , double *, double * , double * , double * , int * , double * , double *);

/*--------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{	
	mxArray *mxtemp;
	double *Z , *PI0 , *A0 , *M0 , *S0;
	double *logl , *PI , *A , *M , *S;
	double *gamma , *alpha , *beta , *dens , *PItemp , *Atemp , *Mtemp , *Stemp, *At ,  *temp_beta , *invsum_alpha;
	double *betadens , *alphabetadenst;
	double *tmp;
	int  *indx;
	double *vect , *vv , *invS  , *temp_S , *temp_invS , *det_S , *res , *resZ;
	const int *dimsZ , *dimsPI0 , *dimsA0 , *dimsM0 , *dimsS0;
	int *dimslogl , *dimsPI , *dimsA , *dimsM , *dimsS;
	int numdimsZ , numdimsPI0 , numdimsA0  , numdimsM0 , numdimsS0;
	int numdimslogl , numdimsPI  , numdimsA , numdimsM , numdimsS;
	int i , d2 , m2;
	int K , m  , d , L = 1 , R=1;
	struct opts options = {30 , 1 , 1 , 1 , 1};


	/*---------------------------------------------------------------*/
	/*---------------------- PARSE INPUT ----------------------------*/	
	/*---------------------------------------------------------------*/
	
	if( (nrhs < 5) )		
	{
		mexErrMsgTxt("At least 5 inputs are required");	
	}
	
	Z          = mxGetPr(prhs[0]);
    numdimsZ   = mxGetNumberOfDimensions(prhs[0]);
	dimsZ      = mxGetDimensions(prhs[0]);
	K          = dimsZ[1];
	for(i = 2 ; i < numdimsZ ; i++)
	{	
		L     *=dimsZ[i];	
	}
	
	PI0        = mxGetPr(prhs[1]);
    numdimsPI0 = mxGetNumberOfDimensions(prhs[1]);
	dimsPI0    = mxGetDimensions(prhs[1]);
	
	if (  dimsPI0[1] != 1 )
	{
		mexErrMsgTxt("PI must be (d x 1 x v1 x ... x vr ) ");			 	
	}	
	for(i = 2 ; i < numdimsPI0 ; i++)	
	{
		R     *=dimsPI0[i];	
	}
	
	A0         = mxGetPr(prhs[2]);
    numdimsA0  = mxGetNumberOfDimensions(prhs[2]);
	dimsA0     = mxGetDimensions(prhs[2]);
	
	if ( (dimsA0[1] != dimsA0[0]) )
	{	
		mexErrMsgTxt("A must be (d x d x v1 x ... x vr)");			 	
	}
	
	d          = dimsA0[0];
	d2         = d*d;
	
	M0         = mxGetPr(prhs[3]);
    numdimsM0  = mxGetNumberOfDimensions(prhs[3]);
	dimsM0     = mxGetDimensions(prhs[3]);
	
	if (  (dimsM0[2] != d)  )
	{	
		mexErrMsgTxt("M0 must be (m x 1 x  d x v1 x ... x vr)");			 	
	}
	
	m          = dimsM0[0];
	m2         = m*m;
	
	S0         = mxGetPr(prhs[4]);
    numdimsS0  = mxGetNumberOfDimensions(prhs[4]);
	dimsS0     = mxGetDimensions(prhs[4]);
	
	if (  (dimsS0[1] != m)  )
	{	
		mexErrMsgTxt("S0 must be (m x m x  d x v1 x ... x vr)");			 	
	}
	
	if(nrhs == 6)	
	{
		mxtemp            = mxGetField(prhs[5] , 0, "nb_ite");
		if(mxtemp != NULL)
		{
			tmp               = mxGetPr(mxtemp);	
			options.nb_ite    = (int) tmp[0];
		}
		
		mxtemp            = mxGetField(prhs[5] , 0, "update_PI");
		if(mxtemp != NULL)
		{
			tmp                = mxGetPr(mxtemp);	
			options.update_PI  = (int) tmp[0];
		}
		
		mxtemp            = mxGetField(prhs[5] , 0, "update_A");
		if(mxtemp != NULL)
		{
			tmp                = mxGetPr(mxtemp);	
			options.update_A   = (int) tmp[0];
		}
		
		mxtemp            = mxGetField(prhs[5] , 0, "update_M");
		if(mxtemp != NULL)
		{
			tmp                = mxGetPr(mxtemp);	
			options.update_M   = (int) tmp[0];
		}

		mxtemp            = mxGetField(prhs[5] , 0, "update_S");
		if(mxtemp != NULL)
		{	
			tmp                = mxGetPr(mxtemp);	
			options.update_S   = (int) tmp[0];
		}
	}
	
	/*---------------------------------------------------------------*/
	/*---------------------- PARSE OUTPUT ---------------------------*/	
	/*---------------------------------------------------------------*/
	
	numdimslogl    = max(numdimsZ - 2 + numdimsPI0 - 2 , 2);
	dimslogl       = (int *)malloc(numdimslogl*sizeof(int));
	dimslogl[0]    = 1;
	dimslogl[1]    = 1;
	
	numdimsPI      = 2 + max(numdimsZ - 2 , 0) + max(numdimsPI0 - 2 , 0);
	dimsPI         = (int *)malloc(numdimsPI*sizeof(int));
	dimsPI[0]      = d;
	dimsPI[1]      = 1;
	
	
	numdimsA       = 2 + max(numdimsZ - 2 , 0) + max(numdimsPI0 - 2 , 0);
	dimsA          = (int *)malloc(numdimsA*sizeof(int));
	dimsA[0]       = d;
	dimsA[1]       = d;
	
	
	numdimsM       = 3 + max(numdimsZ - 2 , 0) + max(numdimsPI0 - 2 , 0);
	dimsM          = (int *)malloc(numdimsM*sizeof(int));
	dimsM[0]       = m;
	dimsM[1]       = 1;
	dimsM[2]       = d;
	
	numdimsS       = 3 + max(numdimsZ - 2 , 0) + max(numdimsPI0 - 2 , 0);
	dimsS          = (int *)malloc(numdimsS*sizeof(int));
	dimsS[0]       = m;
	dimsS[1]       = m;
	dimsS[2]       = d;
	
	for(i = 2 ; i < numdimsZ ; i++)	
	{
		dimslogl[i - 2] = dimsZ[i];	
		dimsPI[i]       = dimsZ[i];
		dimsA[i]        = dimsZ[i];
		dimsM[i+1]      = dimsZ[i];
		dimsS[i+1]      = dimsZ[i];
	}
	for(i = 2 ; i < numdimsPI0 ; i++)	
	{
		dimslogl[i  - 2 + (numdimsZ - 2)] = dimsPI0[i];	
		dimsPI[i + (numdimsZ - 2)]        = dimsPI0[i];
		dimsA[i + (numdimsZ - 2)]         = dimsPI0[i];
		dimsM[i + 1 + (numdimsZ - 2)]     = dimsPI0[i];
		dimsS[i + 1 + (numdimsZ - 2)]     = dimsPI0[i];
	}
	
	
	plhs[0]        = mxCreateNumericArray(numdimslogl , dimslogl , mxDOUBLE_CLASS, mxREAL);
	logl           = mxGetPr(plhs[0]);
	
	plhs[1]        = mxCreateNumericArray(numdimsPI , dimsPI , mxDOUBLE_CLASS, mxREAL);
	PI             = mxGetPr(plhs[1]);
	
	plhs[2]        = mxCreateNumericArray(numdimsA , dimsA , mxDOUBLE_CLASS, mxREAL);
	A              = mxGetPr(plhs[2]);
	
	plhs[3]        = mxCreateNumericArray(numdimsM , dimsM , mxDOUBLE_CLASS, mxREAL);
	M              = mxGetPr(plhs[3]);
	
	plhs[4]        = mxCreateNumericArray(numdimsS , dimsS , mxDOUBLE_CLASS, mxREAL);
	S              = mxGetPr(plhs[4]);
	
	
	gamma          = (double *)malloc(d*K*sizeof(double));	
	alpha          = (double *)malloc(d*K*sizeof(double));
	beta           = (double *)malloc(d*K*sizeof(double));
	dens           = (double *)malloc(d*K*sizeof(double));
	PItemp         = (double *)malloc(d*sizeof(double));
	
	Atemp          = (double *)malloc(d2*sizeof(double));
	At             = (double *)malloc(d2*sizeof(double));
	Mtemp          = (double *)malloc(m*d*sizeof(double));
	Stemp          = (double *)malloc(m2*d*sizeof(double));
	temp_beta      = (double *)malloc(d*sizeof(double));
	invsum_alpha   = (double *)malloc(K*sizeof(double));
	betadens       = (double *)malloc(d*(K-1)*sizeof(double));
	alphabetadenst = (double *)malloc(d2*sizeof(double));
	
	
	vect           = (double *)malloc(m*sizeof(double));
	vv             = (double *)malloc(m*sizeof(double));
	temp_S         = (double *)malloc(m2*sizeof(double));
	temp_invS      = (double *)malloc(m2*sizeof(double));
	det_S          = (double *)malloc(d*sizeof(double));
	invS           = (double *)malloc((m2*d)*sizeof(double));
	indx           = (int *)malloc(m*sizeof(int));
	res            = (double *)malloc(m*sizeof(double));
	resZ           = (double *)malloc(m*d*K*sizeof(double));
	
	
	/*---------------------------------------------------------------*/
	/*------------------------ MAIN CALL ----------------------------*/
	/*---------------------------------------------------------------*/
	
	
	em_ghmm(Z , PI0 , A0  , M0  , S0 , options , logl , PI , A , M , S , m , d , K , L , R , gamma , alpha , beta , dens , PItemp , Atemp , Mtemp , Stemp , At , temp_beta , invsum_alpha , betadens , alphabetadenst ,  invS , temp_S , temp_invS , det_S , vect , vv , indx , res , resZ);
	
	
	/*---------------------------------------------------------------*/
	/*------------------------ FREE MEMORY --------------------------*/
	/*---------------------------------------------------------------*/
	
    free(dimslogl);
    free(dimsPI);
    free(dimsA);
	free(dimsM);
	free(dimsS);
	free(dens);
	free(gamma);
	free(alpha);
	free(beta);
	free(PItemp);
	free(Atemp);
	free(At);
	free(Mtemp);
	free(Stemp);
	free(temp_beta);
	free(invsum_alpha);
	free(betadens);
	free(alphabetadenst);
	
	free(vect);
	free(vv);
	free(temp_S);
	free(temp_invS);
	free(det_S);
	free(invS);
	free(indx);
	free(res);
	free(resZ);
}

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

void em_ghmm(double *Z , double *PI0 , double *A0  , double *M0 , double *S0 , struct opts options , double *logl  , double *PI , double *A , double *M , double *S , int m , int d , int K , int L , int R , double *gamma , double *alpha , double *beta , double *dens , double *PItemp , double *Atemp , double *Mtemp , double *Stemp , double *At , double *temp_beta , double *invsum_alpha , double *betadens , double *alphabetadenst , double *invS , double *temp_S , double *temp_invS , double *det_S , double *vect , double *vv , int *indx , double *res , double *resZ)			 
{
	int i , j , l , k  , h , t , r , d2 = d*d , m2 = m*m, rd , rd2 , md = m*d , mK = m*K , rmd , m2d = m2*d , rm2d, id , rL , dK = K*d;
	int kd , kd1 , lmK , index , index1 , Kd = (K-1)*d , rLd , rLd2 , rLmd , Ld = L*d , Lmd = Ld*m;
	int ld , ld2 , lmd , rLm2d , lm2d , Lm2d = L*m2d , im2 , km , jm , jmK , imK ;
	double  sum_gamma , invsum_gamma  , cted = 1.0/d ;
	double cte = 1.0/pow(2*M_PI , m/2.0);
	double  sum_alpha;
    register double sum , invsum;
	register int v;
	
	for(r = 0 ; r < R ; r++)	
	{
        rd    = r*d;	
        rd2   = r*d2;
		rmd   = r*md;
        rm2d  = r*m2d;
		rL    = r*L;	
		rLd   = r*Ld;
		rLd2  = rLd*d;
		rLmd  = r*Lmd;
		rLm2d = r*Lm2d;
		
		for (l = 0 ; l < L ; l++)	
		{
			lmK  = l*mK;	
			ld   = l*d  + rLd;
			ld2  = l*d2 + rLd2;
			lmd  = l*md + rLmd;
			lm2d = l*m2d + rLm2d;	
			for(i = 0 ; i < d ; i++)	
			{
				PItemp [i] = PI0[i + rd]; 	
			}
			for(i = 0 ; i < d2 ; i++)
			{
				Atemp [i] = A0[i + rd2]; 	
			}
			for(i = 0 ; i < md ; i++)	
			{
				Mtemp [i] = M0[i + rmd]; 	
			}
			
			for(i = 0 ; i < m2d ; i++)	
			{
				Stemp [i] = S0[i + rm2d]; 	
			}
			
			for(t = 0 ; t < options.nb_ite ; t++)	
			{
				
				/* At = Atemp' */

				for(i = 0 ; i < d ; i++)	
				{
					id = i*d;
					for(j = 0 ; j < d ; j++)
					{
						At[j + id]     = Atemp[i + j*d];	
					}
				}
				
				/* invS = inv(S); */
					
					for (i = 0 ; i < d ; i++)						
					{
						im2   = i*m2;				
						for(j = 0 ; j < m2 ; j++)		
						{
							temp_S[j] = Stemp[j + im2];	
						}
						det_S[i]  = inv(temp_S , temp_invS , vect , vv , indx , m);
						for(j = 0 ; j < m2 ; j++)	
						{
							invS[j + im2] = temp_invS[j];	
						}
						det_S[i] = (cte*sqrt(fabs(det_S[i])));
					}
	
				
				/* dens = cte*exp(-0.5res'*invS*res) */
				
				for (k = 0 ; k < K ; k++)			
				{
					km    = k*m + lmK;	
					kd    = k*d;
					for (j = 0 ; j < d ; j++)	
					{
						jm            = j*m;	
						for(i = 0 ; i < m ; i++)
						{
							res[i] = (Z[i + km] - Mtemp[i + jm]);	
						}
						dens[j + kd]  = det_S[j]*exp(-0.5*ytRy(res , invS , m , j*m2)) + NUMERICS_FLOAT_MIN;
					}	
				}			 
						
			    /* Alpha probabilies */				
				
				logl[l + rL] = 0.0;
				sum_alpha    = 0.0;
				
				for (i=0 ; i<d ; i++)
				{
					alpha[i]   = PItemp[i]*dens[i];	
					sum_alpha += alpha[i];
				}
				
				logl[l + rL]     = log(sum_alpha);
				invsum           = 1.0/(sum_alpha + NUMERICS_FLOAT_MIN);
				invsum_alpha[0]  = invsum;
				
				for  (i=0 ; i<d ; i++)
				{ 
					alpha[i] *= invsum;
				}
				
				kd            = 0;
				for(k = 1 ; k < K ; k++)		
				{
					kd        += d;	
					kd1        = kd - d;
					sum_alpha  = 0.0;
					
					for(i = 0 ; i<d ; i++)
					{
						id             = i*d;			
						sum            = 0.0;
						
						for(j = 0 ; j < d ; j++)
						{       
							sum   += Atemp[j + id]*alpha[j + kd1];
						}
						
						v              = i + kd;
						alpha[v]       = sum*dens[v];
						sum_alpha     += alpha[v];	
					}
					logl[l + rL]     += log(sum_alpha);
					invsum           = 1.0/(sum_alpha + NUMERICS_FLOAT_MIN);
					invsum_alpha[k]  = invsum;
					
					for(i = 0 ; i < d ; i++)
					{	
						alpha[i + kd] *= invsum;
					}
				}
											
				for (i=0 ; i<d ; i++)	
				{	
					beta[i + Kd] = cted;
				}
				for(k = K - 2 ; k>=0 ; k--)	
				{
					kd                 = (k + 1)*d;	
					for(i = 0 ; i < d ; i++)
					{
                        v              = i + kd;			
						temp_beta[i]   = beta[v]*dens[v]; 
					}
					
					kd       -= d;
					for(i = 0 ; i<d ; i++)			
					{
						sum       = 0.0;	
						id        = i*d;
						for(j = 0 ; j<d ; j++)
						{       										
							sum   += At[j + id]*temp_beta[j];	 
						}
						beta[i + kd] = sum*invsum_alpha[k];
					}
				}				
				for(j = 0 ; j < d ; j++)
				{
					temp_beta[j] = 0.0;	
				}
				for (k=0 ; k < K ; k++)	
				{
					kd        = k*d;	
					sum_gamma = 0.0;
					for (j = 0 ; j<d ; j++)
					{
						v              = j + kd;	
						gamma[v]       = alpha[v]*beta[v];
						sum_gamma     += gamma[v];
					}
					
					invsum_gamma    = 1.0/(sum_gamma + NUMERICS_FLOAT_MIN);		
					for (j = 0 ; j < d ; j++)
					{
						gamma[j + kd] *= invsum_gamma;	
						temp_beta[j] += gamma[j + kd];
					}
				}
				
				for(i = 0 ; i < d ; i++)				
				{
					temp_beta[i] = 1.0/(temp_beta[i] + NUMERICS_FLOAT_MIN);	
				}
				
				if(options.update_PI)
				{
					for(i = 0 ; i < d ; i++)		
					{
						PItemp[i] = gamma[i];
						
					}
				}
							
				if(options.update_A)
				{	
					for (i = d ; i < dK ; i++)		
					{
						betadens[i - d] = beta[i]*dens[i];	
					}
					
					for(i = 0 ; i < d ; i++)
					{
						id    = i*d;
						for(j = 0 ; j < d ; j++)
						{
							v                   = j + id;	
							alphabetadenst[v]   = 0.0;						
							for (k = 0 ; k < K - 1 ; k++)
							{
								kd                 = k*d;	
								alphabetadenst[v] += alpha[j + kd]*betadens[i + kd];
							}
						}          
					}
					
					for(i = 0 ; i < d ; i++)	
					{
						id  = i*d;	
						sum = 0.0;
						
						for(j = 0 ; j < d ; j++)
						{
							v         = j + id;	
							Atemp[v] *= alphabetadenst[v];
							sum      += Atemp[v];
						}
						
						invsum = 1.0/(sum + NUMERICS_FLOAT_MIN);
						for(j = 0 ; j < d ; j++)
						{
							Atemp[j + id] *= invsum;	
						}				
					}
				}
				
				if(options.update_M)
				{
					for (j = 0 ; j < d ; j++)	
					{
						jm            = j*m;	
						for(i = 0 ; i < m ; i++)
						{
							sum       = 0.0;	
							for (k = 0 ; k < K ; k++)
							{
								sum += Z[i + k*m + lmK]*gamma[j + k*d];	
							}
							Mtemp[i + jm] = sum*temp_beta[j];
						}
					}
				}
					
				if(options.update_S)
				{
					/* resZ(m x K x d) */
					
					for(j = 0 ; j < d ; j++)
					{	
						jmK = j*mK;	
						jm  = j*m;
						for(k = 0 ; k < K ; k++)
						{
							km     = k*m;	
							index  = km + lmK;
							index1 = km + jmK;
							for(i = 0 ; i < m ; i++)
							{	
								resZ[i + index1] = (Z[i + index] - Mtemp[i + jm]);	
							}
						}
					}
					
					/* resZ(m x K x d) */
					
					for(i = 0 ; i < d ; i++)					
					{
						im2 = i*m2;
						imK = i*mK;
						for (j = 0 ; j < m ; j ++)	
						{	
							for (h = 0 ; h <= j ; h++)		
							{
								sum = 0.0;
								for (k = 0 ; k < K ; k++)
								{
									km    = k*m + imK;
									kd    = k*d;				
									sum  += resZ[j + km]*resZ[h + km]*gamma[i + kd];
								}
								Stemp[j + h*m + im2] = sum*temp_beta[i];								
							}		
						}		
						for (j = 0 ; j < m - 1 ; j++)	
						{
							jm     = j*m;	
							for (h = j + 1 ; h < m ; h++)
							{
								Stemp[h*m + j + im2] = Stemp[h + jm + im2];				
							}				
						}
					}	
				}
			}
			
			for(i = 0 ; i < d ; i++)	
			{
				PI[i + ld] = PItemp[i]; 
			}
			for(i = 0 ; i < d2 ; i++)
			{
				A[i + ld2] = Atemp[i]; 	
			}			
			for(i = 0 ; i < md ; i++)	
			{
				M[i + lmd] = Mtemp[i]; 	
			}
			for(i = 0 ; i < m2d ; i++)	
			{
				S[i + lm2d] = Stemp[i]; 	
			}
		}
	}
}

/*----------------------------------------------------------------------------------------------*/
double ytRy(double *y, double *R , int d , int offset)
{
	int  i , j , id;
	register double temp;
	register double Q = 0.0;
	
	for (i = 0 ; i < d ; i++)
	{	
		temp = 0.0;	
		id   = i*d + offset;
		for(j = 0 ; j < d ; j++)
		{
			temp   += y[j]*R[j + id];	
		}
		Q += temp*y[i];
	}
	return Q;
}
/*------------------------------------------------------------------*/
double inv(double *temp , double *invQ  , double *vect , double *vv , int *indx , int d)
{
	int i , j , jd;	
	double dd , det = 1.0;

	if(ludcmp(temp , d , indx , &dd , vv ))
	{	
		for(i = 0 ; i < d ; i++)		
		{			
			det *= temp[i + i*d];		
		}		
		for(j = 0; j < d; j++)
		{            
			for(i = 0; i < d; i++) 				
			{
				vect[i] = 0.0;
			}			
			jd      = j*d;			
			vect[j] = 1.0;			
			lubksb(temp , d , indx , vect);			
			for(i = 0 ; i < d ; i++) 				
			{				
				invQ[jd + i] = vect[i];				
			}
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
        for(j = 0; j < n; j++)
		{
            if((temp = fabs(m[i + j*n])) > big)			
			{
                big = temp;
            }		
		}
        if(big == 0.0)
		{			
            return 0;
        }		
        vv[i] = 1.0 / big;
    }	
    for(j = 0; j < n; j++)
	{
		jn  = j*n;		
        for(i = 0; i < j; i++)			
		{
            sum = m[i + jn];			
            for(k = 0 ; k < i; k++)				
			{
                sum -= m[i + k*n ] * m[k + jn];
            }           
			m[i + jn] = sum;
        }		
        big = 0.0;		
        for(i = j; i < n; i++)			
		{
            sum = m[i + jn];			
            for(k = 0; k < j; k++)				
			{				
				sum -= m[i + k*n] * m[k + jn];				
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
            for(k = 0; k < n; k++)				
			{				
				kn            = k*n;				
                dum           = m[imax + kn];
                m[imax + kn]  = m[j + kn];
                m[j + kn]     = dum;				
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
    }
    return 1;
};

/*-------------------------------------------------------------------------*/
