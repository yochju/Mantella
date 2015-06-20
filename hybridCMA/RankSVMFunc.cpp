#include "math.h"

double DistPow2_Euclidean(double* x1, double* x2,int nx)
{
	double tmp;
	double dist = 0;
	for (int i=0; i<nx; i++)
	{
		tmp = x1[i] - x2[i];
		dist += tmp * tmp;
	}
	return dist;
}

void Encoding(double* x, double* x_Encoded, double* invsqrtC, double* xmean, int ntrain, int nx) //x'(i) = C^(-0.5) * ( x(i) - xmean(i) )
{
double sum;
double* dx = (double*) mxCalloc(nx, sizeof (double));
for (int i = 0; i < ntrain; i++) {
  for (int j = 0; j < nx; j++) {
    dx[j] = x[i* nx + j] - xmean[j];
  }
  for (int j = 0; j < nx; j++) {
    sum = 0;
    for (int k = 0; k < nx; k++) {
      sum += invsqrtC[j*nx + k] * dx[k];
    }
    x_Encoded[i * nx + j] = sum;
  }
}
mxFree(dx);
}


void EvaluatePointRankN(double* Fit, double* x_training, double* x_test, int ntest, int nx, int ntrain, double* p_alpha, double TwoSigmaPow2)
{
	double* x_Encoded = (double*)mxCalloc(nx,sizeof(double));
	double* Kvals = (double*)mxCalloc(ntrain,sizeof(double));

	for(int i=0; i<ntest; i++)
	{
		Encoding(&x_test[i*nx], x_Encoded, p_Cinv, xmean, 1, nx); 
		double* curX = x_Encoded;
	
		for (int j=0; j<ntrain; j++)
		{
          Kvals[j] = exp(- DistPow2_Euclidean(curX, &x_training[j*nx], nx) / TwoSigmaPow2);
		}

		double curFit = 0;
		for (int j=0; j<ntrain-1; j++)
			if (p_alpha[j] != 0)
				curFit += p_alpha[j] *( Kvals[j] - Kvals[j+1] );
		Fit[i] = curFit;
	}	
	mxFree(x_Encoded);
	mxFree(Kvals);
}


//y = RankSVMFunc(model.xtrainEncoded, x_test, npoints, model.N, ... 
//          	model.nTrain, model.alphas, model.TwoSigmaPow2, model.invsqrtC, model.Xmean_model, model.doEncoding, model.kernel);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//input
	double* X_training = mxGetPr(prhs[0]);
	double* X_test = mxGetPr(prhs[1]);
	int ntest = (int)(mxGetScalar(prhs[2]));
	int nx = (int)(mxGetScalar(prhs[3]));				
	int ntraining = (int)(mxGetScalar(prhs[4]));			
	double* p_alpha = mxGetPr(prhs[5]);
	double TwoSigmaPow2 = (double)(mxGetScalar(prhs[6]));	


	int rowLen = mxGetN(prhs[0]);						
	int colLen = mxGetM(prhs[0]);							
	if ((ntraining != rowLen)||(nx != colLen)) 
	{	 // MatLab uses [i*colLen + j] notation, while we use [i*rowLen + j], so .. :) 
		//printf("Error: the matrix 'x_training' should have 'nx' rows and 'ntraining' columns");
		return;
	}
	//output
	plhs[0] = mxCreateDoubleMatrix(1,ntest,  mxREAL);
	double* Fit = mxGetPr(plhs[0]);

	EvaluatePointRankN(Fit, X_training,X_test, ntest,nx, ntraining, p_alpha, TwoSigmaPow2);
}
