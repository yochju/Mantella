#include "mex.h" 
#include "math.h"
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

void CalculateTrainingKernelMatrix(double* x_Encoded, int ntrain, int nx, double* K,
    double* TwoSigmaPow2, double sigma_A, double sigma_Pow) {
  double distPow2, tmp;
  double avrdist = 0;
  for (int i = 0; i < ntrain; i++) {
    K[i * ntrain + i] = 0;

    for (int j = i + 1; j < ntrain; j++) {
      distPow2 = 0;
      for (int k = 0; k < nx; k++) {
        tmp = x_Encoded[i*nx + k] - x_Encoded[j*nx + k];
        distPow2 += tmp * tmp;
      }
      K[i * ntrain + j] = K[j * ntrain + i] = distPow2;
      avrdist += sqrt(distPow2);
    }
  }

  avrdist /= ((ntrain - 1) * ntrain / 2);
  double sigma = sigma_A * pow(avrdist, sigma_Pow);
  *TwoSigmaPow2 = 2.0 * sigma*sigma;

  for (int i = 0; i < ntrain; i++) {
    K[i * ntrain + i] = 1;
    for (int j = i+1; j < ntrain; j++) {
      K[i * ntrain + j] = K[j * ntrain + i] = exp(-K[i * ntrain + j] / *TwoSigmaPow2);
    }
  }
}

void OptimizeL(int ntrain, double* p_Ci, double epsilon, int niter,
    double* p_Kij, double* p_dKij, double* p_alpha, double* p_sumAlphaDKij, double* p_div_dKij) {
  int nAlpha = ntrain - 1;
  double old_alpha, new_alpha, delta_alpha, sumAlpha, dL;
  int i, i1, j;

  double* p_dKii = new double[ntrain];

  for (i = 0; i < nAlpha; i++) {
    for (j = 0; j < nAlpha; j++) {
      p_dKij[i * nAlpha + j] = p_Kij[i * ntrain + j] - p_Kij[i * ntrain + (j + 1)] - p_Kij[(i + 1) * ntrain + j] + p_Kij[(i + 1) * ntrain + (j + 1)];
    }
    p_dKii[i] = p_dKij[i * nAlpha + i];
  }

  for (i = 0; i < nAlpha; i++) {
    p_alpha[i] = p_Ci[i] * (0.95 + 0.05 * rand() / (float) RAND_MAX);
  }

  for (i = 0; i < nAlpha; i++) {
    sumAlpha = 0;
    for (j = 0; j < nAlpha; j++) {
      sumAlpha += p_alpha[j] * p_dKij[i * nAlpha + j];
    }
    p_sumAlphaDKij[i] = (epsilon - sumAlpha) / p_dKij[i * nAlpha + i];
  }

  for (i = 0; i < nAlpha; i++)
    for (j = 0; j < nAlpha; j++)
      p_div_dKij[i * nAlpha + j] = p_dKij[i * nAlpha + j] / p_dKij[j * nAlpha + j];

  for (i = 0; i < niter; i++) {
    i1 = i % nAlpha;
    old_alpha = p_alpha[i1];
    new_alpha = old_alpha + p_sumAlphaDKij[i1];
    if (new_alpha > p_Ci[i1]) {
      new_alpha = p_Ci[i1];
    }
    if (new_alpha < 0) {
      new_alpha = 0;
    }
    delta_alpha = new_alpha - old_alpha;

    dL = delta_alpha * p_dKii[i1] * (p_sumAlphaDKij[i1] - 0.5 * delta_alpha + 0);

    if (dL > 0) {
      for (j = 0; j < nAlpha; j++) {
        p_sumAlphaDKij[j] -= delta_alpha * p_div_dKij[i1 * nAlpha + j];
      }

      p_alpha[i1] = new_alpha;
    }
  }
}

void LearningRankN(double* x_training, double* x_trainingEncoded,
    int nx, int ntrain, int niter, double epsilon, double* p_Ci, double* p_Cinv,
    double sigma_A, double sigma_Pow, double* xmean, double* optAlphas, double* pTwoSigmaPow2) {
  //0. Init Temp Data
  int nAlpha = ntrain - 1;
  double* p_Kij = (double*) mxCalloc(ntrain*ntrain, sizeof (double));
  double* p_dKij = (double*) mxCalloc(nAlpha*nAlpha, sizeof (double));
  double* p_alpha = (double*) mxCalloc(nAlpha, sizeof (double));
  double* p_sumAlphaDKij = (double*) mxCalloc(nAlpha, sizeof (double));
  double* p_div_dKij = (double*) mxCalloc(nAlpha*nAlpha, sizeof (double));

  //1. Transform training points to the new coordinate system and 
  //then compute Euclidean distance instead of 'EXPENSIVE' Mahalanobis distance

  //p_Cinv is the C^-0.5 ;)
  Encoding(x_training, x_trainingEncoded, p_Cinv, xmean, ntrain, nx);

  //2. Calculate the distance between points, then calculate sigma(gamma) and Kernel Matrix
  double TwoSigmaPow2 = 0;
  CalculateTrainingKernelMatrix(x_trainingEncoded, ntrain, nx, p_Kij, &TwoSigmaPow2, sigma_A, sigma_Pow);


  //3. Optimize alpha parameters
  OptimizeL(ntrain, p_Ci, epsilon, niter, p_Kij, p_dKij, p_alpha, p_sumAlphaDKij, p_div_dKij);


  *pTwoSigmaPow2 = TwoSigmaPow2;
  for (int i = 0; i < nAlpha; i++) {
    optAlphas[i] = p_alpha[i];
  }

  mxFree(p_sumAlphaDKij);
  mxFree(p_div_dKij);
  mxFree(p_dKij);
  mxFree(p_Kij);
  mxFree(p_alpha);
}

//[xtrainEncoded, alphas, TwoSigmaPow2] = RankSVMLearn(x_tr, nx, ntrain, niter, epsilon, Ci, kernel, ... 
//                                                        invsqrtC, sigmaA, sigmaPow, xmean, doEncoding, verbose);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  //input
  double* X_training = mxGetPr(prhs[0]);
  int nx = (int) (mxGetScalar(prhs[1]));
  int ntraining = (int) (mxGetScalar(prhs[2]))
  int niter = (int) (mxGetScalar(prhs[3]));
  double epsilon = (double) (mxGetScalar(prhs[4]));
  double* p_Ci = mxGetPr(prhs[5]);
  double* p_Cinv = mxGetPr(prhs[7]);
  double sigma_A = (double) (mxGetScalar(prhs[8]));
  double sigma_Pow = (double) (mxGetScalar(prhs[9]));
  double* xmean = mxGetPr(prhs[10]);
  //
  int rowLen = mxGetN(prhs[0]);
  int colLen = mxGetM(prhs[0]);
  if ((ntraining != rowLen) || (nx != colLen)) { // MatLab uses [i*colLen + j] notation, while we use [i*rowLen + j], so .. :) 
    //printf("Error: the matrix 'x_training' should have 'nx' rows and 'ntraining' columns");
    return;
  }

  //output
  plhs[0] = mxCreateDoubleMatrix(nx, ntraining, mxREAL);
  double* X_trainingEncoded = mxGetPr(plhs[0]);

  plhs[1] = mxCreateDoubleMatrix(ntraining - 1, 1, mxREAL);
  double* optAlphas = mxGetPr(plhs[1]);

  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double* pTwoSigmaPow2 = mxGetPr(plhs[2]);

  LearningRankN(X_training, X_trainingEncoded, nx, ntraining, niter, epsilon, p_Ci,
      p_Cinv, sigma_A, sigma_Pow, xmean, optAlphas, pTwoSigmaPow2);
}
