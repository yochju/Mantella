function model = xacmes_buildModel_RANKSVM(cur_state,coeff)

global algo;

N = cur_state.N;
aSZ = algo.aSZ;
if (0)
    ARX = algo.ARX;
    ARF = algo.ARF;
end;
if (1)
    ARneval = algo.ARneval(1,1:aSZ);
    [ARneval, arindex] = sort(ARneval,2,'ascend');
    ARX = algo.ARX(:,arindex);
    ARF = algo.ARF(1,arindex);
end;
invsqrtC = cur_state.invsqrtC;
Xmean_model = cur_state.xmean;

nTrainMax = aSZ;
npts = floor(coeff(1));
if (nTrainMax > npts)
    nTrainMax = npts;
end;

[xtrain, ftrain, nTrain] = xacmes_selectTrainingPoints(ARX, ARF, aSZ, N, nTrainMax);
if (nTrain < 0)
    disp( 'build model error: Ntrain < 0');
end;
            
nCrossValidation = 0;
CrossValidX = zeros(N,nCrossValidation);
CrossValidF = zeros(1,N);
for i=1:nCrossValidation
    index = 1 + floor( rand() * nTrain );
    nTrain = nTrain - 1;
    CrossValidX(:,i) = xtrain(index,:);
    CrossValidF(i) = ftrain(index);
    xtrain(index,:) = [];
    ftrain(index) = [];
end;
            
niter = floor( 1000*nTrain );
niter = floor( coeff(5)*nTrain );
epsilon = 1;
            
Cval = 10^6;    %default
Cval = 10^coeff(2);
nAlpha = nTrain - 1;
Ci = zeros(nAlpha,1);
z = 1:nAlpha;
           
Ci(z) = Cval*((nAlpha-z).^3.0 );    %default
Ci(z) = Cval*((nAlpha-z).^coeff(3) );
             
sigmaA = 1.0;   %default
sigmaA = coeff(4);
sigmaPow = 1.0;
x_tr = xtrain';   

%xmin = min(x_tr')';
%xmax = max(x_tr')';
%x_tr = (x_tr - repmat(xmin,1,nTrain))./ repmat(xmax-xmin,1,nTrain);
sigma = cur_state.sigma;
%x_tr = x_tr / sigma;
[xtrainEncoded, alphas, TwoSigmaPow2, xmin, xmax] = RankSVMLearn(x_tr, N, nTrain, niter, epsilon, Ci, ... 
                                                        invsqrtC, sigmaA, sigmaPow, Xmean_model);

 model.modelType = 1;
 model.sigma = sigma;
 model.kernelParam1 = TwoSigmaPow2;
 model.xmin = xmin;
 model.xmax = xmax;
 model.N = N;
 model.nTrain = nTrain;
 model.Xmean_model = Xmean_model;
 model.invsqrtC = invsqrtC;
 model.xtrainEncoded = xtrainEncoded;
 model.alphas = alphas;
 model.TwoSigmaPow2 = TwoSigmaPow2;
 
 model.nCrossValidation = nCrossValidation;
 model.CrossValidX = CrossValidX;
 model.CrossValidF = CrossValidF;
 

