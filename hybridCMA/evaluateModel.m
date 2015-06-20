function Fit = evaluateModel(x_test, npoints, model)

if (model.modelType == 1)
%    x_test = (x_test - repmat(model.xmin,1,npoints))./ repmat(model.xmax-model.xmin,1,npoints);
 %   x_test = x_test / model.sigma;
    Fit = -RankSVMFunc(model.xtrainEncoded, x_test, npoints, model.N, ... 
            model.nTrain, model.alphas, model.TwoSigmaPow2); 

end;
