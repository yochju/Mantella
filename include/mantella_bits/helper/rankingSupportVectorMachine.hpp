namespace mant {

  class RankingSVM {
    public:
      explicit RankingSVM(arma::Mat<double> xTraining, unsigned int niter, unsigned int ntraining, arma::Col<double> Ci, arma::Mat<double> invsqrtC, arma::Col<double> xMean) noexcept;
      void LearningRankN();
      
      void SetEpsilon(double epsilon);
      double GetEpsilon() const;
      void SetSigmaPow(double sigmaPow);
      double GetSigmaPow() const;
      void SetSigmaA(double sigmaA);
      double GetSigmaA() const;
      
    protected:
      void Encoding();
      void CalculateTrainingKernelMatrix();
      void OptimizeL();
      
      //RankSVMLearn.cpp Input Variables
      double sigmaA = 1;
      double sigmaPow = 1;
      double epsilon = 1;
      unsigned int dimension; //N
      arma::Mat<double> xTraining;
      unsigned int niter;
      unsigned int ntraining;
      arma::Col<double> Ci;
      arma::Mat<double> invsqrtC;
      arma::Col<double> xMean;
      
      //RankSVMLearn.cpp Output Variables
      arma::Mat<double> xTrainingEncoded;
      arma::Col<double> optAlphas;
      double TwoSigmaPow2 = 0;
      
      //RankSVMLearn.cpp Variables Between Functions
      unsigned int nAlpha; //may need to be int, not sure if ntraining can be 0
      arma::Mat<double> K;
  };

  //
  // Implementation
  //

  RankingSVM::RankingSVM(arma::Mat<double> xTraining, unsigned int niter, unsigned int ntraining, arma::Col<double> Ci, arma::Mat<double> invsqrtC, arma::Col<double> xMean) noexcept {
    this->xTraining = xTraining;
    this->niter = niter;
    this->ntraining = ntraining;
    this->Ci = Ci;
    this->invsqrtC = invsqrtC;
    this->xMean = xMean;
    
    this->dimension = xTraining.n_rows;
    assert(ntraining==xTraining.n_cols);
    
    this->xTrainingEncoded = arma::zeros(this->dimension,ntraining);
    this->optAlphas = arma::zeros(ntraining-1);
    
    K = arma::zeros(ntraining,ntraining);
    nAlpha = ntraining-1;    
  }
  
  void RankingSVM::LearningRankN() {
    Encoding();
    CalculateTrainingKernelMatrix();
    OptimizeL();
  }
  
  void RankingSVM::Encoding() {
    for(int i = 0; i < ntraining; i++) {
      arma::Col<double> dx = xTraining.col(i) - xMean;
      xTrainingEncoded.col(i) = arma::sum(invsqrtC*dx);
    }
  }
  
  void RankingSVM::CalculateTrainingKernelMatrix() {
    double avrdist = 0;
    for(int i=0; i < ntraining; i++) {
      K(i,i) = 0;
      for(int j=i+1; j < ntraining; j++) {
        K(j,i) = std::pow(arma::norm(xTrainingEncoded.col(i) - xTrainingEncoded.col(j)),2);
        K(i,j) = K(j,i);
        avrdist += std::sqrt(K(i,j));
      }
    }
    
    avrdist /= (ntraining-1)*ntraining / 2.0;
    TwoSigmaPow2 = 2 * std::pow((sigmaA * avrdist * sigmaPow),2);
    
    //the for loops can be combined and the upper assignments to K never made.
    //which would be much more efficient. If arma can do Mat minus Mat,
    //the upper for-loop is completely unnecessary
    for(int i=0; i < ntraining; i++) {
      K(i,i) = 1;
      for(int j = i+1; j < ntraining; j++) {
        K(j,i) = exp(-K(j,i) / TwoSigmaPow2);
        K(i,j) = K(j,i);
      }
    }
  }
  
  void RankingSVM::OptimizeL() {
    arma::Col<double> sumAlphaDK = arma::zeros(nAlpha);
    arma::Mat<double> divDK = arma::zeros(nAlpha,nAlpha);
    arma::Mat<double> dK = arma::zeros(nAlpha,nAlpha);
    
    //there is probably a way to combine all these for loops
    for(int i = 0; i < nAlpha; i++) {
      for(int j = 0; j < nAlpha; j++) {
        dK(j,i) = K(j,i) - K(j+1,i) - K(j, i+1) + K(j+1,i+1);
      }
      optAlphas(i) = Ci(i) * (0.95 + 0.05 * arma::randu());
    }
    
    for(int i = 0; i < nAlpha; i++) {
      double sumAlpha = arma::sum(optAlphas % dK.col(i));
      sumAlphaDK(i) = (epsilon - sumAlpha) / dK(i,i);
    }
    
    for(int i = 0; i < nAlpha; i++) {
      for(int j = 0; j < nAlpha; j++) {
        divDK(j,i) = dK(j,i) / dK(j,j);
      }
    }
    
    for(int i = 0; i < niter; i++) {
      int iMod = i%nAlpha;
      double newAlpha = optAlphas(iMod) + sumAlphaDK(iMod);
      if(newAlpha > Ci(iMod)) {
        newAlpha = Ci(iMod);
      }
      if(newAlpha < 0) {
        newAlpha = 0;
      }
      double deltaAlpha = newAlpha - optAlphas(iMod);
      
      double dL = deltaAlpha * dK(iMod,iMod) * (sumAlphaDK(iMod)- 0.5 * deltaAlpha);
      
      if(dL > 0) {
        sumAlphaDK -= deltaAlpha * divDK.col(iMod);
        optAlphas(iMod) = newAlpha;
      }
    }
  }

  void RankingSVM::SetEpsilon(double epsilon) {
    this->epsilon = epsilon;
  }

  double RankingSVM::GetEpsilon() const {
    return epsilon;
  }

  void RankingSVM::SetSigmaPow(double sigmaPow) {
    this->sigmaPow = sigmaPow;
  }

  double RankingSVM::GetSigmaPow() const {
    return sigmaPow;
  }

  void RankingSVM::SetSigmaA(double sigmaA) {
    this->sigmaA = sigmaA;
  }

  double RankingSVM::GetSigmaA() const {
    return sigmaA;
  }
}