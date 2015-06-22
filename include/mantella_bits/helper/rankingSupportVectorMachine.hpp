namespace mant {

  template <typename T>
  class RankingSupportVectorMachine {
  public:
    explicit RankingSupportVectorMachine() noexcept;
    void learn(const arma::Mat<T>& parameters);
    
    void evaluatePoints(const arma::Mat<T>& points, const arma::Mat<T>& parameters, const arma::Col<T>& rankingDirection, T twoSigmaPow2);

    arma::Col<T>& getRankingDirection();
    
    arma::Col<T>& getFitnessForEvaluatedPoints();
    
    T getTwoSigmaPow2();

    void setMaximalNumberOfIterations(const unsigned long long maximalNumberOfIterations);

    void setUpperBound(const arma::Col<T>& upperBound);

    void setKernelParameter(const T kernelParameter);
    
    T getKernelParameter();
    
    void setIterations(const unsigned int iterations);
    
    unsigned int getIterations();

  protected:
    void CalculateTrainingKernelMatrix();
    void OptimizeL();

    //RankSVMLearn.cpp Input Variables
    T kernelParameter = 1.0; //sigmaA
    T sigmaPow = 1.0;
    T epsilon = 1.0;
    unsigned int dimension; //N
    arma::Mat<T> parameters;
    unsigned int niter = 1000;
    unsigned int ntraining;
    arma::Col<T> upperBound; //Ci

    //RankSVMLearn.cpp Output Variables
    arma::Col<T> rankingDirection; //optAlphas
    T twoSigmaPow2 = 0;

    //RankSVMLearn.cpp Variables Between Functions
    unsigned int nAlpha; //may need to be int, not sure if ntraining can be 0
    arma::Mat<T> K;
    
    //RankSVMFunc.cpp
    arma::Col<T> fitness; //Fit
    arma::Mat<T> points;
    T ntest;
  };

  //
  // Implementation
  //

  template <typename T>
  RankingSupportVectorMachine<T>::RankingSupportVectorMachine() noexcept {
    upperBound = -arma::ones(1);
  }
  
  //parameters should already be encoded
  template <typename T>
  void RankingSupportVectorMachine<T>::learn(const arma::Mat<T>& parameters) {
    //init
    this->parameters = parameters;
    this->dimension = parameters.n_rows;
    this->ntraining = parameters.n_cols;
    
    nAlpha = ntraining - 1;
    
    rankingDirection = arma::zeros(ntraining - 1);
    K = arma::zeros(ntraining, ntraining);
    
    //calculate standard upperbound if none was set
    if(upperBound(0) == -1 && upperBound.n_elem == 1) {
      upperBound = 10e6*(arma::pow(nAlpha-arma::linspace(0,nAlpha-1,nAlpha),3));
    }
    
    //calc kernel
    CalculateTrainingKernelMatrix();
    
    //optimize alphas
    OptimizeL();
  }

  template <typename T>
  void RankingSupportVectorMachine<T>::CalculateTrainingKernelMatrix() {
    T avrdist = 0;
    for (int i = 0; i < ntraining; i++) {
      K(i, i) = 0;
      for (int j = i + 1; j < ntraining; j++) {
        K(j, i) = std::pow(arma::norm(parameters.col(i) - parameters.col(j)), 2);
        K(i, j) = K(j, i);
        avrdist += std::sqrt(K(i, j));
      }
    }

    avrdist /= (ntraining - 1) * ntraining / 2.0;
    twoSigmaPow2 = 2 * std::pow((kernelParameter * avrdist * sigmaPow), 2);

    //the for loops can be combined and the upper assignments to K never made.
    //which would be much more efficient. If arma can do Mat minus Mat,
    //the upper for-loop is completely unnecessary
    for (int i = 0; i < ntraining; i++) {
      K(i, i) = 1;
      for (int j = i + 1; j < ntraining; j++) {
        K(j, i) = std::exp(-K(j, i) / twoSigmaPow2);
        K(i, j) = K(j, i);
      }
    }
  }

  template <typename T>
  void RankingSupportVectorMachine<T>::OptimizeL() {
    arma::Col<T> sumAlphaDK = arma::zeros(nAlpha);
    arma::Mat<T> divDK = arma::zeros(nAlpha, nAlpha);
    arma::Mat<T> dK = arma::zeros(nAlpha, nAlpha);

    //there is probably a way to combine all these for loops
    for (int i = 0; i < nAlpha; i++) {
      for (int j = 0; j < nAlpha; j++) {
        dK(j, i) = K(j, i) - K(j + 1, i) - K(j, i + 1) + K(j + 1, i + 1);
      }
      rankingDirection(i) = upperBound(i) * (0.95 + 0.05 * arma::randu());
    }

    for (int i = 0; i < nAlpha; i++) {
      T sumAlpha = arma::sum(rankingDirection % dK.col(i));
      sumAlphaDK(i) = (epsilon - sumAlpha) / dK(i, i);
    }

    for (int i = 0; i < nAlpha; i++) {
      for (int j = 0; j < nAlpha; j++) {
        divDK(j, i) = dK(j, i) / dK(j, j);
      }
    }

    for (int i = 0; i < niter; i++) {
      int iMod = i % nAlpha;
      T newAlpha = rankingDirection(iMod) + sumAlphaDK(iMod);
      if (newAlpha > upperBound(iMod)) {
        newAlpha = upperBound(iMod);
      }
      if (newAlpha < 0) {
        newAlpha = 0;
      }
      T deltaAlpha = newAlpha - rankingDirection(iMod);

      T dL = deltaAlpha * dK(iMod, iMod) * (sumAlphaDK(iMod) - 0.5 * deltaAlpha);

      if (dL > 0) {
        sumAlphaDK -= deltaAlpha * divDK.col(iMod);
        rankingDirection(iMod) = newAlpha;
      }
    }
  }

  //points are ntest*nx
  //last 3 arguments probably unnecessary, but cannot guarantee at this point
  //if HCMA doesn't change them in between.
  //parameters should already be encoded
  template <typename T>
  void RankingSupportVectorMachine<T>::evaluatePoints(const arma::Mat<T>& points, const arma::Mat<T>& parameters, const arma::Col<T>& rankingDirection, T twoSigmaPow2) {
    //init
    this->parameters = parameters;
    this->dimension = parameters.n_rows;
    this->ntraining == parameters.n_cols;
    
    this->points = points;
    
    this->rankingDirection = rankingDirection;
    this->twoSigmaPow2 = twoSigmaPow2;
    
    //TODO: might be rows, but if that's the case all others are also probably wrong
    this->ntest = points.n_cols;
    
    this->fitness = arma::zeros(ntest);
    
    //actual evaluation
    arma::Col<T> Kvals = arma::zeros(ntraining);
    for(int i = 0; i < ntest; i++) {
      arma::Col<T> curPoint = points.col(i);
      
      for(int j = 0; j < ntraining; j++) {
        Kvals(j) = std::exp(-arma::sum(
            arma::pow(curPoint-parameters.col(j),2)
            ));
      }
      
      T curFit = 0;
      for(int j = 0; j<ntraining-1;j++) {
        if(rankingDirection(j) != 0) {
          curFit += rankingDirection(j) * (Kvals(j)-Kvals(j+1));
        }
      }
      fitness(i) = curFit;
    }
  }
  
  template <typename T>
  arma::Col<T>& RankingSupportVectorMachine<T>::getRankingDirection() {
    return rankingDirection;
  }
  
  template <typename T>
  arma::Col<T>& RankingSupportVectorMachine<T>::getFitnessForEvaluatedPoints() {
    return fitness;
  }

  template <typename T>
  T RankingSupportVectorMachine<T>::getTwoSigmaPow2() {
    return twoSigmaPow2;
  }

  template <typename T>
  void RankingSupportVectorMachine<T>::setMaximalNumberOfIterations(const unsigned long long maximalNumberOfIterations) {
    this->niter = maximalNumberOfIterations;
  }

  template <typename T>
  void RankingSupportVectorMachine<T>::setUpperBound(const arma::Col<T>& upperBound) {
    this->upperBound = upperBound;
  }
  
  template <typename T>
  void RankingSupportVectorMachine<T>::setKernelParameter(T kernelParameter) {
    this->kernelParameter = kernelParameter;
  }

  template <typename T>
  T RankingSupportVectorMachine<T>::getKernelParameter() {
    return this->kernelParameter;
  }
  
  template <typename T>
  void RankingSupportVectorMachine<T>::setIterations(const unsigned int iterations) {
    this->niter = iterations;
  }
  
  template <typename T>
  unsigned int RankingSupportVectorMachine<T>::getIterations() {
    return niter;
  }
}