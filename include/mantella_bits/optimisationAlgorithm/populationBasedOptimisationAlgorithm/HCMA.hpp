namespace mant {

  template <typename T>
  class HCMA : public PopulationBasedOptimisationAlgorithm<T> {
  public:
    explicit HCMA(
        const std::shared_ptr<OptimisationProblem<T>> optimisationProblem,
        const uint populationSize) noexcept;

    std::string toString() const noexcept override;
  protected:

    uint lambdaMult; //settings.lambdaMult
    double sigma0 = 2.0; //later found in cur_state.insigma
    uint irun = 0; //cur_state.irun
    uint restarts = 9; //opts.Restarts ; nrestarts
    bool bipopCriterion = 0; //bipop_criterion

    //variables for STEP
    uint STEPevals = 0; //STEPalgo.evals
    uint STEPmaxEvals = 0; //STEPalgo.perc * fgeneric('evaluations')

    //variables for the Surrogate
    bool SURRuseSurrogate = true;
    uint SURRmaxEvals = 0; //settings.MaxEvalsWithSurrogate

    void optimiseImplementation() override;
  };

  //
  //Implementation
  //

  template <typename T>
  HCMA<T>::HCMA(const std::shared_ptr<OptimisationProblem<T>> optimisationProblem,
      const uint populationSize) noexcept
  : PopulationBasedOptimisationAlgorithm<T>(optimisationProblem, populationSize) {
    if (this->numberOfDimensions_ < 10) {
      lambdaMult = 1;
    } else {
      lambdaMult = std::floor(std::pow(10, 1 + std::log2(this->numberOfDimensions_ * 1.0 / 10)));
    }

  }

  template <typename T>
  void HCMA<T>::optimiseImplementation() {
    while ((irun <= restarts || bipopCriterion) && !isFinished()) {
      bool useRealFitFunc = true; //algo.realfunc
      uint iSZ = 0; //algo.iSZ TODO: unsure what iSZ stands for.
      uint aSZ = 0; //algo.aSZ TODO: unsure what aSZ stands for.

      //optModel
      uint numberOfTrainingPoints = 40 + std::floor(4 * (std::pow(numberOfDimensions_, 1.7))); //coeff(1)
      uint cvalCoeff = 6; //coeff(2)
      uint ciCoeff = 3; //coeff(3)
      double sigmaA = 1; //coeff(4)
      uint coeff5 = 1000; //coeff(5)
      arma::Col<double> lowerBoundary = {4 * numberOfDimensions_, 0, 0, 0.5, 100}; //xmin
      arma::Col<double> upperBoundary = {2 * 40 + std::floor(4 * (std::pow(numberOfDimensions_, 1.7))), 10, 6, 2, 1500}; //xmax
      uint maxTrainingPoints = upperBoundary(0); //MaxTrainingPoints

      uint maxArchSize = maxTrainingPoints; //algo.maxArchSize
      arma::Mat<double> ARX = arma::zeros(numberOfDimensions_, maxArchSize); //algo.ARX
      arma::Col<double> ARF = arma::zeros(maxArchSize); //algo.ARF
      arma::Col<double> ARneval = arma::zeros(maxArchSize); //algo.ARneval

      irun++;

      //TODO: other config stuff here

      uint stop = 0; //TODO: bool might suffice
      while ((stop == 0) && !isFinished()) {
        //STEP
        //TODO: insert rest of STEP code
        while (STEPevals < STEPmaxEvals && !isFinished()) {
          for (int i = 0; < numberOfDimensions_; i++) {

          }
        }
        stop = 1; //remove line when code is here
        //STEP finished

        if (this->numberOfIterations_ > SURRmaxEvals) {
          SURRuseSurrogate = false;
        }
        
        if(SURRuseSurrogate == false) {//original CMA
          
        } else {//evaluate lambda points on some iSTEP'th step with adaptation of iSTEP
          
        }
      }
    }
  }

  template <typename T>
  std::string HCMA<T>::toString() const noexcept {
    return "cuckoo_search";
  }
}