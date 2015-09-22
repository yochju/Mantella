#include <mantella_bits/helper/rankingSupportVectorMachine.hpp>

#include <mantella_bits/helper/rng.hpp>

namespace mant {
  
  //parameters should already be encoded
  //Encoding function from HCMAs code will be in HCMA code and not here
  RankingSupportVectorMachine::RankingSupportVectorMachine(const arma::Mat<double>& parameters, arma::uword niter, double epsilon, arma::Col<double> upperBound, double kernelParameter, double sigmaPow) noexcept {  
    //init
    this->parameters = parameters;
    this->dimension = parameters.n_rows;
    this->ntraining = parameters.n_cols;
    
    nAlpha = ntraining - 1;
    
    rankingDirection = arma::zeros(ntraining - 1);
    K = arma::zeros(ntraining, ntraining);

    //calculate standard upperbound if none was set
    if (upperBound(0) == -1 && upperBound.n_elem == 1) {
      upperBound = 10e6 * (arma::pow(nAlpha - arma::linspace(0, nAlpha - 1, nAlpha), 3));
    }
  }

  void RankingSupportVectorMachine::learn() {
    //calc kernel
    //HCMA: CalculateTrainingKernelMatrix()
    double avrdist = 0;
    for (arma::uword i = 0; i < ntraining; i++) {
      K(i, i) = 0;
      for (arma::uword j = i + 1; j < ntraining; j++) {
        double norm = arma::norm(parameters.col(i) - parameters.col(j));
        K(j, i) = std::pow(norm, 2);
        K(i, j) = K(j, i);
        avrdist += norm;
      }
    }

    avrdist /= (ntraining - 1) * ntraining / 2.0;
    twoSigmaPow2 = 2 * std::pow((kernelParameter * avrdist * sigmaPow), 2);

    //the for loops can be combined and the upper assignments to K never made.
    //which would be much more efficient. If arma can do Mat minus Mat,
    //the upper for-loop is completely unnecessary
    for (arma::uword i = 0; i < ntraining; i++) {
      K(i, i) = 1;
      for (arma::uword j = i + 1; j < ntraining; j++) {
        K(j, i) = std::exp(-K(j, i) / twoSigmaPow2);
        K(i, j) = K(j, i);
      }
    }

    //optimize alphas
    //HCMA: OptimizeL()
    arma::Col<double> sumAlphaDK = arma::zeros(nAlpha);
    arma::Mat<double> divDK = arma::zeros(nAlpha, nAlpha);
    arma::Mat<double> dK = arma::zeros(nAlpha, nAlpha);

    //there is probably a way to combine all these for loops
    for (arma::uword i = 0; i < nAlpha; i++) {
      for (arma::uword j = 0; j < nAlpha; j++) {
        dK(j, i) = K(j, i) - K(j + 1, i) - K(j, i + 1) + K(j + 1, i + 1);
      }
      rankingDirection(i) = upperBound(i) * (0.95 + 0.05 * std::uniform_real_distribution<double>(0.0, 1.0)(Rng::getGenerator()));
    }

    for (arma::uword i = 0; i < nAlpha; i++) {
      double sumAlpha = arma::dot(rankingDirection, dK.col(i));
      sumAlphaDK(i) = (epsilon - sumAlpha) / dK(i, i);
    }

    for (arma::uword i = 0; i < nAlpha; i++) {
      for (arma::uword j = 0; j < nAlpha; j++) {
        divDK(j, i) = dK(j, i) / dK(j, j);
      }
    }

    for (arma::uword i = 0; i < niter; i++) {
      int iMod = i % nAlpha;
      double newAlpha = rankingDirection(iMod) + sumAlphaDK(iMod);
      if (newAlpha > upperBound(iMod)) {
        newAlpha = upperBound(iMod);
      }
      if (newAlpha < 0) {
        newAlpha = 0;
      }
      double deltaAlpha = newAlpha - rankingDirection(iMod);

      double dL = deltaAlpha * dK(iMod, iMod) * (sumAlphaDK(iMod) - 0.5 * deltaAlpha);

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
  void RankingSupportVectorMachine::evaluatePoints(const arma::Mat<double>& points, const arma::Mat<double>& parameters, const arma::Col<double>& rankingDirection, double twoSigmaPow2) {
    //init
    this->parameters = parameters;
    this->dimension = parameters.n_rows;
    this->ntraining = parameters.n_cols;

    this->points = points;

    this->rankingDirection = rankingDirection;
    this->twoSigmaPow2 = twoSigmaPow2;

    //TODO: might be rows, but if that's the case all others are also probably wrong
    this->ntest = points.n_cols;

    this->fitness = arma::zeros(ntest);

    //actual evaluation
    arma::Col<double> Kvals = arma::zeros(ntraining);
    for (int i = 0; i < ntest; i++) {
      arma::Col<double> curPoint = points.col(i);

      for (arma::uword j = 0; j < ntraining; j++) {
        Kvals(j) = std::exp(-arma::sum(
            arma::pow(curPoint - parameters.col(j), 2)
            ));
      }

      double curFit = 0;
      for (arma::uword j = 0; j < ntraining - 1; j++) {
        if (rankingDirection(j) != 0) {
          curFit += rankingDirection(j) * (Kvals(j) - Kvals(j + 1));
        }
      }
      fitness(i) = curFit;
    }
  }

  arma::Col<double>& RankingSupportVectorMachine::getRankingDirection() {
    return rankingDirection;
  }

  arma::Col<double>& RankingSupportVectorMachine::getFitnessForEvaluatedPoints() {
    return fitness;
  }

  double RankingSupportVectorMachine::getTwoSigmaPow2() {
    return twoSigmaPow2;
  }

  void RankingSupportVectorMachine::setMaximalNumberOfIterations(const unsigned long long maximalNumberOfIterations) {
    this->niter = maximalNumberOfIterations;
  }

  void RankingSupportVectorMachine::setUpperBound(const arma::Col<double>& upperBound) {
    this->upperBound = upperBound;
  }

  void RankingSupportVectorMachine::setKernelParameter(double kernelParameter) {
    this->kernelParameter = kernelParameter;
  }

  double RankingSupportVectorMachine::getKernelParameter() {
    return this->kernelParameter;
  }

  void RankingSupportVectorMachine::setIterations(const arma::uword iterations) {
    this->niter = iterations;
  }

  arma::uword RankingSupportVectorMachine::getIterations() {
    return niter;
  }
}