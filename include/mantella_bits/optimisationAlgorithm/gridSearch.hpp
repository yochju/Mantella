#pragma once

// C++ standard library
#include <memory>
#include <string>

// Armadillo
#include <armadillo>

// Mantella
#include "mantella_bits/optimisationAlgorithm.hpp"
namespace mant {
  class OptimisationProblem;
}

namespace mant {
  class GridSearch : public OptimisationAlgorithm {
   public:
    explicit GridSearch();

    void optimise(
        OptimisationProblem& optimisationProblem,
        const arma::Mat<double>& initialParameters) override;

    void optimise(
        OptimisationProblem& optimisationProblem);
    
    void setNumberOfSamplesPerDimension(
        const arma::Col<arma::uword>& numberOfSamplesPerDimension);

   protected:
    arma::Col<arma::uword> numberOfSamplesPerDimension_;
    
    std::vector<arma::Col<double>> samples_;
    arma::Col<arma::uword> sampleIndicies_;
  };
}