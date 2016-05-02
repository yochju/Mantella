#pragma once

// Armadillo
#include <armadillo>

namespace mant {
  namespace itd {
    class OrbitalProblem {
     protected:
      virtual arma::Col<double> objectiveFunction(const arma::Col<double>& parameter) = 0;
    };
  }
}