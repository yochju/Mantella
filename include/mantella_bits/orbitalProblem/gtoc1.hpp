#pragma once

// Armadillo
#include <armadillo>

// Mantella
#include "mantella_bits/orbitalProblem.hpp"

namespace mant {
  namespace itd {

    class GTOC1 : public OrbitalProblem {
     public:
      GTOC1(
          const std::vector<arma::Col<double>::fixed<7>>& orbitalTargetSequence);

      arma::Col<double> objectiveFunction(
          const arma::Col<double>& parameter);

     private:
      double rocketMass_ = 1500.0; // Satellite initial mass [Kg]
      double rocketSpecificImpulse_ = 2500.0; // Satellite specific impulse [s]
      double rocketDVlaunch_ = 2500; // Launcher DV in m/s

      std::vector<arma::Col<double>::fixed<7>> orbitalTargetSequence_;
    };
  }
}