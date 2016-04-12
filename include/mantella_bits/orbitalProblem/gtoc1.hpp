#pragma once

// C++ standard library
#include <vector>

// Armadillo
#include <armadillo>

// Mantella
#include "mantella_bits/orbitalProblem.hpp"

namespace mant {
  namespace itd {
    
    class GTOC1 : public OrbitalProblem {
     public:
      using OrbitalProblem::OrbitalProblem;    
      GTOC1(std::vector<arma::Col<double>::fixed<6>> orbitalTargetSequence);
      
      arma::Col<double> problemFunction(arma::Col<double> parameter);
      
      private:
        double rocketMass_ = 1500.0; // Satellite initial mass [Kg]
        double rocketSpecificImpulse_ = 2500.0; // Satellite specific impulse [s]
        double rocketDVlaunch_ = 2.5; // Launcher DV in km/s
        
        arma::Col<double>::fixed<3> asteroidKeplerValues_ = {2.5897261, 0.2734625, 6.40734, 128.34711, 264.78691, 320.479555};
        
        std::vector<arma::Col<double>::fixed<6>> orbitalTargetSequence_;
    };
  }
}