#pragma once

// Armadillo
#include <armadillo>

namespace mant {
  namespace itd {
    
    class OrbitalProblem {
      protected:
        explicit OrbitalProblem();
        
        virtual arma::Col<double> problemFunction(arma::Col<double> parameter) = 0;
      
    };
  }
}