#pragma once

namespace mant {
  namespace itd {
    class OrbitalProblem {
      protected:
        explicit OrbitalProblem();
        
        virtual arma::Col<double> problemFunction() = 0;
      
    };
  }
}