#include <hop_bits/propertiesAnalysis/functionModelAnalysis.hpp>

namespace hop {
  arma::Col<double> FunctionModelAnalysis::getResiduals() const {
    return residuals_;
  }
}
