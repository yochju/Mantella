#pragma once

// C++ standard library
#include <utility>

// Armadillo
#include <armadillo>

namespace mant {
  namespace itd {

    const double standardGravitationalParameterOfSun = 1.32712440018e20;

    std::pair<double, double> gravityAssist(
        const arma::Col<double>::fixed<3>& inboundVelocity,
        const arma::Col<double>::fixed<3>& outboundVelocity);

    std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> lambert(
        const arma::Col<double>::fixed<3>& departurePosition,
        const arma::Col<double>::fixed<3>& arrivalPosition,
        const double transferTime);

    std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> orbitOnPosition(
        const double modifiedJulianDay2000,
        const arma::Mat<double>::fixed<2, 6>& keplerianElements);
  }
}
