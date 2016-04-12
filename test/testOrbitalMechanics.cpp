// Armadillo
#include <armadillo>

// Catch
#include <catch.hpp>
#include "catchExtension.hpp"

// Mantella
#include <mantella>

SCENARIO("gravityAssist", "[orbitalMechanics][gravityAssist]") {
  GIVEN("An inbound and outbound velocity vector") {
    WHEN("Two different vectors are given") {
      THEN("Return std::pair with velocity difference and rp") {
        arma::Col<double>::fixed<3> inboundVector = {123.0, 246.0, 369.0};
        arma::Col<double>::fixed<3> outboundVector = {234.0, 135.0, 470.0};
        std::pair<double, double> resultPair = mant::itd::gravityAssist(inboundVector, outboundVector);
        CHECK(resultPair.first == Approx(69.2261));
        CHECK(resultPair.second == Approx(2.0046e-05));
      }
    }

    WHEN("Inbound and outbound vector are equal") {
      THEN("Throw a std::logic_error") {
        arma::Col<double>::fixed<3> inboundVector = {123.0, 246.0, 369.0};
        arma::Col<double>::fixed<3> outboundVector = {123.0, 246.0, 369.0};

        CHECK_THROWS_AS(mant::itd::gravityAssist(inboundVector, outboundVector), std::logic_error);
      }
    }

    WHEN("Inbound and outbound vector have infinty values") {
      THEN("Throw a std::logic_error") {
        arma::Col<double>::fixed<3> inboundVector = {arma::datum::inf, arma::datum::inf, arma::datum::inf};
        arma::Col<double>::fixed<3> outboundVector = {arma::datum::inf, arma::datum::inf, arma::datum::inf};

        CHECK_THROWS_AS(mant::itd::gravityAssist(inboundVector, outboundVector), std::logic_error);
      }
    }
  }
}

SCENARIO("lambert", "[orbitalMechanics][lambert]") {
  GIVEN("") {
    WHEN("") {
      THEN("") {
      }
    }

    WHEN("") {
      THEN("Throw a std::logic_error") {
      }
    }
  }
}

SCENARIO("orbitOnPosition", "[orbitalMechanics][orbitOnPosition]") {
  GIVEN("Timestamp in mjd2000 format and (2,6)-Matrix with keplerian elements") {
    arma::Mat<double>::fixed<2, 6> keplerianElementsVenus =
        {{0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255},
         {0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418}};

    WHEN("Valid timestamp and keplerian elements are given") {
      THEN("Return position and velocity vector") {
        arma::Col<double>::fixed<3> expectedPosition = {1.08442e+11, -4.42187e+09, -6.31965e+09};
        arma::Col<double>::fixed<3> expectedVelocity = {1272.54, 34831.5, 402.894};

        std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> resultPair = mant::itd::orbitOnPosition(1234.0, keplerianElementsVenus);

        IS_EQUAL(resultPair.first, expectedPosition);
        IS_EQUAL(resultPair.second, expectedVelocity);
      }
    }

    WHEN("Timestamp is not in range (-73048.0, 18263.0)") {
      THEN("Throw a std::logic_error") {
        double mjd2000 = -73048.0 - std::abs(continuousRandomNumber());
        CHECK_THROWS_AS(mant::itd::orbitOnPosition(mjd2000, keplerianElementsVenus), std::logic_error);

        mjd2000 = 18263.0 + std::abs(continuousRandomNumber());
        CHECK_THROWS_AS(mant::itd::orbitOnPosition(mjd2000, keplerianElementsVenus), std::logic_error);

        mjd2000 = arma::datum::inf;
        CHECK_THROWS_AS(mant::itd::orbitOnPosition(mjd2000, keplerianElementsVenus), std::logic_error);
      }
    }

    WHEN("Eccentricity is 1.0 or greater") {
      THEN("Throw a std::logic_error") {
        arma::Mat<double>::fixed<2, 6> keplerianElementsVenusHighEccentricity =
            {{0.72333566, 1.20677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255},
             {0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418}};

        double mjd2000 = std::abs(continuousRandomNumber());
        CHECK_THROWS_AS(mant::itd::orbitOnPosition(mjd2000, keplerianElementsVenusHighEccentricity), std::logic_error);
      }
    }
  }
}
