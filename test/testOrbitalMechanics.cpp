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
  GIVEN("Timestamp in mjd format and a 7-element vector with keplerian elements and a reference date in mjd format") {
    arma::Col<double>::fixed<7> keplerianElementsVenusData = {7.233268496749391e-01, 6.755697267164094e-03, 3.394589632336535e+00, 5.518541455452200e+01, 7.667837563371675e+01, 4.931425178852966e+01, 51544.0};

    WHEN("Valid timestamp and keplerian elements are given") {
      THEN("Return position and velocity vector") {
        double mjd = 51544.0;
        arma::Col<double>::fixed<3> expectedPosition = {-1.07559e+11, -3.41113e+09, 5.12246e+09};
        arma::Col<double>::fixed<3> expectedVelocity = {874.506, -35141.9, -1232.67};

        std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> resultPair = mant::itd::orbitOnPosition(mjd, keplerianElementsVenusData);

        IS_EQUAL(resultPair.first, expectedPosition);
        IS_EQUAL(resultPair.second, expectedVelocity);
      }
    }

    WHEN("The modifiedJulianDay is infinite") {
      THEN("Throw a std::logic_error") {
        double mjd = arma::datum::inf;

        CHECK_THROWS_AS(mant::itd::orbitOnPosition(mjd, keplerianElementsVenusData), std::logic_error);
      }
    }

    WHEN("Eccentricity is 1.0 or greater") {
      THEN("Throw a std::logic_error") {
        arma::Col<double>::fixed<7> keplerianElementsVenusDataHighEccentricity = {7.233268496749391e-01, 1.6755697267164094, 3.394589632336535e+00, 5.518541455452200e+01, 7.667837563371675e+01, 4.931425178852966e+01, 51544.0};
        double mjd = 54000.0 + continuousRandomNumber();

        CHECK_THROWS_AS(mant::itd::orbitOnPosition(mjd, keplerianElementsVenusDataHighEccentricity), std::logic_error);
      }
    }
  }
}
