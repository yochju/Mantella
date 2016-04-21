#include "mantella_bits/orbitalMechanics.hpp"

// C++ standard library
#include <cmath>

// Mantella
#include "mantella_bits/assert.hpp"
#include "mantella_bits/numericalAnalysis.hpp"
#include "mantella_bits/optimisationProblem.hpp"
#include "mantella_bits/optimisationAlgorithm/hookeJeevesAlgorithm.hpp"

namespace mant {
  namespace itd {

    std::pair<double, double> gravityAssist(
        const arma::Col<double>::fixed<3>& inboundVelocity,
        const arma::Col<double>::fixed<3>& outboundVelocity) {
      double inboundVelocityLength = arma::norm(inboundVelocity);
      double outboundVelocityLength = arma::norm(outboundVelocity);

      double velocitiesDotProduct = arma::dot(inboundVelocity, outboundVelocity);
      double alpha = std::acos(velocitiesDotProduct / (inboundVelocityLength * outboundVelocityLength));

      double inboundAcceleration = 1.0 / std::pow(inboundVelocityLength, 2.0);
      double outboundAcceleration = 1.0 / std::pow(outboundVelocityLength, 2.0);

      double rp = mant::brent(
          [&inboundAcceleration, &outboundAcceleration, &alpha](double parameter) { 
              return std::asin(inboundAcceleration / (inboundAcceleration + parameter)) + std::asin(outboundAcceleration / (outboundAcceleration + parameter)) - alpha;
          },
          0.0, 1.0, 100, 1e-10);

      double DV = std::abs(std::sqrt(std::pow(outboundVelocityLength, 2) + 2.0 / rp) - std::sqrt(std::pow(inboundVelocityLength, 2) + 2.0 / rp));

      return {DV, rp};
    }

    std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> lambert(
        const arma::Col<double>::fixed<3>& departurePosition,
        const arma::Col<double>::fixed<3>& arrivalPosition,
        const double transferTime) {
      verify(transferTime > 0.0, "lambert: The transfer time must be greater than zero.");

      //std::cout << "START LAMBERT" << std::endl;

      double departureDistanceFromSun = arma::norm(departurePosition);
      double arrivalDistanceFromSun = arma::norm(arrivalPosition);
      double depatureArrivalDotProduct = arma::dot(departurePosition, arrivalPosition);

      arma::uword maximalNumberOfRevolutions = 0;
      arma::uword numberOfRevolutions = 0;

      double t_m = 1.0;
      double A = t_m * std::sqrt(departureDistanceFromSun * arrivalDistanceFromSun * (1.0 + depatureArrivalDotProduct / (departureDistanceFromSun * arrivalDistanceFromSun)));
      double B = 0.0;

      auto timeOfFlightFunction = [&A, &B, &t_m, &departureDistanceFromSun, &arrivalDistanceFromSun, &transferTime, &numberOfRevolutions](
          const double parameter) {
        if((numberOfRevolutions == 0 && (std::abs(parameter + 4.0 * arma::datum::pi) < 1e-5 || std::abs(parameter - 4.0 * std::pow(arma::datum::pi, 2.0)) < 1e-5))
        || (numberOfRevolutions > 0 && (std::abs(parameter - 4.0 * std::pow(numberOfRevolutions * arma::datum::pi, 2.0)) < 1e-5 || std::abs(parameter - 4.0 * std::pow((numberOfRevolutions + 1.0) * arma::datum::pi, 2.0)) < 1e-5))) { 
          return arma::datum::inf;
        }
          
        A = t_m * A;
        //std::cout << "A: " << A << "  -  ";
                
        // Calculates the values for the second and third Stumpff function c2 and c3.
        double c2;
        double c3;
        if (parameter > 0) {
          c2 = (1.0 - std::cos(std::sqrt(parameter))) / parameter;
          c3 = (1.0 - std::sin(std::sqrt(parameter)) / std::sqrt(parameter)) / parameter;
        } else if (parameter < 0) {
          c2 = (1.0 - std::cosh(std::sqrt(-parameter))) / parameter;
          c3 = (1.0 - std::sinh(std::sqrt(-parameter)) / std::sqrt(-parameter)) / parameter;
        } else {
          c2 = 1.0/2.0;
          c3 = 1.0/6.0;
        }
        
        //std::cout << "c2: " << c2 << "  -  ";
        //std::cout << "c3: " << c3 << "  -  ";
        
        B = departureDistanceFromSun + arrivalDistanceFromSun + 1.0 / std::sqrt(c2) * (A * (parameter * c3 - 1.0));
        
        //std::cout << "B: " << B << "  -  ";
        //std::cout << "val: " << 1.0 / std::sqrt(standardGravitationalParameterOfSun) * (std::pow(std::sqrt(B / c2), 3.0) * c3 + A * std::sqrt(B)) / 86400.0 - transferTime << std::endl;
                 
        return (1.0 / std::sqrt(standardGravitationalParameterOfSun)) * (std::pow(std::sqrt(B / c2), 3.0) * c3 + A * std::sqrt(B)) - transferTime;
      };

      //std::cout << "function: " << timeOfFlightFunction(3.0 * std::pow(arma::datum::pi, 2.0)) << std::endl;

      auto calculateVelocityVectorsFunction = [&A, &B, &departurePosition, &departureDistanceFromSun, &arrivalPosition, &arrivalDistanceFromSun]() -> std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> {
        //std::cout << "begin calculateVelocityVectorsFunction" << std::endl;
        double F = 1.0 - B / departureDistanceFromSun;
        double G = A * std::sqrt(B / standardGravitationalParameterOfSun);
        double Gdot = 1.0 - B / arrivalDistanceFromSun;
        //std::cout << "end calculateVelocityVectorsFunction" << std::endl;
        return {(1.0 / G) * (arrivalPosition - departurePosition * F), (1.0 / G) * (Gdot * arrivalPosition - departurePosition)};
      };
      //std::cout << "pre decl" << std::endl;
      std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> velocityPairs;
      velocityPairs.reserve(maximalNumberOfRevolutions * 4 + 2);
      //std::cout << "post decl" << std::endl;
      double lowerBound = -4.0 * arma::datum::pi;
      double upperBound = 4.0 * std::pow(arma::datum::pi, 2.0);

      // For zero revolutions
      //std::cout << "0" << std::endl;
      double brentShortWaySolution = mant::brent(timeOfFlightFunction, lowerBound, upperBound, 100, 1e-10);
      //std::cout << "pre push_back" << std::endl;
      //std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> varr = calculateVelocityVectorsFunction();
      //std::cout << "varr type: " << typeid(varr).name() << std::endl;
      velocityPairs.push_back(calculateVelocityVectorsFunction());
      //std::cout << "post push_back" << std::endl;
      t_m = -t_m;
      //std::cout << "1" << std::endl;
      double brentLongWaySolution = mant::brent(timeOfFlightFunction, lowerBound, upperBound, 100, 1e-10);
      velocityPairs.push_back(calculateVelocityVectorsFunction());
      t_m = -t_m;

      mant::OptimisationProblem timeOfFlightMinimumProblem(1);
      timeOfFlightMinimumProblem.setObjectiveFunction([&timeOfFlightFunction](const arma::Col<double>& parameter) {
        //std::cout << "in setObjectiveFunction" << std::endl;
        return timeOfFlightFunction(parameter(0));
      });

      mant::HookeJeevesAlgorithm timeOfFlightMinimumAlgorithm;
      timeOfFlightMinimumAlgorithm.setMaximalNumberOfIterations(100);

      for (numberOfRevolutions = 1; numberOfRevolutions <= maximalNumberOfRevolutions; numberOfRevolutions++) {
        lowerBound = 4.0 * std::pow((numberOfRevolutions)*arma::datum::pi, 2.0);
        upperBound = 4.0 * std::pow((numberOfRevolutions + 1) * arma::datum::pi, 2.0);
        //std::cout << "(lower, upper) = (" << lowerBound << ", " << upperBound << ")" << std::endl;
        timeOfFlightMinimumProblem.setLowerBounds({lowerBound});
        //std::cout << "between bounds" << std::endl;
        timeOfFlightMinimumProblem.setUpperBounds({upperBound});
        //std::cout << "post bounds" << std::endl;
        timeOfFlightMinimumAlgorithm.optimise(timeOfFlightMinimumProblem, {(upperBound - lowerBound) / 2.0});
        //std::cout << "after optimise" << std::endl;
        arma::Col<double> timeOfFlightMinimum = timeOfFlightMinimumAlgorithm.getBestParameter();
        //std::cout << "2 - " << numberOfRevolutions << std::endl;
        double brentShortWaySolutionLeft = mant::brent(timeOfFlightFunction, lowerBound, timeOfFlightMinimum(0), 100, 1e-10);

        if (std::isnan(brentShortWaySolutionLeft)) {
          break;
        }

        std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> brentShortWaySolutionLeftVelocities = calculateVelocityVectorsFunction();
        //std::cout << "3 - " << numberOfRevolutions << std::endl;
        double brentShortWaySolutionRight = mant::brent(timeOfFlightFunction, timeOfFlightMinimum(0), upperBound, 100, 1e-10);
        std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> brentShortWaySolutionRightVelocities = calculateVelocityVectorsFunction();
        t_m = -t_m;

        timeOfFlightMinimumAlgorithm.optimise(timeOfFlightMinimumProblem, {(upperBound - lowerBound) / 2.0});
        timeOfFlightMinimum = timeOfFlightMinimumAlgorithm.getBestParameter();
        //std::cout << "4 - " << numberOfRevolutions << std::endl;
        double brentLongWaySolutionLeft = mant::brent(timeOfFlightFunction, lowerBound, timeOfFlightMinimum(0), 100, 1e-10);

        if (std::isnan(brentLongWaySolutionLeft)) {
          break;
        }

        std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> brentLongWaySolutionLeftVelocities = calculateVelocityVectorsFunction();
        //std::cout << "5 - " << numberOfRevolutions << std::endl;
        double brentLongWaySolutionRight = mant::brent(timeOfFlightFunction, timeOfFlightMinimum(0), upperBound, 100, 1e-10);
        std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> brentLongWaySolutionRightVelocities = calculateVelocityVectorsFunction();

        t_m = -t_m;

        velocityPairs.push_back(brentShortWaySolutionLeftVelocities);
        velocityPairs.push_back(brentShortWaySolutionRightVelocities);
        velocityPairs.push_back(brentLongWaySolutionLeftVelocities);
        velocityPairs.push_back(brentLongWaySolutionRightVelocities);
      }

      return velocityPairs;
    }

    std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> orbitOnPosition(
        const double modifiedJulianDay,
        const arma::Col<double>::fixed<7>& keplerianElements) {
      if (!std::isfinite(modifiedJulianDay)) {
        // TODO
        verify(std::isfinite(modifiedJulianDay), "orbitOnPosition: The modifiedJulianDay must be finite.");
      }

      double timePassedSincePeriapsis = (modifiedJulianDay - keplerianElements(6)) * 86400.0; // in sec

      double semiMajorAxis = keplerianElements(0) * 149597870691.0; //in m
      double eccentricity = keplerianElements(1);
      double inclination = (arma::datum::pi / 180.0) * keplerianElements(2);
      double omg = (arma::datum::pi / 180.0) * keplerianElements(3);
      double omp = (arma::datum::pi / 180.0) * keplerianElements(4);
      double ea = (arma::datum::pi / 180.0) * keplerianElements(5);

      if (eccentricity >= 1.0) {
        // TODO
        verify(eccentricity < 1.0, "orbitOnPosition: The eccentricity must be lesser than 1.0.");
      }
      double meanMotion = std::sqrt(standardGravitationalParameterOfSun / std::pow(semiMajorAxis, 3.0));
      ea += timePassedSincePeriapsis * meanMotion;

      //m2e begin
      ea = mant::brent([&ea, &eccentricity](
                           double parameter) { 
        return parameter - eccentricity * std::sin(parameter) - ea;
      },
          ea + eccentricity * std::cos(ea) - 1.0, ea + eccentricity * std::cos(ea) + 1.0, 100, 1e-10); //TODO bounds round about (ea + eccentricity * std::cos(ea)) +- 1.0 okay?
      //m2e end

      //par2ic begin
      double b = semiMajorAxis * std::sqrt(1.0 - std::pow(eccentricity, 2.0));
      double n = std::sqrt(standardGravitationalParameterOfSun / std::pow(semiMajorAxis, 3.0));
      double xper = semiMajorAxis * (std::cos(ea) - eccentricity);
      double yper = b * std::sin(ea);
      double xdotper = -(semiMajorAxis * n * std::sin(ea)) / (1.0 - eccentricity * std::cos(ea));
      double ydotper = (b * n * std::cos(ea)) / (1.0 - eccentricity * std::cos(ea));

      double cosomg = std::cos(omg);
      double cosomp = std::cos(omp);
      double sinomg = std::sin(omg);
      double sinomp = std::sin(omp);
      double cosi = std::cos(inclination);
      double sini = std::sin(inclination);

      arma::Mat<double>::fixed<3, 3> rotationMatrix = {{cosomg * cosomp - sinomg * sinomp * cosi, -cosomg * sinomp - sinomg * cosomp * cosi, sinomg * sini},
          {sinomg * cosomp + cosomg * sinomp * cosi, -sinomg * sinomp + cosomg * cosomp * cosi, -cosomg * sini},
          {sinomp * sini, cosomp * sini, cosi}};

      arma::Col<double>::fixed<3> temp = {xper, yper, 0.0};
      arma::Col<double>::fixed<3> temp2 = {xdotper, ydotper, 0.0};

      arma::Col<double>::fixed<3> positionVector = rotationMatrix * temp;
      arma::Col<double>::fixed<3> velocityVector = rotationMatrix * temp2;

      return std::make_pair(positionVector, velocityVector);
      //par2ic end
    }
  }
}
