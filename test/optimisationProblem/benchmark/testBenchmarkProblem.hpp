#include <catch/catch.hpp>

#include <cstdlib>
#include <array>
#include <memory>
#include <string>

#include <armadillo>

#include <hop>

TEST_CASE("Black box optimisation benchmark problem", "[benchmark]") {

  SECTION("Test impementation correctness") {
    for(auto numberOfDimensions : {2, 40}) {
      std::array<std::unique_ptr<hop::BenchmarkProblem>, 24> bechmarkProblems = {
        std::unique_ptr<hop::BenchmarkProblem>(new hop::SphereFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::EllipsoidalFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::RastriginFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::BuecheRastriginFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::LinearSlope(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::AttractiveSectorFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::StepEllipsoidalFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::RosenbrockFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::RosenbrockFunctionRotated(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::EllipsoidalFunctionRotated(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::DiscusFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::BentCigarFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::SharpRidgeFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::DifferentPowersFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::RastriginFunctionRotated(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::WeierstrassFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::SchaffersF7Function(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::SchaffersF7FunctionIllConditioned(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::CompositeGriewankRosenbrockFunctionF8F2(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::SchwefelFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::GallaghersGaussian101mePeaksFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::GallaghersGaussian21hiPeaksFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::KatsuuraFunction(numberOfDimensions)),
        std::unique_ptr<hop::BenchmarkProblem>(new hop::LunacekBiRastriginFunction(numberOfDimensions))
      };

      arma::Mat<double> parameters;
      parameters.load("../test/data/parameters," + std::to_string(numberOfDimensions) +".mat");

      arma::Col<double> translation;
      translation.load("../test/data/translation," + std::to_string(numberOfDimensions) +".mat");

      arma::Col<double> one;
      one.load("../test/data/one," + std::to_string(numberOfDimensions) +".mat");

      arma::Mat<double> rotationR;
      rotationR.load("../test/data/rotationR," + std::to_string(numberOfDimensions) +".mat");

      arma::Mat<double> rotationQ;
      rotationQ.load("../test/data/rotationQ," + std::to_string(numberOfDimensions) +".mat");

      arma::Mat<double> deltaC101;
      deltaC101.load("../test/data/deltaC101," + std::to_string(numberOfDimensions) +".mat");

      arma::Mat<double> localOptimaY101;
      localOptimaY101.load("../test/data/localOptimaY101," + std::to_string(numberOfDimensions) +".mat");

      arma::Mat<double> deltaC21;
      deltaC21.load("../test/data/deltaC21," + std::to_string(numberOfDimensions) +".mat");

      arma::Mat<double> localOptimaY21;
      localOptimaY21.load("../test/data/localOptimaY21," + std::to_string(numberOfDimensions) +".mat");

      for (std::size_t n = 0; n < bechmarkProblems.size(); n++) {
        arma::Col<double> expected;
        expected.load("../test/data/expectedF" + std::to_string(n + 1) + "," + std::to_string(numberOfDimensions) +".mat");

        bechmarkProblems.at(n)->setMaximalNumberOfEvaluations(parameters.n_cols);
        bechmarkProblems.at(n)->setObjectiveValueTranslation(0);
        bechmarkProblems.at(n)->setTranslation(translation);
        bechmarkProblems.at(n)->setOne(one);
        bechmarkProblems.at(n)->setRotationR(rotationR);
        bechmarkProblems.at(n)->setRotationQ(rotationQ);
        bechmarkProblems.at(n)->setDeltaC101(deltaC101);
        bechmarkProblems.at(n)->setLocalOptimaY101(localOptimaY101);
        bechmarkProblems.at(n)->setDeltaC21(deltaC21);
        bechmarkProblems.at(n)->setLocalOptimaY21(localOptimaY21);

        for (std::size_t j = 0; j < parameters.n_cols; j++) {
          REQUIRE(bechmarkProblems.at(n)->getObjectiveValue(parameters.col(j)) == Approx(expected.at(j)));
        }
      }
    }
  }
}

