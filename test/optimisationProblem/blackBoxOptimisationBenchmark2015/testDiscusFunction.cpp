// Catch
#include <catch.hpp>

// C++ Standard Library
#include <cstdlib>
#include <string>

// Armadillo
#include <armadillo>

// Boost
#include <boost/filesystem.hpp>

// HOP
#include <hop>

extern boost::filesystem::path testDirectory;

TEST_CASE("BBOB2015-DiscusFunction", "") {
  for (const auto& numberOfDimensions : {2, 40}) {
    hop::bbob2015::DiscusFunction discusFunction(numberOfDimensions);

    arma::Mat<double> parameters;
    parameters.load(testDirectory.string() + "/data/optimisationProblem/blackBoxOptimisationBenchmark2015/parameters,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Col<double> translation;
    translation.load(testDirectory.string() + "/data/optimisationProblem/blackBoxOptimisationBenchmark2015/translation,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Mat<double> rotationR;
    rotationR.load(testDirectory.string() + "/data/optimisationProblem/blackBoxOptimisationBenchmark2015/rotationR,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Col<double> expected;
    expected.load(testDirectory.string() + "/data/optimisationProblem/blackBoxOptimisationBenchmark2015/expectedDiscusFunction,dim" + std::to_string(numberOfDimensions) +".mat");

    discusFunction.setObjectiveValueTranslation(0);
    discusFunction.setTranslation(translation);
    discusFunction.setRotationR(rotationR);

    for (std::size_t n = 0; n < parameters.n_cols; ++n) {
      CHECK(discusFunction.getObjectiveValue(parameters.col(n)) == Approx(expected.at(n)));
    }
  }
}