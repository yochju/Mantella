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

TEST_CASE("SphereFunction", "") {
  for (const auto& numberOfDimensions : {2, 40}) {
    hop::bbob2013::SphereFunction sphereFunction(numberOfDimensions);

    arma::Mat<double> parameters;
    parameters.load(testDirectory.string() + "/data/optimisationProblem/benchmark/blackBoxOptimisationBenchmark2013/parameters,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Col<double> translation;
    translation.load(testDirectory.string() + "/data/optimisationProblem/benchmark/blackBoxOptimisationBenchmark2013/translation,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Col<double> expected;
    expected.load(testDirectory.string() + "/data/optimisationProblem/benchmark/blackBoxOptimisationBenchmark2013/expectedSphereFunction,dim" + std::to_string(numberOfDimensions) +".mat");

    sphereFunction.setObjectiveValueTranslation(0);
    sphereFunction.setTranslation(translation);

    for (std::size_t n = 0; n < parameters.n_cols; ++n) {
      CHECK(sphereFunction.getObjectiveValue(parameters.col(n)) == Approx(expected.at(n)));
    }
  }
}