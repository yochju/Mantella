#pragma once

#include <optimisationProblem/benchmark/benchmarkProblem.hpp>

namespace hop {
  class RosenbrockFunctionRotated : public BenchmarkProblem {
    public:
      RosenbrockFunctionRotated(const unsigned int& numberOfDimensions);

    protected:
      const double _max;

      double getObjectiveValueImplementation(const arma::Col<double>& parameter) const override;
  };
}
