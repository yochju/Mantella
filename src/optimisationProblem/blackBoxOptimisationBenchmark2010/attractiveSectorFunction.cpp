#include <hop_bits/optimisationProblem/blackBoxOptimisationBenchmark2010/attractiveSectorFunction.hpp>

// Cereal
#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>

namespace hop {
  namespace bbob2010 {
    double AttractiveSectorFunction::getObjectiveValueImplementation(
        const arma::Col<double>& parameter) const noexcept {
      arma::Col<double> z = rotationQ_ * (delta_ % (rotationR_ * (parameter - translation_)));
      z.elem(arma::find(z % translation_ > 0.0)) *= 100.0;

      return std::pow(getOscillationTransformation(std::pow(norm(z), 2.0)), 0.9);
    }

    std::string AttractiveSectorFunction::to_string() const noexcept {
      return "AttractiveSectorFunction";
    }
  }
}

CEREAL_REGISTER_TYPE(hop::bbob2010::AttractiveSectorFunction)