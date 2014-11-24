#include <hop_bits/optimisationProblem/benchmark/blackBoxOptimisationBenchmark2013/bentCigarFunction.hpp>

// Cereal
#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>

namespace hop {
  namespace bbob2013 {
    double BentCigarFunction::getObjectiveValueImplementation(
        const arma::Col<double>& parameter) const noexcept {
      arma::Col<double> z = arma::square(rotationR_ * getAsymmetricTransformation(0.5, rotationR_ * (parameter - translation_)));
      return z.at(0) + 1000000.0 * arma::accu(z.subvec(1, z.n_elem - 1));
    }

    std::string BentCigarFunction::to_string() const noexcept {
      return "BentCigarFunction";
    }
  }
}

CEREAL_REGISTER_TYPE(hop::bbob2013::BentCigarFunction)
