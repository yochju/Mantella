#pragma once

// C++ Standard Library
#include <cmath>

// HOP
#include <hop_bits/optimisationProblem/blackBoxOptimisationBenchmark2015.hpp>

namespace hop {
  namespace bbob2015 {
    class RastriginFunction : public BlackBoxOptimisationBenchmark2015 {
      public:
        using BlackBoxOptimisationBenchmark2015::BlackBoxOptimisationBenchmark2015;

        RastriginFunction(const RastriginFunction&) = delete;
        RastriginFunction& operator=(const RastriginFunction&) = delete;

        std::string to_string() const noexcept override;

      protected:
        const arma::Col<double> delta_ = getScaling(std::sqrt(10.0));

        double getObjectiveValueImplementation(
            const arma::Col<double>& parameter) const noexcept override;

        friend class cereal::access;

        template <typename Archive>
        void serialize(
            Archive& archive) noexcept {
          archive(cereal::make_nvp("BlackBoxOptimisationBenchmark2015", cereal::base_class<BlackBoxOptimisationBenchmark2015>(this)));
          archive(cereal::make_nvp("numberOfDimensions", numberOfDimensions_));
          archive(cereal::make_nvp("translation", translation_));
        }

        template <typename Archive>
        static void load_and_construct(
            Archive& archive,
            cereal::construct<RastriginFunction>& construct) noexcept {
          unsigned int numberOfDimensions;
          archive(cereal::make_nvp("numberOfDimensions", numberOfDimensions));
          construct(numberOfDimensions);

          archive(cereal::make_nvp("BlackBoxOptimisationBenchmark2015", cereal::base_class<BlackBoxOptimisationBenchmark2015>(construct.ptr())));
          archive(cereal::make_nvp("translation", construct->translation_));
        }
    };
  }
}