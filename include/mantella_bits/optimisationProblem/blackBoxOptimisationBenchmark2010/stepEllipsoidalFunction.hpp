#pragma once

// C++ Standard Library
#include <cmath>

// Mantella
#include <mantella_bits/optimisationProblem/blackBoxOptimisationBenchmark2010.hpp>

namespace mant {
  namespace bbob2010 {
    class StepEllipsoidalFunction : public BlackBoxOptimisationBenchmark2010 {
      public:
        using BlackBoxOptimisationBenchmark2010::BlackBoxOptimisationBenchmark2010;

        StepEllipsoidalFunction(const StepEllipsoidalFunction&) = delete;
        StepEllipsoidalFunction& operator=(const StepEllipsoidalFunction&) = delete;

        std::string to_string() const noexcept override;

      protected:
        const arma::Col<double> scaling_ = getScaling(100.0);
        const arma::Col<double> delta_ = getScaling(std::sqrt(10.0));

        double getObjectiveValueImplementation(
            const arma::Col<double>& parameter) const noexcept override;

        friend class cereal::access;

        template <typename Archive>
        void serialize(
            Archive& archive) noexcept {
          archive(cereal::make_nvp("BlackBoxOptimisationBenchmark2010", cereal::base_class<BlackBoxOptimisationBenchmark2010>(this)));
          archive(cereal::make_nvp("numberOfDimensions", numberOfDimensions_));
          archive(cereal::make_nvp("translation", translation_));
          archive(cereal::make_nvp("rotationR", rotationR_));
          archive(cereal::make_nvp("rotationQ", rotationQ_));
        }

        template <typename Archive>
        static void load_and_construct(
            Archive& archive,
            cereal::construct<StepEllipsoidalFunction>& construct) noexcept {
          unsigned int numberOfDimensions;
          archive(cereal::make_nvp("numberOfDimensions", numberOfDimensions));
          construct(numberOfDimensions);

          archive(cereal::make_nvp("BlackBoxOptimisationBenchmark2010", cereal::base_class<BlackBoxOptimisationBenchmark2010>(construct.ptr())));
          archive(cereal::make_nvp("translation", construct->translation_));
          archive(cereal::make_nvp("rotationR", construct->rotationR_));
          archive(cereal::make_nvp("rotationQ", construct->rotationQ_));
        }
    };
  }
}