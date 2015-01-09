#pragma once

// Mantella
#include <mantella_bits/optimisationProblem/blackBoxOptimisationBenchmark2010.hpp>

namespace mant {
  namespace bbob2010 {
    class EllipsoidalFunction : public BlackBoxOptimisationBenchmark2010 {
      public:
        using BlackBoxOptimisationBenchmark2010::BlackBoxOptimisationBenchmark2010;

        EllipsoidalFunction(const EllipsoidalFunction&) = delete;
        EllipsoidalFunction& operator=(const EllipsoidalFunction&) = delete;

        std::string to_string() const noexcept override;

      protected:
        const arma::Col<double> scaling_ = getScaling(1000000.0);

        double getObjectiveValueImplementation(
            const arma::Col<double>& parameter) const noexcept override;

        friend class cereal::access;

        template <typename Archive>
        void serialize(Archive& archive) noexcept {
          archive(cereal::make_nvp("BlackBoxOptimisationBenchmark2010", cereal::base_class<BlackBoxOptimisationBenchmark2010>(this)));
          archive(cereal::make_nvp("numberOfDimensions", numberOfDimensions_));
          archive(cereal::make_nvp("translation", translation_));
        }

        template <typename Archive>
        static void load_and_construct(
            Archive& archive,
            cereal::construct<EllipsoidalFunction>& construct) noexcept {
          unsigned int numberOfDimensions;
          archive(cereal::make_nvp("numberOfDimensions", numberOfDimensions));
          construct(numberOfDimensions);

          archive(cereal::make_nvp("BlackBoxOptimisationBenchmark2010", cereal::base_class<BlackBoxOptimisationBenchmark2010>(construct.ptr())));
          archive(cereal::make_nvp("translation", construct->translation_));
        }
    };
  }
}