#include "mantella_bits/optimisationProblem.hpp"

// C++ standard library
#include <stdexcept>
#include <utility>

// Mantella
#include "mantella_bits/assert.hpp"
#include "mantella_bits/config.hpp"
#include "mantella_bits/mpi.hpp"

namespace mant {
  OptimisationProblem::OptimisationProblem(
      const arma::uword numberOfDimensions)
      : numberOfDimensions_(synchronise(numberOfDimensions)) {
    if (numberOfDimensions_ < 1) {
      throw std::domain_error("OptimisationProblem: The number of dimensions must be greater than 0.");
    }

    // Initialises the counters.
    reset();

    // Sets the default settings.
    setLowerBounds(arma::zeros<arma::vec>(numberOfDimensions_) - 10);
    setUpperBounds(arma::zeros<arma::vec>(numberOfDimensions_) + 10);

    // Sets all parameter and objective value space modifiers with their neutral element.
    setParameterPermutation(arma::regspace<arma::uvec>(0, numberOfDimensions_ - 1));
    setParameterScaling(arma::ones<arma::vec>(numberOfDimensions_));
    setParameterTranslation(arma::zeros<arma::vec>(numberOfDimensions_));
    if (numberOfDimensions_ > 1) {
      setParameterRotation(arma::eye<arma::mat>(numberOfDimensions_, numberOfDimensions_));
    }
    setMinimalParameterDistance(arma::zeros<arma::vec>(numberOfDimensions_));

    setObjectiveValueScaling(1.0);
    setObjectiveValueTranslation(0.0);
  }

  void OptimisationProblem::setObjectiveFunctions(
      const std::vector<std::pair<std::function<double(const arma::vec& parameter_)>, std::string>>& objectiveFunctions) {
    if (objectiveFunctions.size() == 0) {
      throw std::invalid_argument("OptimisationProblem.setObjectiveFunctions: At least one objective function must be defined.");
    }

    for (const auto& objectiveFunction : objectiveFunctions) {
      if (!static_cast<bool>(objectiveFunction.first)) {
        throw std::invalid_argument("OptimisationProblem.setObjectiveFunctions: All objective functions must be callable.");
      }
    }

    objectiveFunctions_ = objectiveFunctions;

    // Resets all counters and caches, as the problem could have changed.
    reset();
  }

  std::vector<std::pair<std::function<double(const arma::vec& parameter_)>, std::string>> OptimisationProblem::getObjectiveFunctions() const {
    return objectiveFunctions_;
  }

  double OptimisationProblem::getObjectiveValue(
      arma::vec parameter) {
    if (parameter.n_elem != numberOfDimensions_) {
      throw std::invalid_argument("OptimisationProblem.getObjectiveValue: The number of elements must be equal to the optimisation problem's number of dimensions.");
    } else if (objectiveFunctions_.size() == 0) {
      throw std::invalid_argument("OptimisationProblem.getObjectiveValue: At least one objective function must be defined.");
    }

    ++usedNumberOfEvaluations_;

    // Discretises the parameter
    const arma::uvec& elementsToDiscretise = arma::find(minimalParameterDistance_ > 0);
    parameter.elem(elementsToDiscretise) = arma::floor(parameter.elem(elementsToDiscretise) / minimalParameterDistance_.elem(elementsToDiscretise)) % minimalParameterDistance_.elem(elementsToDiscretise);
    if (numberOfDimensions_ > 1) {
      parameter = parameter.elem(parameterPermutation_);
    }
    parameter = (parameterScaling_ % parameter - parameterTranslation_);
    if (numberOfDimensions_ > 1) {
      parameter = parameterRotation_ * parameter;
    }

    double objectiveValue = 0.0;
    const auto n = cachedSamples_.find(parameter);
    if (n == cachedSamples_.cend()) {
      // Increases the number of distinct evaluations only, if we actually compute the value.
      ++usedNumberOfDistinctEvaluations_;

      for (const auto& objectiveFunction : objectiveFunctions_) {
        objectiveValue += objectiveFunction.first(parameter);
      }

      objectiveValue = objectiveValueScaling_ * objectiveValue + objectiveValueTranslation_;

      if (::mant::isCachingSamples) {
        cachedSamples_.insert({parameter, objectiveValue});
      }
    } else {
      objectiveValue = n->second;
    }

    return objectiveValue;
  }

  double OptimisationProblem::getObjectiveValueOfNormalisedParameter(
      const arma::vec& normalisedParameter) {
    if (normalisedParameter.n_elem != numberOfDimensions_) {
      throw std::invalid_argument("OptimisationProblem.getObjectiveValueOfNormalisedParameter: The number of elements must be equal to the optimisation problem's number of dimensions.");
    } else if (arma::any(lowerBounds_ > upperBounds_)) {
      throw std::logic_error("OptimisationProblem.getObjectiveValueOfNormalisedParameter: All upper bounds must be greater than their lower one.");
    }

    return getObjectiveValue(lowerBounds_ + normalisedParameter % (upperBounds_ - lowerBounds_));
  }

  void OptimisationProblem::setLowerBounds(
      const arma::vec& lowerBounds) {
    if (lowerBounds.n_elem != numberOfDimensions_) {
      throw std::invalid_argument("OptimisationProblem.setLowerBounds: The lower bounds' number of elements must be equal to the optimisation problem's number of dimensions.");
    }

    lowerBounds_ = synchronise(lowerBounds);
  }

  arma::vec OptimisationProblem::getLowerBounds() const {
    return lowerBounds_;
  }

  void OptimisationProblem::setUpperBounds(
      const arma::vec& upperBounds) {
    if (upperBounds.n_elem != numberOfDimensions_) {
      throw std::invalid_argument("OptimisationProblem.setUpperBounds: The upper bounds' number of elements must be equal to the optimisation problem's number of dimensions.");
    }

    upperBounds_ = synchronise(upperBounds);
  }

  arma::vec OptimisationProblem::getUpperBounds() const {
    return upperBounds_;
  }

  void OptimisationProblem::setParameterPermutation(
      const arma::uvec& parameterPermutation) {
    if (parameterPermutation.n_elem != numberOfDimensions_) {
      throw std::invalid_argument("OptimisationProblem.setParameterPermutation: The parameter permutation's number of elements must be equal to the optimisation problem's number of dimensions.");
    } else if (!isPermutationVector(parameterPermutation, numberOfDimensions_, numberOfDimensions_)) {
      throw std::invalid_argument("OptimisationProblem.setParameterPermutation: The parameter permutation must be an actual permutation on the optimisation problem.");
    }

    parameterPermutation_ = synchronise(parameterPermutation);

    // Resets all counters and caches, as the problem could have changed.
    reset();
  }

  arma::uvec OptimisationProblem::getParameterPermutation() const {
    return parameterPermutation_;
  }

  void OptimisationProblem::setParameterScaling(
      const arma::vec& parameterScaling) {
    if (parameterScaling.n_elem != numberOfDimensions_) {
      throw std::invalid_argument("OptimisationProblem.setParameterScaling: The parameter scaling's number of elements must be equal to the optimisation problem's number of dimensions.");
    }

    parameterScaling_ = synchronise(parameterScaling);

    // Resets all counters and caches, as the problem could have changed.
    reset();
  }

  arma::vec OptimisationProblem::getParameterScaling() const {
    return parameterScaling_;
  }

  void OptimisationProblem::setParameterTranslation(
      const arma::vec& parameterTranslation) {
    if (parameterTranslation.n_elem != numberOfDimensions_) {
      throw std::invalid_argument("OptimisationProblem.setParameterTranslation: The parameter translation's number of elements must be equal to the optimisation problem's number of dimensions.");
    }

    parameterTranslation_ = synchronise(parameterTranslation);

    // Resets all counters and caches, as the problem could have changed.
    reset();
  }

  arma::vec OptimisationProblem::getParameterTranslation() const {
    return parameterTranslation_;
  }

  void OptimisationProblem::setParameterRotation(
      const arma::mat& parameterRotation) {
    if (parameterRotation.n_rows != numberOfDimensions_) {
      throw std::invalid_argument("OptimisationProblem.setParameterRotation: The parameter rotation's number of rows must be equal to the optimisation problem's number of dimensions.");
    } else if (!isRotationMatrix(parameterRotation)) {
      throw std::invalid_argument("OptimisationProblem.setParameterRotation: The parameter rotation must be an actual rotation matrix.");
    }

    parameterRotation_ = synchronise(parameterRotation);

    // Resets all counters and caches, as the problem could have changed.
    reset();
  }

  arma::mat OptimisationProblem::getParameterRotation() const {
    return parameterRotation_;
  }

  void OptimisationProblem::setMinimalParameterDistance(
      const arma::vec& minimalParameterDistance) {
    if (minimalParameterDistance.n_elem != numberOfDimensions_) {
      throw std::invalid_argument("OptimisationProblem.setMinimalParameterDistance: The minimal parameter distance's number of elements must be equal to the optimisation problem's number of dimensions.");
    } else if (arma::any(minimalParameterDistance < 0)) {
      throw std::domain_error("OptimisationProblem.setMinimalParameterDistance: Each minimal parameter distance must be positive (including 0).");
    }

    minimalParameterDistance_ = synchronise(minimalParameterDistance);
  }

  arma::vec OptimisationProblem::getMinimalParameterDistance() const {
    return minimalParameterDistance_;
  }

  void OptimisationProblem::setObjectiveValueScaling(
      const double objectiveValueScaling) {
    objectiveValueScaling_ = synchronise(objectiveValueScaling);

    // Resets all counters and caches, as the problem could have changed.
    reset();
  }

  double OptimisationProblem::getObjectiveValueScaling() const {
    return objectiveValueScaling_;
  }

  void OptimisationProblem::setObjectiveValueTranslation(
      const double objectiveValueTranslation) {
    objectiveValueTranslation_ = synchronise(objectiveValueTranslation);

    // Resets all counters and caches, as the problem could have changed.
    reset();
  }

  double OptimisationProblem::getObjectiveValueTranslation() const {
    return objectiveValueTranslation_;
  }

  std::unordered_map<arma::vec, double, Hash, IsEqual> OptimisationProblem::getCachedSamples() const {
    return cachedSamples_;
  }

  arma::uword OptimisationProblem::getUsedNumberOfEvaluations() const {
    return usedNumberOfEvaluations_;
  }

  arma::uword OptimisationProblem::getUsedNumberOfDistinctEvaluations() const {
    return usedNumberOfDistinctEvaluations_;
  }

  void OptimisationProblem::reset() {
    usedNumberOfEvaluations_ = 0;
    usedNumberOfDistinctEvaluations_ = 0;

    cachedSamples_.clear();
  }
}
