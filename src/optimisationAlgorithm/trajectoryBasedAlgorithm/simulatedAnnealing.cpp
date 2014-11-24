#include <hop_bits/optimisationAlgorithm/trajectoryBasedAlgorithm/simulatedAnnealing.hpp>

// C++ Standard Library
#include <cmath>

// HOP
#include <hop_bits/helper/rng.hpp>

namespace hop {
  SimulatedAnnealing::SimulatedAnnealing(
      const std::shared_ptr<OptimisationProblem> optimisationProblem) noexcept
    : TrajectoryBasedAlgorithm(optimisationProblem),
      candidateObjectiveValue_(std::numeric_limits<double>::infinity()),
      candidateSoftConstraintValue_(std::numeric_limits<double>::infinity()) {
    setMaximalStepSize((optimisationProblem->getUpperBounds() - optimisationProblem->getLowerBounds()) / 10);
  }

  void SimulatedAnnealing::optimiseImplementation() noexcept {
    ++numberOfIterations_;

    bestParameter_ = initialParameter_;
    bestSoftConstraintValue_ = optimisationProblem_->getSoftConstraintsValue(initialParameter_);
    bestObjectiveValue_ = optimisationProblem_->getObjectiveValue(initialParameter_);

    state_ = bestParameter_;
    while(!isFinished() && !isTerminated()) {
      ++numberOfIterations_;

      candidateParameter_ = state_ + maximalStepSize_ % getVelocity();
      candidateSoftConstraintValue_ = optimisationProblem_->getSoftConstraintsValue(candidateParameter_);
      candidateObjectiveValue_ = optimisationProblem_->getObjectiveValue(candidateParameter_);

      if(candidateSoftConstraintValue_ < bestSoftConstraintValue_ || candidateSoftConstraintValue_ == bestSoftConstraintValue_ && candidateObjectiveValue_ < bestObjectiveValue_) {
        state_ = candidateParameter_;

        bestParameter_ = candidateParameter_;
        bestSoftConstraintValue_ = candidateSoftConstraintValue_;
        bestObjectiveValue_ = candidateObjectiveValue_;
      } else if(isAcceptableState()) {
        state_ = candidateParameter_;
      }
    }
  }

  bool SimulatedAnnealing::isAcceptableState() noexcept {
    return std::exp((bestObjectiveValue_ - candidateObjectiveValue_) / (numberOfIterations_ / maximalNumberOfIterations_)) < std::uniform_real_distribution<double>(0, 1)(Rng::generator);
  }

  arma::Col<double> SimulatedAnnealing::getVelocity() noexcept {
    return arma::normalise(arma::randn<arma::Col<double>>(optimisationProblem_->getNumberOfDimensions())) * std::uniform_real_distribution<double>(0, 1)(Rng::generator);
  }

  void SimulatedAnnealing::setMaximalStepSize(
      const arma::Col<double>& maximalStepSize) {
    if(maximalStepSize.n_rows != optimisationProblem_->getNumberOfDimensions()) {
      throw std::logic_error("The dimension of the maximal step size (" + std::to_string(maximalStepSize.n_elem) + ") must match the dimension of the optimisation problem (" + std::to_string(optimisationProblem_->getNumberOfDimensions()) + ").");
    }
    maximalStepSize_ = maximalStepSize;
  }

  std::string SimulatedAnnealing::to_string() const noexcept {
    return "SimulatedAnnealing";
  }
}
