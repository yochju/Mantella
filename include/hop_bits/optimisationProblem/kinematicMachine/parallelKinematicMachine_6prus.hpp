#pragma once

// HOP
#include <hop_bits/optimisationProblem/kinematicMachine.hpp>

namespace hop {
  namespace robotic {
    class ParallelKinematicMachine_6PRUS : public KinematicMachine {
      public:
        explicit ParallelKinematicMachine_6PRUS() noexcept;
        explicit ParallelKinematicMachine_6PRUS(
            const arma::Mat<double>::fixed<2, 6>& baseJointsAngles,
            const arma::Mat<double>::fixed<3, 6>& relativeEndEffectorJoints,
            const arma::Mat<double>::fixed<2, 6>& linkLengths,
            const arma::Mat<double>::fixed<3, 6>& redundantJointStarts,
            const arma::Mat<double>::fixed<3, 6>& redundantJointEnds) noexcept;

        // Copy constructors are not used in this library and deleted to avoid unintended/any usage.
        ParallelKinematicMachine_6PRUS(const ParallelKinematicMachine_6PRUS&) = delete;
        ParallelKinematicMachine_6PRUS& operator=(const ParallelKinematicMachine_6PRUS&) = delete;

        std::vector<arma::Mat<double>> getModelCharacterisation(
            const arma::Col<double>& endEffectorPose,
            const arma::Mat<double>& redundantJointActuations) const noexcept override;

        arma::Mat<double> getActuation(
            const arma::Col<double>& endEffectorPose,
            const arma::Mat<double>& redundantJointActuations) const noexcept override;

        double getPositionError(
            const arma::Col<double>& endEffectorPose,
            const arma::Mat<double>& redundantJointActuations) const noexcept override;

      protected:
        arma::Mat<double>::fixed<3, 6> endEffectorJointsRelative_;
        arma::Mat<double>::fixed<3, 6> linkLengths_;

        arma::Cube<double>::fixed<3, 3, 6> baseJointsRotation_;
        arma::Mat<double>::fixed<3, 6> baseJointsNormal_;

        arma::Mat<double>::fixed<3, 6> redundantJointStarts_;
        arma::Mat<double>::fixed<3, 6> redundantJointEnds_;
        arma::Mat<double>::fixed<3, 6> redundantJointsStartToEnd_;
        arma::Col<arma::uword> redundantJointIndicies_;
        arma::Mat<double> redundantJointAngles_;
    };
  }
}