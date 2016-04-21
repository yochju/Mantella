#include "mantella_bits/orbitalProblem/gtoc1.hpp"

// Mantella
#include "mantella_bits/assert.hpp"
#include "mantella_bits/orbitalMechanics.hpp"

namespace mant {
  namespace itd {

    GTOC1::GTOC1(const std::vector<arma::Col<double>::fixed<7>>& orbitalTargetSequence) {
      verify(orbitalTargetSequence.size() >= 2, "gtoc1: Sequence must have at least size of two.");
      orbitalTargetSequence_ = orbitalTargetSequence;
    }

    arma::Col<double> GTOC1::problemFunction(const arma::Col<double>& parameter) {
      double totalTimePassed = 0.0;
      
      std::cout << "t as input: ";
      for(double el : parameter){
        std::cout << el << " ";
      }
      std::cout << std::endl;

      verify(parameter.n_elem == orbitalTargetSequence_.size(), "gtoc1.problemFunction: Number of parameter and sequence length must be equal.");

      std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> orbitBodyPositionsAndVelocities;

      for (int i = 0; i < parameter.n_elem; i++) {
        totalTimePassed += parameter(i);
        orbitBodyPositionsAndVelocities.push_back(orbitOnPosition(totalTimePassed, orbitalTargetSequence_.at(i)));
      }
      std::cout << "All positions: ";
      for(auto pair : orbitBodyPositionsAndVelocities){  
        arma::Col<double>::fixed<3> elem = pair.first;
        std::cout << "(" << elem(0) << ", " << elem(1) << ", " << elem(2) << ") ";
      }
      std::cout << std::endl;
    
      std::cout << "All velocities: ";
      for(auto pair : orbitBodyPositionsAndVelocities){  
        arma::Col<double>::fixed<3> elem = pair.second;
        std::cout << "(" << elem(0) << ", " << elem(1) << ", " << elem(2) << ") ";
      }
      std::cout << std::endl;

      double deltaVelocities = 0;

      std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> lambertVelocities = lambert(orbitBodyPositionsAndVelocities.at(0).first, orbitBodyPositionsAndVelocities.at(1).first, parameter(1) * 86400.0);
      
      std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> bestVelocities;
      double currentDV = std::numeric_limits<double>::max();

      for (auto pair : lambertVelocities) {
        if (currentDV > arma::norm(pair.first) - arma::norm(pair.second)) {
          bestVelocities = pair;
        }
      }
      
      std::cout << "lambert solution: (" << bestVelocities.first(0) << ", " << bestVelocities.first(1) << ", " << bestVelocities.first(2) << ") ~ " << "(" << bestVelocities.second(0) << ", " << bestVelocities.second(1) << ", " << bestVelocities.second(2) << ") " << std::endl; 
      
      // Launcher Constraint
      if(arma::norm(bestVelocities.first - orbitBodyPositionsAndVelocities.at(0).first) > rocketDVlaunch_){
        deltaVelocities += arma::norm(bestVelocities.first) - arma::norm(bestVelocities.second) - rocketDVlaunch_;
      }
      
      std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> pastBestVelocities;

      for (int i = 2; i < parameter.n_elem; i++) {
        pastBestVelocities = bestVelocities;

        lambertVelocities = lambert(orbitBodyPositionsAndVelocities.at(i - 1).first, orbitBodyPositionsAndVelocities.at(i).first, parameter(i) * 86400.0);

        currentDV = std::numeric_limits<double>::max();

        for (auto pair : lambertVelocities) {
          if (currentDV > arma::norm(pair.first) - arma::norm(pair.second)) {
            bestVelocities = pair;
          }
        }
        
        std::cout << "lambert solution: (" << bestVelocities.first(0) << ", " << bestVelocities.first(1) << ", " << bestVelocities.first(2) << ") ~ " << "(" << bestVelocities.second(0) << ", " << bestVelocities.second(1) << ", " << bestVelocities.second(2) << ") " << std::endl; 

        //deltaVelocities += arma::norm(bestVelocities.first - bestVelocities.second);

        std::pair<double, double> dvAndRp = gravityAssist(pastBestVelocities.second - orbitBodyPositionsAndVelocities.at(i - 1).second, bestVelocities.first - orbitBodyPositionsAndVelocities.at(i - 1).second);
        deltaVelocities += dvAndRp.first;
        
        printf("Vin, Vout, DV, rp - %f, %f, %f, %f\n", arma::norm(pastBestVelocities.second), arma::norm(bestVelocities.first), dvAndRp.first, dvAndRp.second);

        /* ????????
        	for (i_count = 0; i_count < 3; i_count++)
          Dum_Vec[i_count] = v[n-1][i_count] - V_Lamb[1][1][i_count];
          DVrel = norm2(Dum_Vec);      
          DVarr = DVrel;     
        */
      }

      // penalties
      /* ??????????? TODO
        for (i_count = 0;i_count < n-2; i_count++)
          if (rp[i_count] < penalty[sequence[i_count+1]])
            DVtot += penalty_coeffs[sequence[i_count+1]]*fabs(rp[i_count] - penalty[sequence[i_count+1]]);
      */

      // Launcher Constraint
      /*
      if (DV[0] > DVlaunch){
        DVtot += (DV[0] - rocketDVlaunch_);
      }
      */
      std::cout << "stuff: " << rocketSpecificImpulse_ * 9.80665 << " " << -deltaVelocities / (rocketSpecificImpulse_ * 9.80665) << " " << std::exp(-deltaVelocities / (rocketSpecificImpulse_ * 9.80665)) << std::endl;
      double finalMass = rocketMass_ * std::exp(-deltaVelocities / (rocketSpecificImpulse_ * 9.80665));
      //std::cout << "finalMass: " << finalMass << std::endl;
      arma::Col<double>::fixed<3> relativeVelocityToLastTarget = orbitBodyPositionsAndVelocities.back().second - bestVelocities.second;
      //std::cout << "relativeVelocityToLastTarget: " << relativeVelocityToLastTarget << std::endl;
      
      //std::cout << "stuff2: " << arma::dot(relativeVelocityToLastTarget, bestVelocities.second) << " " << -finalMass * std::abs(arma::dot(relativeVelocityToLastTarget, bestVelocities.second)) << std::endl;
      return {-finalMass * std::abs(arma::dot(relativeVelocityToLastTarget, bestVelocities.second))};
    }
  }
}
