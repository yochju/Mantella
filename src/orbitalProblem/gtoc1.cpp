// Mantella
#include "mantella_bits/assert.hpp"
#include "mantella_bits/orbitalMechanics.hpp"
#include "mantella_bits/orbitalProblem/gtoc1.hpp"

namespace mant {
  namespace itd {
    
    GTOC1::GTOC1(std::vector<arma::Col<double>::fixed<6>> orbitalTargetSequence) {
        
        verify(orbitalTargetSequence.size() >= 2, "gtoc1: Sequence must have at least size of two.");
        orbitalTargetSequence_ = orbitalTargetSequence;
      }
    
    arma::Col<double> GTOC1::problemFunction(arma::Col<double> parameter){
      
      double totalTimePassed = 0.0;
      
      verify(parameter.n_cols == orbitalTargetSequence_.size(), "gtoc1.problemFunction: Number of parameter and sequence length must be equal.");
      
      std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> orbitBodyPositionsAndVelocities;
      
      for (int i = 0; i < parameter.n_cols; i++) {
        totalTimePassed += parameter(i);
        orbitBodyPositionsAndVelocities.push_back(orbitOnPosition(totalTimePassed, orbitalTargetSequence_.at(i)));
      }
      
      double deltaVelocities = 0;
      
      std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> lambertVelocities = lambert(orbitBodyPositionsAndVelocities.at(0).first, orbitBodyPositionsAndVelocities.at(1).first, parameter(1));
      
      std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> bestVelocities;
      double currentDV = std::numeric_limits<double>::max();
        
      for (auto pair : lambertVelocities) {
                 
        if(currentDV > arma::norm(pair.first - pair.second)){
          bestVelocities = pair;
        }
      }
      
      deltaVelocities += arma::norm(bestVelocities.first - bestVelocities.second) - rocketDVlaunch_;
      
      std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> pastBestVelocities;
      
      for (int i = 2; i < parameter.n_cols; i++) {
              
        pastBestVelocities = bestVelocities;
        
        lambertVelocities = lambert(orbitBodyPositionsAndVelocities.at(i-1).first, orbitBodyPositionsAndVelocities.at(i).first, parameter(i));
        
        currentDV = std::numeric_limits<double>::max();
        
        for (auto pair : lambertVelocities) {
                   
          if(currentDV > arma::norm(pair.first - pair.second)){
            bestVelocities = pair;
          }
        }
        
        deltaVelocities += arma::norm(bestVelocities.first - bestVelocities.second);
        
        std::pair<double, double> dvAndRp = gravityAssist(pastBestVelocities.second, bestVelocities.first);
        deltaVelocities += dvAndRp.first;
        
        /* ????????
        	for (i_count = 0; i_count < 3; i_count++)
          Dum_Vec[i_count] = v[n-1][i_count] - V_Lamb[1][1][i_count];
          DVrel = norm2(Dum_Vec);      
          DVarr = DVrel;     
        */ 
      }
      
      // penalties
      /* ???????????
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
      
      double finalMass = rocketMass_ * std::exp(-deltaVelocities / (rocketSpecificImpulse_ * 0.00980665));
      arma::Col<double>::fixed<3> relativeVelocityToLastTarget = orbitBodyPositionsAndVelocities.back().second - bestVelocities.second;
 
      return {-finalMass * std::abs(arma::dot(relativeVelocityToLastTarget, bestVelocities.second))};
    }
  }
}
