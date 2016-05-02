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

    //0: parameter wird in MJD angegeben (nicht JD und nicht MJD2000)
    arma::Col<double> GTOC1::objectiveFunction(const arma::Col<double>& parameter) {
      double totalTimePassed = 0.0;
      
      /*
      std::cout << "t as input: ";
      for(double el : parameter){
        std::cout << el << " ";
      }
      std::cout << std::endl;
      */
      
      verify(parameter.n_elem == orbitalTargetSequence_.size(), "gtoc1.objectiveFunction: Number of parameter and sequence length must be equal.");

      //1: Ab hier werden erstmal alle Positionen und Beschleunigungen der Planeten passend zur Zeit berechnet und zwischengespeichert
      std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> orbitBodyPositionsAndVelocities;

      for (arma::uword i = 0; i < parameter.n_elem; i++) {
        totalTimePassed += parameter(i);
        orbitBodyPositionsAndVelocities.push_back(orbitOnPosition(totalTimePassed, orbitalTargetSequence_.at(i)));
      }
      /*
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
      */
      
      double deltaVelocities = 0;

      //2: Initiale Flugstrecke
      //2.1: Lambert vom Startplaneten zum ersten Ziel. Hier werden die Positionen der vorher berechneten Planeten übergeben.
      std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> lambertVelocities = lambert(orbitBodyPositionsAndVelocities.at(0).first, orbitBodyPositionsAndVelocities.at(1).first, parameter(1) * 86400.0); // in seconds
      
      //2.2: Ab hier werden die besten Beschleunigungen bestimmt (Min. Differenz der jeweiligen Längen)
      std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> bestVelocities = lambertVelocities.at(0);
      double currentDV = std::abs(arma::norm(bestVelocities.first) - arma::norm(bestVelocities.second));
      
      /*std::cout << "\t" << currentDV << " - (" << bestVelocities.first(0) << ", " << bestVelocities.first(1) << ", " << bestVelocities.first(2) << ")" << std::endl;
      std::cout << "\ttest " << std::abs(arma::norm(lambertVelocities.at(1).first) - arma::norm(lambertVelocities.at(1).second)) << " - (" << lambertVelocities.at(1).first(0) << ", " << lambertVelocities.at(1).first(1) << ", " << lambertVelocities.at(1).first(2) << ")" << std::endl;
      */
      
      //2.3: Prüfung aller möglichen Lösungen
      for (arma::uword iter = 1; iter < lambertVelocities.size(); iter++) {
        if (currentDV > std::abs(arma::norm(lambertVelocities.at(iter).first) - arma::norm(lambertVelocities.at(iter).second))) {
          currentDV = std::abs(arma::norm(lambertVelocities.at(iter).first) - arma::norm(lambertVelocities.at(iter).second));
          bestVelocities = lambertVelocities.at(iter);
        }
        //std::cout << "\t" << currentDV << " - (" << bestVelocities.first(0) << ", " << bestVelocities.first(1) << ", " << bestVelocities.first(2) << ")" << std::endl;
      }
      
      /*
      std::cout << "lambert solution1: (" << bestVelocities.first(0) << ", " << bestVelocities.first(1) << ", " << bestVelocities.first(2) << ") ~ " << "(" << bestVelocities.second(0) << ", " << bestVelocities.second(1) << ", " << bestVelocities.second(2) << ") " << std::endl; 
      */
      
      //2.4: Bestimmung der zusätzlichen Startbeschleunigung (penalty aus gtopToolbox)
      // Launcher Constraint
      if(arma::norm(bestVelocities.first - orbitBodyPositionsAndVelocities.at(0).first) > rocketDVlaunch_){
        deltaVelocities += arma::norm(bestVelocities.first) - arma::norm(bestVelocities.second) - rocketDVlaunch_;
      }
      
      //Zum Vergleich mit alten Beschleunigungen
      std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> pastBestVelocities;

      //3: Folgende Flugstrecken bis zum Zielplaneten. Ab hier Iteration über alle möglichen Flugstrecken zwischen den Planeten der Sequenz inklusive gravityAssist
      for (arma::uword i = 2; i < parameter.n_elem; i++) {
        pastBestVelocities = bestVelocities;
        //3.1: Lambert zwischen zwei Planeten. Hier werden die Positionen der vorher berechneten Planeten übergeben.
        lambertVelocities = lambert(orbitBodyPositionsAndVelocities.at(i - 1).first, orbitBodyPositionsAndVelocities.at(i).first, parameter(i) * 86400.0);

        //3.2: Ab hier werden die besten Beschleunigungen bestimmt (Min. Differenz der jeweiligen Längen)
        bestVelocities = lambertVelocities.at(0);
        
        currentDV = std::abs(arma::norm(bestVelocities.first) - arma::norm(bestVelocities.second));
        
        //std::cout << "\t" << currentDV << " - (" << bestVelocities.first(0) << ", " << bestVelocities.first(1) << ", " << bestVelocities.first(2) << ")" << std::endl;
        
        //3.3: Prüfung aller möglichen Lösungen
        for (arma::uword iter = 1; iter < lambertVelocities.size(); iter++) {
          if (currentDV > std::abs(arma::norm(lambertVelocities.at(iter).first) - arma::norm(lambertVelocities.at(iter).second))) {
            currentDV = std::abs(arma::norm(lambertVelocities.at(iter).first) - arma::norm(lambertVelocities.at(iter).second));
            bestVelocities = lambertVelocities.at(iter);
          }
          //std::cout << "\t" << currentDV << " - (" << bestVelocities.first(0) << ", " << bestVelocities.first(1) << ", " << bestVelocities.first(2) << ")" << std::endl;
        }
        
        //std::cout << "lambert solution" << i << ": (" << bestVelocities.first(0) << ", " << bestVelocities.first(1) << ", " << bestVelocities.first(2) << ") ~ " << "(" << bestVelocities.second(0) << ", " << bestVelocities.second(1) << ", " << bestVelocities.second(2) << ") " << std::endl; 

        //deltaVelocities += arma::norm(bestVelocities.first - bestVelocities.second);

        //3.4 Vorherige und aktuelle Beschleunigung wird übergeben und der Swing-by berechnet
        std::pair<double, double> dvAndRp = gravityAssist(pastBestVelocities.second - orbitBodyPositionsAndVelocities.at(i - 1).second, bestVelocities.first - orbitBodyPositionsAndVelocities.at(i - 1).second);
        deltaVelocities += dvAndRp.first;
        
        //printf("Vin, Vout, DV, rp - %f, %f, %f, %f\n", arma::norm(pastBestVelocities.second), arma::norm(bestVelocities.first), dvAndRp.first, dvAndRp.second);

        
        //---Fehlende penalty (?): Bei dem Codeabschnitt in der gtopToolbox bin ich mir nicht so sicher, was es mathematisch bedeuten mag...
        /* 
        	for (i_count = 0; i_count < 3; i_count++)
          Dum_Vec[i_count] = v[n-1][i_count] - V_Lamb[1][1][i_count];
          DVrel = norm2(Dum_Vec);      
          DVarr = DVrel;     
        */
      }

      //---Fehlende penalty (?): Bei dem Codeabschnitt in der gtopToolbox bin ich mir nicht so sicher, was es mathematisch bedeuten mag...
      /* ??????????? TODO
        for (i_count = 0;i_count < n-2; i_count++)
          if (rp[i_count] < penalty[sequence[i_count+1]])
            DVtot += penalty_coeffs[sequence[i_count+1]]*fabs(rp[i_count] - penalty[sequence[i_count+1]]);
      */

      //std::cout << "tmp output: " << rocketSpecificImpulse_ * 9.80665 << " " << -deltaVelocities / (rocketSpecificImpulse_ * 9.80665) << " " << std::exp(-deltaVelocities / (rocketSpecificImpulse_ * 9.80665)) << std::endl;
      double finalMass = rocketMass_ * std::exp(-deltaVelocities / (rocketSpecificImpulse_ * 9.80665));
      //std::cout << "finalMass: " << finalMass << std::endl;
      arma::Col<double>::fixed<3> relativeVelocityToLastTarget = orbitBodyPositionsAndVelocities.back().second - bestVelocities.second;
      //std::cout << "relativeVelocityToLastTarget: " << relativeVelocityToLastTarget << std::endl;
      
      //std::cout << "stuff2: " << arma::dot(relativeVelocityToLastTarget, bestVelocities.second) << " " << -finalMass * std::abs(arma::dot(relativeVelocityToLastTarget, bestVelocities.second)) << std::endl;
      return {-finalMass * std::abs(arma::dot(relativeVelocityToLastTarget, bestVelocities.second))};
    }
  }
}
