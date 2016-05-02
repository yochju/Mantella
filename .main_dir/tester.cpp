#include <iostream>
#include <mantella>
#include <keplerian_toolbox/keplerian_toolbox.h>
#include <gtop>

#include <fstream>

//#include "planetData.h"
//#include "gridAnalysis.h"

int main() {
  
  // brent
  /*
  std::cout << "finite: " << std::isfinite(arma::datum::inf) << std::endl;
  
  double output = mant::brent([](double param){ 
      if(std::abs(param + 1.5 * arma::datum::pi) < 1e-3) {
        return arma::datum::inf;        
      } else if(std::abs(param + 2.5 * arma::datum::pi) < 1e-3) {
        return arma::datum::inf;
      }
      return std::tan(param); 
    }, -2.5 * arma::datum::pi, -1.5 * arma::datum::pi, 100, 1e-3);
  
  std::cout << "brent result: " << output << std::endl;
  */ 
  
  // GTOC1
  std::cout << "--Mantella GTOC1" << std::endl;
  arma::Col<double>::fixed<7> venus = {7.233268496749391E-01, 6.755697267164094E-03, 3.394589632336535E+00, 5.518541455452200E+01, 7.667837563371675E+01, 4.931425178852966E+01, 51544.0};  
  arma::Col<double>::fixed<7> earth = {1.000371833989169E+00, 1.704239716781501E-02, 2.669113820737183E-04, 2.977668064579176E+02, 1.639752443600624E+02, 3.581891404220149E+02, 51544.0};  
  //arma::Col<double>::fixed<7> uranus = {1.922994520785802E+01, 4.439340361752395E-02, 7.723573015712479E-01, 9.661124460900974E+01, 7.395999806487224E+01, 1.429075754495578E+02, 51544.0};  
  arma::Col<double>::fixed<7> astroid2001TW229 = {2.5897261, 0.2734625, 6.40734, 128.34711, 264.78691, 320.479555, 53600.0};  

  std::vector<arma::Col<double>::fixed<7>> seq;

  seq.push_back(earth);
  seq.push_back(venus);
  seq.push_back(earth);
  seq.push_back(venus);
  seq.push_back(astroid2001TW229);
  
  mant::itd::GTOC1 gtoc1Object(seq);
  
  //MJD Format
  arma::Col<double> param = {58353.0, 350.5, 1100.3, 330.5, 200.0};
  arma::Col<double> result = gtoc1Object.problemFunction(param);
  
  std::cout << "gtoc1 result: " << result(0) << std::endl << std::endl;
 

  std::cout << "--GTOP GTOC1" << std::endl;

  //23/08/2018 00:00:00.000
  //MJD2000 Format
  std::vector<double> variables = {58353.0 - 51544.0, 350.5, 1100.3, 330.5, 200.0};
  std::vector<double> rp = {6352, 6779, 6352, 6779, 600001, 70001};
  
  std::cout << "gtoc1 result: " << gtoc1(variables, rp) << std::endl;
  
  /*
  // Orbit on position
  std::cout << "---GTOC for Venus" << std::endl;
  double pos[3], vel[3];
  Planet_Ephemerides_Analytical(1234.0, 2, pos, vel);
  std::cout << "pos = (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
  std::cout << "vel = (" << vel[0] << ", " << vel[1] << ", " << vel[2] << ")" << std::endl << std::endl;

  
  std::cout << "---PYKEP Mat for venus" << std::endl;
  kep_toolbox::array3D r,v;
  kep_toolbox::planet::jpl_lp planet("venus");

  planet.eph(0.0, r, v);
  std::cout << "pos = (" << r.at(0) << ", " << r.at(1) << ", " << r.at(2) << ")" << std::endl;
  std::cout << "vel = (" << v.at(0) << ", " << v.at(1) << ", " << v.at(2) << ")" << std::endl << std::endl;


  
  std::cout << "---MANTELLA Vec for venus" << std::endl;
  arma::Col<double>::fixed<7> keplerianElementsVenusVec = {7.233268496749391E-01, 6.755697267164094E-03, 3.394589632336535E+00, 5.518541455452200E+01, 7.667837563371675E+01, 4.931425178852966E+01, 51544.0};
  arma::Col<double>::fixed<7> keplerianElementsUranusVec = {1.922994520785802E+01, 4.439340361752395E-02, 7.723573015712479E-01, 9.661124460900974E+01, 7.395999806487224E+01, 1.429075754495578E+02, 51544.0};
  //arma::Mat<double>::fixed<2,6> keplerianElements =  {{arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf},
  //                                                    {arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf}};
  std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> posVelPair = mant::itd::orbitOnPosition(51544.0, keplerianElementsVenusVec);

  std::cout << "pos = (" << posVelPair.first(0) << ", " << posVelPair.first(1) << ", " << posVelPair.first(2) << ")" << std::endl;
  std::cout << "vel = (" << posVelPair.second(0) << ", " << posVelPair.second(1) << ", " << posVelPair.second(2) << ")" << std::endl << std::endl;
  
  
  std::cout << "---PYKEP Vec for Astr" << std::endl;
  //kep_toolbox::array3D r,v;
  kep_toolbox::planet::gtoc2 asteroid(0);

  asteroid.eph(1234.0, r, v);
  //std::cout << asteroid.human_readable_extra() << std::endl;
  std::cout << "pos = (" << r.at(0) << ", " << r.at(1) << ", " << r.at(2) << ")" << std::endl;
  std::cout << "vel = (" << v.at(0) << ", " << v.at(1) << ", " << v.at(2) << ")" << std::endl << std::endl;

  
  std::cout << "---MANTELLA Vec for Astr" << std::endl;
  arma::Col<double>::fixed<7> keplerianElementsVec =  {3.9501468,0.2391642,6.87574,16.88982,48.9603,229.49648, 54000};
  //arma::Mat<double>::fixed<2,6> keplerianElements =  {{arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf},
  //                                                    {arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf}};
  posVelPair = mant::itd::orbitOnPosition(1234.0, keplerianElementsVec);

  std::cout << "pos = (" << posVelPair.first(0) << ", " << posVelPair.first(1) << ", " << posVelPair.first(2) << ")" << std::endl;
  std::cout << "vel = (" << posVelPair.second(0) << ", " << posVelPair.second(1) << ", " << posVelPair.second(2) << ")" << std::endl << std::endl;
*/
/*

  // gravity assist

  arma::Col<double>::fixed<3> in = {123.0, 246.0, 369.0};
  arma::Col<double>::fixed<3> out = {234.0, 135.0, 470.0};


  double vIn = arma::norm(in);
  double vOut = arma::norm(out);
  double alpha = std::acos(arma::dot(in, out) / (vIn * vOut));
  double DV;
  double rp;
  
  PowSwingByInv(vIn, vOut, alpha, DV, rp);
  
  std::cout << "gtoc : [DV , rp] = [" << DV << " , " << rp << "]" << std::endl;


  std::pair<double, double> res = mant::itd::gravityAssist(in, out);

  std::cout << "mantella : [DV , rp] = [" << res.first << " , " << res.second << "]" << std::endl;
*/

/*
  // lambert
  std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> lambert(
    const arma::Col<double>::fixed<3>& departurePosition,
    const arma::Col<double>::fixed<3>& arrivalPosition,
    const double transferTime) {


  std::vector<std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>>> resVec;
  
  resVec = mant::itd::lambert({-1.4349114384e11, -4.268153824e10, 0.0}, {8.77537144e9, 4.50685708e10, 2.85978352e9}, 100.0);
  
  for(auto pair : resVec) {
    std::cout << "lambert solutions" << std::endl;
    std::cout << pair.first << " - " << pair.second << std::endl;
  }
*/

  return 0;
}


