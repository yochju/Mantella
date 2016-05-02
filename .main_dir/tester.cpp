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
/*
  arma::Mat<double>::fixed<2, 6> earth 
    = { {1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0},
        {0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0}};
  
  arma::Mat<double>::fixed<2, 6> venus 
    = { {0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255},
        {0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418}};
                                          
  arma::Mat<double>::fixed<2, 6> jupiter 
    = { {5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909},
        {-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106}};
                                          
  arma::Mat<double>::fixed<2, 6> saturn 
    = { {9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448},
        {-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794}};
                                          
  arma::Mat<double>::fixed<2, 6> astroid2001TW229 
    = { {1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891},
        {0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343}};
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
  // eph test for venus
  //general
  ofstream mantPos, mantVel, gtocPos, gtocVel;
  mantPos.open("mantPos.txt");
  mantVel.open("mantVel.txt");
  gtocPos.open("gtocPos.txt");
  gtocVel.open("gtocVel.txt");
  
  double pos[3], vel[3];
  std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> pair;
  
  for(double i = 0.0; i < 1000.0; i = i + 10.0){
    //mant
    arma::Col<double>::fixed<7> venus = {7.233268496749391E-01, 6.755697267164094E-03, 3.394589632336535E+00, 5.518541455452200E+01, 7.667837563371675E+01, 4.931425178852966E+01, 51544.0};  
    
    pair = mant::itd::orbitOnPosition(i + 51544.0, venus);
    
    mantPos << i << " " << pair.first(0) << " " << pair.first(1) << " " << pair.first(2) << "\n";
    mantVel << i << " " << pair.second(0) << " " << pair.second(1) << " " << pair.second(2) << "\n";
      
    //gtoc   
    Planet_Ephemerides_Analytical (i, 2, pos, vel);
    gtocPos << i << " " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    gtocVel << i << " " << vel[0] << " " << vel[1] << " " << vel[2] << "\n";      
  }
  
  mantPos.close();
  mantVel.close();
  gtocPos.close();
  gtocVel.close();

  */ 
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

  
  
  std::cout << "---MANTELLA Mat for Venus" << std::endl;
  arma::Mat<double>::fixed<2,6> keplerianElementsMat =  { {0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255},
                                                          {0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418}};
  //arma::Mat<double>::fixed<2,6> keplerianElements =  {{arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf},
  //                                                    {arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf, arma::datum::inf}};
  std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> posVelPair = mant::itd::orbitOnPosition(1234.0, keplerianElementsMat);

  std::cout << "pos = (" << posVelPair.first(0) << ", " << posVelPair.first(1) << ", " << posVelPair.first(2) << ")" << std::endl;
  std::cout << "vel = (" << posVelPair.second(0) << ", " << posVelPair.second(1) << ", " << posVelPair.second(2) << ")" << std::endl << std::endl;
  
  
  
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
    // Mantella test
    mant::OptimisationProblem optProb(8);

    arma::Col<double> lower = {3000, 14, 14, 14, 14, 100, 366, 300};
    arma::Col<double> upper = {10000, 2000, 2000, 2000, 2000, 9000, 9000, 9000};
    optProb.setLowerBounds(lower);
    optProb.setUpperBounds(upper);


    std::vector<double> rp = {6352, 6779, 6352, 6779, 600001, 70001};
    optProb.setObjectiveFunction([&rp](const arma::Col<double>& parameter) {

        std::vector<double> variables = {
                parameter(0),
                parameter(1),
                parameter(2),
                parameter(3),
                parameter(4),
                parameter(5),
                parameter(6),
                parameter(7)
        };
        if(arma::accu(parameter) - parameter(0) > 30*365)
            return std::numeric_limits<double>::max();

        return gtoc1(variables, rp);
    });

    mant::RandomSearch algo;
    algo.setMaximalNumberOfIterations(1000);
    algo.optimise(optProb);
*/
/*
    arma::Col<double> bestParameter = algo.getBestParameter();

    std::cout << "Best value:\n\t" << std::defaultfloat << algo.getBestObjectiveValue() << std::endl;
    std::cout << "Best parameters:\n" << std::defaultfloat << bestParameter << std::endl; //%(upper-lower) + lower << std::endl;

    //Final data in outputPos.dat and outputVel.dat
    std::vector<double> bestSolutionVariables = {bestParameter(0), bestParameter(1),bestParameter(2),bestParameter(3),bestParameter(4),bestParameter(5),bestParameter(6),bestParameter(7)};
    gtoc1(bestSolutionVariables, rp);


    //std::vector<double> x = {6809.476683160, 169.598512787, 1079.375156244, 56.53776494142, 1044.014046276, 3824.160968179, 1042.885114734, 3393.057868710};
    //std::vector<double> x = {9000, 14, 14, 14, 14, 100, 366, 300};
    //std::cout << "Competition solution: " << gtoc1(x, rp) << std::endl;

//    allPlanetDataToFile(1, 40000, "planetData", 1);

    std::ofstream parameterFileHandler;
    parameterFileHandler.open("bestParameter.dat", std::ofstream::out | std::ofstream::trunc);

    parameterFileHandler << (upper(0)-lower(0))*bestParameter(0)+lower(0) << " ";
    parameterFileHandler << (upper(1)-lower(1))*bestParameter(1)+lower(1) << " ";
    parameterFileHandler << (upper(2)-lower(2))*bestParameter(2)+lower(2) << " ";
    parameterFileHandler << (upper(3)-lower(3))*bestParameter(3)+lower(3) << " ";
    parameterFileHandler << (upper(4)-lower(4))*bestParameter(4)+lower(4) << " ";
    parameterFileHandler << (upper(5)-lower(5))*bestParameter(5)+lower(5) << " ";
    parameterFileHandler << (upper(6)-lower(6))*bestParameter(6)+lower(6) << " ";
    parameterFileHandler << (upper(7)-lower(7))*bestParameter(7)+lower(7);

    parameterFileHandler.close();


    arma::Col<double> lin1 = arma::linspace<arma::Col<double>>(14, 2000, 100);
    arma::Col<double> lin2 = arma::linspace<arma::Col<double>>(300, 9000, 100);

    GridAnalysis analyser;
    analyser.linAnalysis(lin1, lin2);
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


