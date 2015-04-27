namespace mant {

  class RankingSVM {
    public:     
      
    protected:
      DistPow2_Euclidean(arma::Col<double> x1, arma::Col<double> x2);
      
      double sigmaA = 1;
      double sigmaPow = 1;
      double epsilon = 1;
  };

  //
  // Implementation
  //

}