namespace mant {
  template <typename T>
  class HCMA : public PopulationBasedOptimisationAlgorithm<T> {
    public:
      explicit HCMA(
          const std::shared_ptr<OptimisationProblem<T>> optimisationProblem,
          const unsigned int populationSize) noexcept;
      
          std::string toString() const noexcept override;
      protected:
        
        unsigned int lambdaMult; //settings.lambdaMult
        double sigma0 = 2.0; //later found in cur_state.insigma

      void optimiseImplementation() override;
  };
  
  //
  //Implementation
  //
  
  template <typename T>
  HCMA<T>::HCMA(const std::shared_ptr<OptimisationProblem<T>> optimisationProblem,
      const unsigned int populationSize) noexcept
    : PopulationBasedOptimisationAlgorithm<T>(optimisationProblem, populationSize) {
    if(this->numberOfDimensions_ < 10) {
      lambdaMult = 1;
    } else {
      lambdaMult = floor(std::pow(10,1+std::log2(this->numberOfDimensions_*1.0/10)));
    }
    
  }
  
  template <typename T>
  void HCMA<T>::optimiseImplementation() {
    
  }
  
  template <typename T>
  std::string HCMA<T>::toString() const noexcept {
    return "cuckoo_search";
  }
}