namespace mant {
  template <typename T>
  class HCMA : public PopulationBasedOptimisationAlgorithm<T> {
    public:
      explicit HCMA(
          const std::shared_ptr<OptimisationProblem<T>> optimisationProblem,
          const unsigned int populationSize) noexcept;
      
          std::string toString() const noexcept override;
      protected:

      void optimiseImplementation() override;
  };
  
  template <typename T>
  void HCMA<T>::optimiseImplementation() {
    
  }
  
  template <typename T>
  std::string HCMA<T>::toString() const noexcept {
    return "cuckoo_search";
  }
}