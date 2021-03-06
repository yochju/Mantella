#include "mantella_bits/numberTheory.hpp"

// C++ standard library
#include <limits>
#include <stdexcept>

// Mantella
#include "mantella_bits/assert.hpp"

namespace mant {
  arma::mat hammersleySet(
      const arma::uvec& bases,
      const arma::uvec& seeds,
      const arma::uword numberOfElements) {
    if (bases.is_empty()) {
      throw std::invalid_argument("hammersleySet: The bases must not be empty.");
    } else if (arma::size(bases) != arma::size(seeds)) {
      throw std::logic_error("hammersleySet: The number of bases must be equal to the number of seeds.");
    } else if (arma::any(bases < 2)) {
      throw std::domain_error("hammersleySet: All bases must be greater than 1.");
    }

    if (numberOfElements == 0) {
      return {};
    }

    if (!isRepresentableAsFloatingPoint(seeds.n_elem)) {
      throw std::overflow_error("hammersleySet: The number of seeds must be representable as a floating point.");
    } else if (!isRepresentableAsFloatingPoint(numberOfElements)) {
      throw std::overflow_error("hammersleySet: The number of elements must be representable as a floating point.");
    }

    // A Hammersley set is generated by producing a Van der Corput sequence for each dimension and skipping the last element of each dimension.
    // The last element is than set to n / `numberOfElements` (n = the dimension).
    arma::mat sequence(seeds.n_elem, numberOfElements);
    for (arma::uword n = 0; n < sequence.n_rows; ++n) {
      if (numberOfElements > 1) {
        sequence.submat(n, 0, n, numberOfElements - 2) = vanDerCorputSequence(bases(n), seeds(n), numberOfElements - 1).t();
      }
      sequence(n, numberOfElements - 1) = static_cast<double>(n + 1) / static_cast<double>(numberOfElements);
    }

    return sequence;
  }

  arma::mat haltonSequence(
      const arma::uvec& bases,
      const arma::uvec& seeds,
      const arma::uword numberOfElements) {
    if (bases.is_empty()) {
      throw std::invalid_argument("hammersleySet: The bases must not be empty.");
    } else if (arma::size(bases) != arma::size(seeds)) {
      throw std::logic_error("haltonSequence: The number of bases must be equal to the number of seeds.");
    } else if (arma::any(bases < 2)) {
      throw std::domain_error("haltonSequence: All bases must be greater than 1.");
    }

    if (numberOfElements == 0) {
      return {};
    }

    // A Halton sequence is generated by producing a Van der Corput sequence for each dimension.
    arma::mat sequence(seeds.n_elem, numberOfElements);
    for (arma::uword n = 0; n < sequence.n_rows; ++n) {
      sequence.row(n) = vanDerCorputSequence(bases(n), seeds(n), numberOfElements).t();
    }

    return sequence;
  }

  arma::vec vanDerCorputSequence(
      const arma::uword base,
      const arma::uword seed,
      const arma::uword numberOfElements) {
    if (base < 2) {
      throw std::domain_error("vanDerCorputSequence: The base must be greater than 1.");
    }

    if (numberOfElements == 0) {
      return {};
    }

    /* The *n*-th number of a van der Corput sequence can be calculated as follows:
     *
     * infinity                    
     *   sum    d_i * `base`^(-i)
     *  i = 1
     *
     * with *d_i* being the *i*-th digit of the *n*-th `base`-adic number.
     */

    if (std::numeric_limits<arma::uword>::max() - numberOfElements < seed) {
      throw std::overflow_error("vanDerCorputSequence: The seed plus the number of elements must be less than or equal to the largest supported integer.");
    } else if (!isRepresentableAsFloatingPoint(base)) {
      throw std::overflow_error("vanDerCorputSequence: The base must be representable as a floating point.");
    }

    arma::vec sequence(numberOfElements);
    for (arma::uword n = 0; n < numberOfElements; ++n) {
      double number = 0.0;

      double inverseDigitBase = 1.0 / static_cast<decltype(inverseDigitBase)>(base);
      arma::uword remaining = n + seed;
      while (remaining > 0) {
        number += static_cast<decltype(number)>(remaining % base) * inverseDigitBase;

        inverseDigitBase /= static_cast<decltype(inverseDigitBase)>(base);
        remaining /= base;
      }

      sequence(n) = number;
    }

    return sequence;
  }
}
