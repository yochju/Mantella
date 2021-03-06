#pragma once

// Armadillo
#include <armadillo>

namespace mant {
  arma::mat hammersleySet(
      const arma::uvec& base,
      const arma::uvec& seed,
      const arma::uword numberOfElements);

  arma::mat haltonSequence(
      const arma::uvec& base,
      const arma::uvec& seed,
      const arma::uword numberOfElements);

  arma::vec vanDerCorputSequence(
      const arma::uword base,
      const arma::uword seed,
      const arma::uword numberOfElements);
}
