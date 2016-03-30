//
// Created by Michael Linderman on 12/16/15.
//

#pragma once

#include <boost/range/iterator_range.hpp>

#include "aseq/util/flyweight.hpp"

namespace aseq {
namespace model {

struct allele_tag {};

class Allele : public util::flyweight_string_no_track<allele_tag> {
  typedef util::flyweight_string_no_track<allele_tag> BaseType;

 public:
  typedef typename BaseType::initializer initializer;

  static const Allele A, G, C, T, N;
  static const Allele MISSING;
  static const Allele NON_REF;

  using BaseType::BaseType;

  Allele() : BaseType(MISSING) {}
  Allele(const Allele& allele) = default;
  Allele(Allele&& allele) = default;

  Allele& operator=(const Allele& allele) = default;
  Allele& operator=(Allele&& allele) = default;

  bool IsSymbolic() const;
  Allele SubAllele(size_t pos = 0, size_t len = std::string::npos) const;
};

inline Allele operator+(const std::string& lhs, Allele rhs) { return Allele(lhs + rhs.get()); }
inline Allele operator+(Allele lhs, const std::string& rhs) { return Allele(lhs.get() + rhs); }

}  // namespace model
}  // namespace aseq