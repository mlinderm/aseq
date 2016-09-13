//
// Created by Michael Linderman on 12/13/15.
//

#include <ostream>

#include "aseq/io/vcf.hpp"

namespace aseq {
namespace io {

const VCFHeader::Field::Number VCFHeader::Field::R, VCFHeader::Field::A, VCFHeader::Field::G,
    VCFHeader::Field::UNBOUNDED;

#define FILTER_FIELD(id, description) \
  const VCFHeader::Field VCFHeader::FILTER::id(#id, description);
#define INFO_FIELD(id, number, type, description) \
  const VCFHeader::Field VCFHeader::INFO::id(#id, number, VCFHeader::type, description);
#define FORMAT_FIELD(id, number, type, description) \
  const VCFHeader::Field VCFHeader::FORMAT::id(#id, number, VCFHeader::type, description);

#include "aseq/io/vcf_fields.def"

#undef FORMAT_FIELD
#undef INFO_FIELD
#undef FILTER_FIELD

VCFHeader &VCFHeader::AddFields(const VCFHeader &other) {
#define FIELD_COPY(KIND)                               \
  for (auto &field : other.KIND##Values()) {           \
    auto r = Add##KIND##Field(field);                  \
    if (!r.second && !r.first.IsCompatible(field)) {   \
      throw new util::incompatible_header_attribute(); \
    }                                                  \
  }

  FIELD_COPY(INFO);
  FIELD_COPY(FORMAT);
  FIELD_COPY(FILTER);

#undef FIELD_COPY

  return *this;
}

std::pair<const VCFHeader::Field &, bool> VCFHeader::AddField(VCFHeader::Fields &fields,
                                                              const VCFHeader::Field &field) {
  auto r = fields.emplace(field.id_, field);
  return std::make_pair(r.first->second, r.second);
}

bool VCFHeader::HasField(const Fields &fields, const Fields::key_type &key) {
  return fields.find(key) != fields.end();
}

std::ostream &operator<<(std::ostream &os, const VCFHeader::Field &field) {
  return (os << field.id_);
}

void VCFHeader::SetSitesOnly() {
  FORMAT_.clear();
  samples_.clear();
}

}  // namespace io
}  // namespace aseq
