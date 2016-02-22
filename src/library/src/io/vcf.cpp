//
// Created by Michael Linderman on 12/13/15.
//

#include "aseq/io/vcf.hpp"

namespace aseq {
namespace io {

const VCFHeader::Field::Number VCFHeader::Field::R, VCFHeader::Field::A, VCFHeader::Field::G,
    VCFHeader::Field::UNBOUNDED;

#define INFO_FIELD(id, number, type, description) \
  const VCFHeader::Field VCFHeader::INFO::id(#id, number, VCFHeader::type, description);
#define FILTER_FIELD(id, description) \
  const VCFHeader::Field VCFHeader::FILTER::id(#id, description);
#include "aseq/io/vcf_fields.def"
#undef FILTER_FIELD
#undef INFO_FIELD

std::pair<const VCFHeader::Field &, bool> VCFHeader::AddField(VCFHeader::Fields &fields,
                                                              const VCFHeader::Field &field) {
  auto r = fields.emplace(field.id_, field);
  return std::make_pair(r.first->second, r.second);
}

bool VCFHeader::HasField(const Fields &fields, const Fields::key_type &key) {
  return fields.find(key) != fields.end();
}
}  // namespace io
}  // namespace aseq
