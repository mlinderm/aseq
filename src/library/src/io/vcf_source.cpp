//
// Created by Michael Linderman on 12/13/15.
//

#include "aseq/io/vcf.hpp"

namespace aseq {
namespace io {

VCFSource::VCFSource(FileFormat format, VCFSource::Reader&& reader)
    : header_(format), reader_(std::forward<Reader>(reader)) {}

FileFormat VCFSource::file_format() const { return header_.file_format(); }

VariantSourceInterface::NextResult VCFSource::NextVariant() {
  return boost::optional<model::VariantContext>();
}
}  // namespace io
}  // namespace aseq