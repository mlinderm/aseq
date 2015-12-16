//
// Created by Michael Linderman on 12/16/15.
//

#include "aseq/io/vcf.hpp"

namespace aseq {
namespace io {

VCFSink::VCFSink(const VCFHeader &header, Writer &&writer)
    : header_(header), writer_(std::forward<Writer>(writer)) {}

void VCFSink::PushVariant(const model::VariantContext &context) {}

}  // namespace io
}  // namespace aseq