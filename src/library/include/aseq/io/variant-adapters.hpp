//
// Created by Michael Linderman on 4/5/16.
//

#pragma once

#include "aseq/util/attributes.hpp"
#include "aseq/io/variant.hpp"
#include "aseq/io/vcf.hpp"

namespace aseq {
namespace io {

class VariantMergeSourceInterface : public VariantSourceInterface {
 public:
  typedef std::unique_ptr<VariantMergeSourceInterface> FactoryResult;
  typedef util::Attributes::key_type SourceLabel;

  VariantMergeSourceInterface()
      : source_key_("SET"),
        source1_label_("source1"),
        source2_label_("source2"),
        merged_label_("source1-source2") {}

  static FactoryResult MakeMergeVariantSource(VariantSourceInterface::FactoryResult &&source1,
                                              VariantSourceInterface::FactoryResult &&source2);

  const util::Attributes::key_type &source_key() const { return source_key_; }
  const SourceLabel &source1_label() const { return source1_label_; }
  const SourceLabel &source2_label() const { return source2_label_; }
  const SourceLabel &merged_label() const { return merged_label_; }

 protected:
  util::Attributes::key_type source_key_;
  SourceLabel source1_label_, source2_label_, merged_label_;
};
}
}