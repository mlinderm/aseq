//
// Created by Michael Linderman on 4/5/16.
//

#include <algorithm>
#include <iostream>
#include <unordered_set>

#include <aseq/io/variant.hpp>
#include <aseq/io/vcf.hpp>
#include "aseq/io/variant-adapters.hpp"
#include "aseq/algorithm/variant.hpp"

namespace aseq {
namespace io {

namespace impl {
class VCFMergeSource : public VariantMergeSourceInterface {
  using Attributes = util::Attributes;

 public:
  VCFMergeSource() = delete;
  VCFMergeSource(VariantSourceInterface::FactoryResult &&source1,
                 VariantSourceInterface::FactoryResult &&source2)
      : source1_(std::move(source1)), source2_(std::move(source2)) {
    if (!source1_ || !source2_)
      throw util::invalid_argument() << util::error_message("Invalid source supplied as argument");

    const VCFHeader &header1 = static_cast<const VCFHeader &>(source1_->header()),
                    &header2 = static_cast<const VCFHeader &>(source2_->header());

    header_ = header1;
    header_.AddFields(header2);
    {  // Add set key
      auto r = header_.AddINFOField(VCFHeader::Field(source_key_, 1, VCFHeader::Field::Type::STRING,
                                                     "Source of merged variant"));
      if (!r.second) {
        throw util::incompatible_header_attribute()
            << util::error_message("Source key already exists in files to be merged");
      }
    }
    {  // Merge samples in new header
      std::vector<model::Sample> samples1(header1.samples()), samples2(header2.samples());
      std::sort(samples1.begin(), samples1.end());
      std::sort(samples2.begin(), samples2.end());

      std::set_difference(samples1.begin(), samples1.end(), samples2.begin(), samples2.end(),
                          std::back_inserter(source1only_samples_));
      std::set_difference(samples2.begin(), samples2.end(), samples1.begin(), samples1.end(),
                          std::back_inserter(source2only_samples_));

      std::vector<model::Sample> merged_samples;
      std::set_union(samples1.begin(), samples1.end(), samples2.begin(), samples2.end(),
                     std::back_inserter(merged_samples));
      header_.SetSamples(merged_samples.begin(), merged_samples.end());
    }

    // Initialize initial variants
    variant1_ = source1_->NextVariant();
    variant2_ = source2_->NextVariant();
  }

  FileFormat file_format() const override { return header_.file_format(); }
  const VCFHeader &header() const override { return header_; }

  virtual bool IsIndexed() const override { return source1_->IsIndexed() && source2_->IsIndexed(); }

  virtual void SetRegion(model::Contig contig, model::Pos pos, model::Pos end) override {
    source1_->SetRegion(contig, pos, end);
    source2_->SetRegion(contig, pos, end);

    // Need to reinitialize variants after setting region
    variant1_ = source1_->NextVariant();
    variant2_ = source2_->NextVariant();
  }

  virtual NextResult NextVariant() override {
    if (variant1_ && !variant2_) {
      return ModifyAndGetNext(variant1_, source1_, source1_label_, source2only_samples_);
    } else if (!variant1_ && variant2_) {
      return ModifyAndGetNext(variant2_, source2_, source2_label_, source1only_samples_);
    } else if (variant1_ && variant2_) {
      model::CompareVariants cmp;
      using result_type = model::CompareVariants::result_type;
      switch (cmp(*variant1_, *variant2_)) {
        default:
          throw util::invalid_argument()
              << util::error_message("Merging multi-allelic variants is not supported");
        case result_type::EQUAL: {
          model::VariantContext result = std::move(*variant1_);

          {  // Set merged QUAL as the lesser of the defined QUALs
            const auto &qual2 = variant2_->qual();
            if (qual2 && (!result.qual() || qual2 < result.qual())) result.SetQual(qual2);
          }
          // TODO: Merge FILT

          MergeINFOAttributes(result.attributes(), std::move(variant2_->attributes()));
          result.SetAttribute(source_key_, merged_label_);

          result.MergeGenotypes(variant2_->genotypes());

          variant1_ = source1_->NextVariant();
          variant2_ = source2_->NextVariant();
          return NextResult(std::move(result));
        }
        case result_type::BEFORE:
          return ModifyAndGetNext(variant1_, source1_, source1_label_, source2only_samples_);
        case result_type::AFTER:
          return ModifyAndGetNext(variant2_, source2_, source2_label_, source1only_samples_);
      }
    } else
      return NextResult();
  }

  static FactoryResult MakeMergeVariantSource(FactoryResult &&source1, FactoryResult &&source2);

 private:
  VariantSourceInterface::FactoryResult source1_, source2_;
  NextResult variant1_, variant2_;
  VCFHeader header_;
  std::vector<model::Sample> source1only_samples_, source2only_samples_;

  NextResult ModifyAndGetNext(NextResult &variant, VariantSourceInterface::FactoryResult &source,
                              const SourceLabel &label,
                              const std::vector<model::Sample> &no_call_samples) {
    model::VariantContext result = std::move(*variant);
    variant = source->NextVariant();

    result.SetAttribute(source_key_, label);
    for (auto &s : no_call_samples) {  // Add "no call" samples (from other source)
      result.AddGenotype(s, model::Genotype::kNone);
    }

    return NextResult(std::move(result));
  }

  void MergeINFOAttributes(util::Attributes &dst, util::Attributes &&src) {
    for (auto d = dst.begin(); d != dst.end();) {
      auto s = src.find(d->first);
      if (s != src.end()) {
        // Delete "dst" attributes with differing values
        d = (d->second != s->second) ? dst.erase(d) : ++d;
        src.erase(s);  // "src" always deleted since even if equal it is duplicated
      }
    }
    dst.insert(std::make_move_iterator(src.begin()), std::make_move_iterator(src.end()));
  }
};
}

VariantMergeSourceInterface::FactoryResult VariantMergeSourceInterface::MakeMergeVariantSource(
    VariantSourceInterface::FactoryResult &&source1,
    VariantSourceInterface::FactoryResult &&source2) {
  if (IsVCF(source1->file_format()) && IsVCF(source2->file_format())) {
    return std::make_unique<impl::VCFMergeSource>(std::move(source1), std::move(source2));
  } else {
    throw util::invalid_argument() << util::error_message("merging only supported for VCF sources");
  }
}
}
}