//
// Created by Michael Linderman on 12/13/15.
//

#include <fstream>

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include "aseq/model/variant_context.hpp"
#include "aseq/model/genotype.hpp"
#include "aseq/io/variant.hpp"
#include "aseq/io/vcf.hpp"

namespace fs = boost::filesystem;
using namespace aseq::io;

extern fs::path test_inputs_g;

class VariantSourceTest : public ::testing::TestWithParam<const char*> {
 protected:
  VariantSourceTest() : file_(test_inputs_g) {}

  virtual void SetUp() {
    file_ /= GetParam();
    ASSERT_TRUE(fs::exists(file_));
  }

  fs::path file_;
};

class VCFSitesOnlySourceTest : public VariantSourceTest {};

TEST_P(VCFSitesOnlySourceTest, DetectsFileFormat) {
  auto source = VariantSourceInterface::MakeVariantSource(file_);
  ASSERT_TRUE(source);
  ASSERT_EQ(FileFormat::VCF4_2, source->file_format());

  auto* vcf_source = dynamic_cast<VCFSource*>(source.get());
  auto& header = vcf_source->header();

  EXPECT_EQ(3, boost::size(header.INFOValues()));

  // Check for standard fields
  EXPECT_TRUE(header.HasFILTERField(VCFHeader::FILTER::PASS));
  EXPECT_TRUE(header.HasFORMATField(VCFHeader::FORMAT::GT));
  EXPECT_TRUE(header.HasFORMATField(VCFHeader::FORMAT::FT));
}

TEST_P(VCFSitesOnlySourceTest, IteratesAndParsesVariants) {
  auto source = VariantSourceInterface::MakeVariantSource(file_);
  ASSERT_TRUE(source);

  VariantSourceInterface::NextResult r;

  // 1        11916414        rs149774989     C       T       1765.16 PASS    AC=3;AN=48
  r = source->NextVariant();
  ASSERT_TRUE(r);

  // 1        11916594        rs141879773     G       A       1699.71 PASS    AC=1;AN=48
  r = source->NextVariant();
  ASSERT_TRUE(r);

  // 1        11916764        rs79387574      C       G       5780.16 PASS    AC=4;AN=48
  r = source->NextVariant();
  ASSERT_TRUE(r);

  r = source->NextVariant();
  EXPECT_FALSE(r);  // There are three variants in the file
}

INSTANTIATE_TEST_CASE_P(VCFSitesOnly, VCFSitesOnlySourceTest,
                        ::testing::Values("sites_only.vcf", "sites_only.vcf.gz"));

class VCFSingleSampleSourceTest : public VariantSourceTest {};

TEST_P(VCFSingleSampleSourceTest, IteratesAndParsesVariants) {
  using aseq::model::Genotype;
  using aseq::model::Allele;
  auto source = VariantSourceInterface::MakeVariantSource(file_);
  ASSERT_TRUE(source);

  VariantSourceInterface::NextResult r;

  // chr1     100     rs1     A       T       100.0   PASS    AC=1;AN=2       GT      0/1
  r = source->NextVariant();
  ASSERT_TRUE(r);
  {
    EXPECT_TRUE(r->HasQual());
    EXPECT_TRUE(r->HasFilter() && r->IsPASSing());
    EXPECT_TRUE(r->IsBiallelic());
    EXPECT_EQ(Allele::A, r->ref());
    EXPECT_EQ(Allele::T, r->alt(0));

    auto& gt = r->GetGenotype("NA12878");
    EXPECT_EQ(Genotype::kRefAlt, gt.alleles());
  }

  r = source->NextVariant();
  EXPECT_FALSE(r);  // There is one variant in the file
}

INSTANTIATE_TEST_CASE_P(VCFSingleSample, VCFSingleSampleSourceTest,
                        ::testing::Values("single_sample.vcf", "single_sample.vcf.gz"));

class VCFSpecificationSourceTest : public VariantSourceTest {};

TEST_P(VCFSpecificationSourceTest, IteratesAndParsesVariants) {
  auto source = VariantSourceInterface::MakeVariantSource(file_);
  ASSERT_TRUE(source);

  auto* vcf_source = dynamic_cast<VCFSource*>(source.get());
  auto& header = vcf_source->header();
  EXPECT_TRUE(header.HasFILTERField("q10"));
  EXPECT_TRUE(header.HasFORMATField(VCFHeader::FORMAT::GQ));

  VariantSourceInterface::NextResult r;

  for (size_t i = 0; i < 5; i++) {
    r = source->NextVariant();
    ASSERT_TRUE(r);
    EXPECT_EQ(3, r->NumGenotypes());
  }

  r = source->NextVariant();
  EXPECT_FALSE(r);  // There are five variants in the file
}

INSTANTIATE_TEST_CASE_P(VCFSpecification, VCFSpecificationSourceTest,
                        ::testing::Values("vcf_specification.vcf"));
