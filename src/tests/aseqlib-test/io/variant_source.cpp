//
// Created by Michael Linderman on 12/13/15.
//

#include <fstream>

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

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

  EXPECT_TRUE(header.FILTER_HasField(VCFHeader::FILTER::PASS));

  EXPECT_EQ(3, boost::size(header.INFO_Values()));
  ASSERT_TRUE(header.INFO_HasField(VCFHeader::INFO::AN));
}

INSTANTIATE_TEST_CASE_P(VCFSitesOnly, VCFSitesOnlySourceTest,
                        ::testing::Values("sites_only.vcf", "sites_only.vcf.gz"));