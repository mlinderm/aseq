//
// Created by Michael Linderman on 12/16/15.
//

#include <sstream>

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include "aseq/io/variant.hpp"

using namespace aseq::io;
namespace fs = boost::filesystem;

TEST(VCFSitesOnlySinkTest, WritesHeader) {
  std::stringstream sink_stream;
  auto sink = VariantSinkInterface::MakeVariantSink(FileFormat::VCF4_2, sink_stream);
  ASSERT_TRUE(sink);
}

class TabixVCFWriterTest : public ::testing::Test {
 protected:
  TabixVCFWriterTest() : directory_(fs::unique_path()) {
    fs::create_directory(directory_);
    file_ = directory_ / "test.vcf.gz";
  }

  virtual void SetUp() { ASSERT_TRUE(fs::exists(directory_)); }
  virtual void TearDown() { fs::remove_all(directory_); }

  fs::path directory_;
  fs::path file_;
};

TEST_F(TabixVCFWriterTest, WritesTabixIndex) {
  {
    auto sink = VariantSinkInterface::MakeVariantSink(FileFormat::VCF4_2, file_);
    ASSERT_TRUE(sink);
  }
  auto index = fs::path(file_) += ".tbi";
  ASSERT_TRUE(fs::exists(index));

  EXPECT_NO_THROW({
    auto source = VariantSourceInterface::MakeVariantSource(file_);
    ASSERT_TRUE(source);
    EXPECT_EQ(FileFormat::VCF4_2, source->file_format());
  });
}