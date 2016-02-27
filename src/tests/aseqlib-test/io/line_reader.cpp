//
// Created by Michael Linderman on 12/12/15.
//

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include <cppformat/format.h>

#include "aseq/io/line.hpp"

namespace fs = boost::filesystem;
using namespace aseq::io;

extern fs::path test_inputs_g;

TEST(ASCIIStreamLineReaderTest, ReadsLinesFromASCIIStream) {
  std::stringstream content("line1\nline2\nline3\n");

  auto reader = ASCIILineReaderInterface::MakeLineReader(content);
  ASSERT_TRUE(reader);

  for (int i = 1; i <= 3; i++) {
    auto line = reader->ReadNextLine();
    EXPECT_TRUE(line);
    EXPECT_EQ(fmt::format("line{}", i), boost::copy_range<std::string>(*line));
  }
  EXPECT_FALSE(reader->ReadNextLine());
}

class TabixLineReaderTest : public ::testing::Test {
 protected:
  TabixLineReaderTest() : file_(test_inputs_g) { file_ /= "three_line.vcf.gz"; }

  virtual void SetUp() { ASSERT_TRUE(fs::exists(file_)); }

  fs::path file_;
};

TEST_F(TabixLineReaderTest, ReadsLinesFromTabixFile) {
  auto reader = ASCIILineReaderInterface::MakeLineReader(file_);
  ASSERT_TRUE(reader);

  for (int i = 1; i <= 3; i++) {
    auto line = reader->ReadNextLine();
    EXPECT_TRUE(line);
    EXPECT_GT(line->size(), 0);
  }
  EXPECT_FALSE(reader->ReadNextLine());
}