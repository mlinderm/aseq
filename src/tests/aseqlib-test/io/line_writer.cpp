//
// Created by Michael Linderman on 12/16/15.
//

#include <fstream>
#include <istream>

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include <cppformat/format.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "aseq/io/line.hpp"

namespace fs = boost::filesystem;
using namespace aseq::io;

TEST(ASCIIStreamWriterTest, WritesLinesToASCIIStream) {
  std::stringstream sink_stream;

  auto writer = ASCIILineWriterInterface::MakeLineWriter(sink_stream);
  ASSERT_TRUE(writer);

  for (int i = 1; i <= 3; i++) {
    writer->Write(fmt::format("line{}\n", i));
  }
  EXPECT_EQ("line1\nline2\nline3\n", sink_stream.str());
}

class BGZipLineWriterTest : public ::testing::Test {
 protected:
  BGZipLineWriterTest() : directory_(fs::unique_path()) {
    fs::create_directory(directory_);
    file_ = directory_ / "test.gz";
  }

  virtual void SetUp() { ASSERT_TRUE(fs::exists(directory_)); }
  virtual void TearDown() { fs::remove_all(directory_); }

  fs::path directory_;
  fs::path file_;
};

TEST_F(BGZipLineWriterTest, WriteLinesToBGZipFile) {
  using namespace boost::iostreams;

  {
    auto writer = ASCIILineWriterInterface::MakeLineWriter(file_);
    ASSERT_TRUE(writer);
    for (int i = 1; i <= 3; i++) {
      writer->Write(fmt::format("line{}\n", i));
    }
  }

  std::ifstream file(file_.native(), std::ios_base::in | std::ios_base::binary);
  filtering_istream in;
  in.push(gzip_decompressor());
  in.push(file);

  for (int i = 1; i <= 3; i++) {
    std::string line;
    std::getline(in, line);
    EXPECT_EQ(fmt::format("line{}", i), line);
  }
}