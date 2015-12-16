//
// Created by Michael Linderman on 12/16/15.
//

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include "aseq/io/line.hpp"

namespace fs = boost::filesystem;
using namespace aseq::io;

TEST(ASCIIStreamWriterTest, WritesLinesToASCIIStream) {
  std::stringstream sink_stream;

  auto writer = ASCIILineWriterInterface::MakeLineWriter(sink_stream);
  ASSERT_TRUE(writer);

  for (int i = 1; i <= 3; i++) {
    writer->WriteLine(std::string("line") + std::to_string(i));
  }
  EXPECT_EQ("line1\nline2\nline3\n", sink_stream.str());
}