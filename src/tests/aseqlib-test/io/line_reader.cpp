//
// Created by Michael Linderman on 12/12/15.
//

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include "aseq/io/line.hpp"

namespace fs = boost::filesystem;
using namespace aseq::io;

TEST(ASCIIStreamLineReaderTest, ReadsLinesFromASCIIStream) {
  std::stringstream content("line1\nline2\nline3\n");

  auto reader = ASCIILineReaderInterface::MakeLineReader(content);
  ASSERT_TRUE(reader);

  for (int i = 1; i <= 3; i++) {
    auto line = reader->ReadNextLine();
    EXPECT_TRUE(line);
    EXPECT_EQ(std::string("line") + std::to_string(i), boost::copy_range<std::string>(*line));
  }
  EXPECT_FALSE(reader->ReadNextLine());
}