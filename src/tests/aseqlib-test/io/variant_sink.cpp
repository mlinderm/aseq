//
// Created by Michael Linderman on 12/16/15.
//

#include <sstream>
#include <gtest/gtest.h>

#include "aseq/io/variant.hpp"

using namespace aseq::io;

TEST(VCFSitesOnlySinkTest, WritesHeader) {
  std::stringstream sink_stream;
  auto sink = VariantSinkInterface::MakeVariantSink(FileFormat::VCF4_2, sink_stream);
  ASSERT_TRUE(sink);
}
