//
// Created by Michael Linderman on 2/28/16.
//

#include <gtest/gtest.h>

#include "aseq/io/fasta.hpp"

using namespace aseq::io;

TEST(FastaSinkTest, WritesFastaEntryWithLineWrapping) {
  std::stringstream sink_stream;

  FastaSink sink(sink_stream, 5);
  sink.PushSequence("ref", "ABCDEFGHIJKL");

  EXPECT_EQ(std::string(">ref\nABCDE\nFGHIJ\nKL\n"), sink_stream.str());
}
