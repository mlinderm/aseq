//
// Created by Michael Linderman on 2/22/16.
//

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include "aseq/io/reference.hpp"

using namespace aseq::io;

namespace fs = boost::filesystem;
extern fs::path test_inputs_g;

class ReferenceSourceTest : public ::testing::Test {
 protected:
  ReferenceSourceTest() : file_(test_inputs_g) { file_ /= "reference.fa"; }

  virtual void SetUp() { ASSERT_TRUE(fs::exists(file_)); }

  fs::path file_;
};

TEST_F(ReferenceSourceTest, ReadsSequenceFromFa) {
  ASSERT_NO_THROW({
    ReferenceSource source(file_);
    auto seq = source.Sequence("chr1", 2, 10);
    EXPECT_EQ("ATCACAGGT", seq);
  });
}
