#include <gmock/gmock.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>

boost::filesystem::path test_inputs_g;

int main(int argc, char* argv[]) {
  ::testing::InitGoogleMock(&argc, argv);

  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;

  // Path to input test files
  if (argc >= 2) {
    test_inputs_g = argv[1];
  }

  return RUN_ALL_TESTS();
}
