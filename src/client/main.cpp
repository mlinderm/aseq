#include <iostream>

#include <docopt.h>
#include <glog/logging.h>
#include <boost/filesystem/path.hpp>

#include "aseq-version.h"
#include "aseq/aseq.hpp"
#include "aseq/io/variant.hpp"

namespace fs = boost::filesystem;

static const char USAGE[] = R"(aseq Sequencing analysis toolkit

Usage:
  aseq <file1>
  aseq (-h | --help)

Options:
-h --help              Show this screen.
-R=<ref>, --ref=<ref>  Indexed reference FASTA file
)";

int main(int argc, char* argv[]) {
  std::map<std::string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc},
                                                             true,  // show help if requested
                                                             ASEQ_VERSION,  // version string
                                                             false          // options first
                                                             );
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;

  // Always log version
  LOG(INFO) << "aseq version: " << ASEQ_VERSION;

  using namespace aseq::io;

  auto source1 = VariantSourceInterface::MakeVariantSource(args["<file1>"].asString());
  auto sink = VariantSinkInterface::MakeVariantSink(FileFormat::VCF4_2, std::cout);

  while (auto v = source1->NextVariant()) {
    sink->PushVariant(v.get());
  }

  return 0;
}
