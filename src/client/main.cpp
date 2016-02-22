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
  aseq <file1> <file2>
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
  auto source2 = VariantSourceInterface::MakeVariantSource(args["<file2>"].asString());
  auto sink = VariantSinkInterface::MakeVariantSink(FileFormat::VCF4_2, std::cout);
  aseq::model::CompareVariants cmp;

  auto v1 = source1->NextVariant();
  auto v2 = source2->NextVariant();
  while (true) {
    if (!v1) {
      while (auto v = source2->NextVariant()) {
        std::cerr << v.get() << std::endl;
      }
      break;
    }
    if (!v2) {
      while (auto v = source1->NextVariant()) {
        std::cerr << v.get() << std::endl;
      }
      break;
    }
    switch (cmp(v1.get(), v2.get())) {
      case aseq::model::CompareVariants::result_type::BEFORE:
        std::cerr << v1.get() << std::endl;
        v1 = source1->NextVariant();
        break;
      case aseq::model::CompareVariants::result_type::AFTER:
        std::cerr << v2.get() << std::endl;
        v2 = source2->NextVariant();
        break;
      default:
        std::cerr << v1.get() << "," << v2.get() << std::endl;
        v1 = source1->NextVariant();
        v2 = source2->NextVariant();
        break;
    }
  }

  return 0;
}
