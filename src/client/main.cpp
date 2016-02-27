#include <iostream>

#include <docopt.h>
#include <glog/logging.h>
#include <boost/filesystem/path.hpp>
#include <cppformat/format.h>

#include "aseq-version.h"
#include "aseq/aseq.hpp"
#include "aseq/io/variant.hpp"
#include "aseq/io/reference.hpp"

namespace fs = boost::filesystem;

static const char USAGE[] = R"(aseq Sequencing analysis toolkit

Usage: aseq [-h|--help] <command> [<args>...]

Options:
  -h --help              Show this screen.

The aseq commands are:
   variants              Analyze variants

See 'aseq help <command>' for more information on a specific command.

)";

class CommandInterface {
 public:
  virtual int Main(const std::vector<std::string>& argv) = 0;
};

class VariantsCommands : public CommandInterface {
 public:
  virtual const char* Usage() const {
    return R"(aseq Variant analysis commands

Usage:
  aseq variants (-h | --help)
  aseq variants consensus [--flank F] -R <ref> <file>

Options:
  -h --help                  Show this screen.
  -R <ref>, --ref <ref>      Reference fasta (indexed)
  --flank F                  Length of flanking sequence [default: 1000]
)";
  }

  virtual int Main(const std::vector<std::string>& argv) override;
};

int main(int argc, char* argv[]) {
  std::map<std::string, docopt::value> args =
      docopt::docopt(USAGE, {argv + 1, argv + argc}, true, ASEQ_VERSION, true);

  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;

  // Always log version
  LOG(INFO) << "aseq version: " << ASEQ_VERSION;

  if (args["<command>"].asString() == "variants") {
    std::vector<std::string> sub_argv({std::string("variants")});
    auto a = args["<args>"].asStringList();
    sub_argv.insert(sub_argv.end(), a.begin(), a.end());
    return VariantsCommands().Main(sub_argv);
  } else {
    std::cerr << fmt::format("{} is not a valid command", args["<command>"]) << std::endl;
    std::cerr << USAGE;
    return 1;
  }

  /*
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
  */
  return 0;
}

int VariantsCommands::Main(const std::vector<std::string>& argv) {
  std::map<std::string, docopt::value> args = docopt::docopt(Usage(), argv, true, ASEQ_VERSION);
  using namespace aseq::io;

  try {
    ReferenceSource ref(args["--ref"].asString());
    auto source = VariantSourceInterface::MakeVariantSource(args["<file>"].asString());
    auto v = source->NextVariant();
    if (!v) {
      LOG(ERROR) << "No variants found in input";
      return 1;
    }

    std::cout << ">ref" << std::endl << ref.Sequence(v->contig(), v->pos(), v->end()) << std::endl;
    if (source->NextVariant()) {
      LOG(WARNING) << "Multiple variants in input, but FASTA only emitted for the first";
    }
  } catch (aseq::util::exception_base& e) {
    LOG(ERROR) << e.what();
    return 1;
  }

  return 0;
}
