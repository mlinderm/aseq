#include <iostream>

#include <docopt.h>
#include <glog/logging.h>
#include <cppformat/format.h>

#include "aseq-version.h"

#include "commands.hpp"

static const char USAGE[] = R"(aseq Sequencing analysis toolkit

Usage: aseq [-h|--help] <command> [<args>...]

Options:
  -h --help              Show this screen.

The aseq commands are:
   variants              Analyze variants

See 'aseq help <command>' for more information on a specific command.

)";

// Create map of commands
std::map<std::string, int (*)(const std::vector<std::string>&)> kCommands{
    {"variants", &VariantsMain}};

int main(int argc, char* argv[]) {
  std::map<std::string, docopt::value> args =
      docopt::docopt(USAGE, {argv + 1, argv + argc}, true, ASEQ_VERSION, true);

  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;

  std::string command(args["<command>"].asString());

  auto itr = kCommands.find(command);
  if (itr != kCommands.end()) {
    std::vector<std::string> sub_argv(args["<args>"].asStringList());
    sub_argv.insert(sub_argv.begin(), itr->first);
    return itr->second(sub_argv);
  } else {
    fmt::print(std::cerr, "'{}' is not a valid command\n", command);
    std::cerr << USAGE;
  }

  return 1;

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
}
