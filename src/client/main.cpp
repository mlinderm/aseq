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
}
