//
// Created by Michael Linderman on 12/12/15.
//

#include <istream>
#include <fstream>

#include <boost/filesystem.hpp>

#include "aseq/util/exception.hpp"
#include "aseq/io/line.hpp"

using namespace aseq::util;
namespace fs = boost::filesystem;

namespace aseq {
namespace io {

namespace impl {

class ASCIIStreamLineReader : public ASCIILineReaderInterface {
 public:
  ASCIIStreamLineReader() = delete;

  explicit ASCIIStreamLineReader(std::istream &istream) : istream_(istream) {}

  explicit ASCIIStreamLineReader(const fs::path &path)
      : owned_istream_(std::make_unique<std::ifstream>(path.native())), istream_(*owned_istream_) {}

  virtual NextResult ReadNextLine() override {
    if (!std::getline(istream_, line_))
      return NextResult();
    else
      return NextResult(
          boost::make_iterator_range<Line::iterator>(line_.data(), line_.data() + line_.size()));
  }

 private:
  std::unique_ptr<std::istream> owned_istream_;
  std::istream &istream_;
  std::string line_;
};

}  // impl namespace

ASCIILineReaderInterface::FactoryResult ASCIILineReaderInterface::MakeLineReader(
    std::istream &istream) {
  return std::make_unique<impl::ASCIIStreamLineReader>(istream);
}

ASCIILineReaderInterface::FactoryResult ASCIILineReaderInterface::MakeLineReader(
    const fs::path &path) {
  if (!fs::exists(path)) {  // Catch a common error case
    throw file_parse_error() << error_message("input file does not exist");
  }

  return std::make_unique<impl::ASCIIStreamLineReader>(path);
}
}
}