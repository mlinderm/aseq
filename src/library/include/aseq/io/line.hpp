//
// Created by Michael Linderman on 12/12/15.
//

#pragma once

#include <memory>
#include <iosfwd>

#include <boost/optional.hpp>
#include <boost/range/iterator_range.hpp>

namespace boost {
namespace filesystem {

class path;
}
}

namespace aseq {
namespace io {

class ASCIILineReaderInterface {
 public:
  typedef boost::iterator_range<const char*> Line;
  typedef std::unique_ptr<ASCIILineReaderInterface> FactoryResult;
  typedef boost::optional<Line> NextResult;

 public:
  virtual ~ASCIILineReaderInterface() {}

  virtual NextResult ReadNextLine() = 0;

  static FactoryResult MakeLineReader(std::istream& istream);
  static FactoryResult MakeLineReader(const boost::filesystem::path& file);
};

class ASCIILineWriterInterface {
 public:
  typedef boost::iterator_range<const char*> Line;
  typedef std::unique_ptr<ASCIILineWriterInterface> FactoryResult;

 public:
  virtual ~ASCIILineWriterInterface() {}

  virtual void WriteLine(const Line&) = 0;
  void WriteLine(const char* line) {
    WriteLine(boost::make_iterator_range(line, line + strlen(line)));
  }
  void WriteLine(const std::string& line) {
    WriteLine(boost::make_iterator_range(line.data(), line.data() + line.size()));
  }

  static FactoryResult MakeLineWriter(std::ostream& ostream);
  static FactoryResult MakeLineWriter(const boost::filesystem::path& file);
};

}  // namespace io
}  // namespace aseq