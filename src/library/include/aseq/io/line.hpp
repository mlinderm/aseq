//
// Created by Michael Linderman on 12/12/15.
//

#pragma once

#include <memory>
#include <iosfwd>

#include <boost/optional.hpp>
#include <boost/range/iterator_range.hpp>

#include "aseq/util/exception.hpp"
#include "aseq/io/file_format.hpp"
#include "aseq/model/region.hpp"

namespace boost {
namespace filesystem {

class path;
}
}

namespace aseq {
namespace io {

typedef boost::iterator_range<const char*> Line;

class ASCIILineReaderInterface {
 public:
  typedef std::unique_ptr<ASCIILineReaderInterface> FactoryResult;
  typedef boost::optional<Line> NextResult;

 public:
  virtual ~ASCIILineReaderInterface() {}

  virtual bool IsIndexed() const { return false; }
  virtual void SetRegion(model::Contig contig, model::Pos pos, model::Pos end) {
    throw util::indexed_access_not_supported();
  }
  virtual NextResult ReadNextLine() = 0;

  static FactoryResult MakeLineReader(std::istream& istream);
  static FactoryResult MakeLineReader(const boost::filesystem::path& file);
};

class ASCIILineWriterInterface {
 public:
  typedef std::unique_ptr<ASCIILineWriterInterface> FactoryResult;

 public:
  virtual ~ASCIILineWriterInterface() {}

  virtual void Write(const Line&) = 0;
  void Write(const char* line) { Write(boost::make_iterator_range(line, line + strlen(line))); }
  void Write(const std::string& line) {
    Write(boost::make_iterator_range(line.data(), line.data() + line.size()));
  }

  static FactoryResult MakeLineWriter(std::ostream& ostream);
  static FactoryResult MakeLineWriter(const boost::filesystem::path& file,
                                      FileFormat format = FileFormat::UNKNOWN);
};

}  // namespace io
}  // namespace aseq