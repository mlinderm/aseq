//
// Created by Michael Linderman on 2/28/16.
//

#pragma once

#include <iosfwd>
#include <string>
#include <memory>

namespace boost {
namespace filesystem {

class path;
}
}

namespace aseq {
namespace io {
class FastaSink {
 public:
  FastaSink(std::ostream& ostream) : FastaSink(ostream, 80) {}
  FastaSink(const boost::filesystem::path& path) : FastaSink(path, 80) {}
  FastaSink(std::ostream& ostream, size_t columns) : ostream_(ostream), columns_(columns) {}
  FastaSink(const boost::filesystem::path& path, size_t columns);

  FastaSink& PushSequence(const std::string& name, const std::string& sequence);

 private:
  std::unique_ptr<std::ostream> owned_ostream_;
  std::ostream& ostream_;
  size_t columns_;
};
}
}
