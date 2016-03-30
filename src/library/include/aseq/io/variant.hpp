//
// Created by Michael Linderman on 12/13/15.
//

#pragma once

#include <memory>

#include <boost/optional.hpp>

#include "aseq/io/file_format.hpp"
#include "aseq/io/line.hpp"
#include "aseq/model/variant_context.hpp"

namespace boost {
namespace filesystem {
class path;
}
}

namespace aseq {
namespace io {

class VariantHeaderInterface {
 public:
  virtual size_t NumSamples() const = 0;
  virtual const model::Sample& Sample(size_t idx) const = 0;

  virtual void SetSitesOnly() = 0;
};

class VariantSourceInterface {
 public:
  typedef std::unique_ptr<VariantSourceInterface> FactoryResult;
  typedef boost::optional<model::VariantContext> NextResult;

 public:
  virtual ~VariantSourceInterface() {}

  virtual FileFormat file_format() const = 0;
  virtual const VariantHeaderInterface& header() const = 0;

  virtual NextResult NextVariant() = 0;

  static FactoryResult MakeVariantSource(std::istream& istream);
  static FactoryResult MakeVariantSource(const boost::filesystem::path& path);
};

class VariantSinkInterface {
 public:
  typedef std::unique_ptr<VariantSinkInterface> FactoryResult;

 public:
  virtual ~VariantSinkInterface() {}

  virtual void PushVariant(const model::VariantContext&) = 0;

  static FactoryResult MakeVariantSink(FileFormat format, std::ostream& ostream);
  static FactoryResult MakeVariantSink(FileFormat format, const boost::filesystem::path& path);
  static FactoryResult MakeVariantSink(const VariantSourceInterface& source, std::ostream& ostream,
                                       bool sites_only = false);
  static FactoryResult MakeVariantSink(const VariantSourceInterface& source,
                                       const boost::filesystem::path& path,
                                       bool sites_only = false);
};
}
}
