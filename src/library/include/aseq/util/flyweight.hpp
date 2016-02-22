//
// Created by Michael Linderman on 12/16/15.
//

#pragma once

#include <string>
#include <functional>

#include <boost/functional/hash.hpp>
#include <boost/flyweight.hpp>
#include <boost/flyweight/tag.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <boost/range/iterator_range.hpp>

namespace aseq {
namespace util {

template <class Tag>
class flyweight_string_no_track {
  typedef typename boost::flyweight<std::string, boost::flyweights::tag<Tag>,
                                    boost::flyweights::no_tracking> flyweight_type;

 public:
  typedef typename flyweight_type::initializer initializer;

  flyweight_string_no_track() {}
  flyweight_string_no_track(char c) : value_(std::to_string(c)) {}
  flyweight_string_no_track(const std::string& s) : value_(s) {}
  flyweight_string_no_track(const char* s) : value_(s) {}
  template <typename I>
  flyweight_string_no_track(const boost::iterator_range<I>& r)
      : value_(r.begin(), r.end()) {}
  template <typename I>
  flyweight_string_no_track(I begin, I end)
      : value_(begin, end) {}
  flyweight_string_no_track(const flyweight_string_no_track&) = default;
  flyweight_string_no_track(flyweight_string_no_track&&) = default;

  operator const std::string&() const { return value_.get(); }
  const std::string& get() const { return value_.get(); }

  /* String Interface */
  typedef std::string::value_type value_type;
  typedef std::string::const_reference const_reference;
  typedef std::string::const_iterator const_iterator;
  typedef std::string::const_reverse_iterator reverse_const_iterator;

  const_iterator begin() const { return this->get().begin(); }
  const_iterator end() const { return this->get().end(); }
  reverse_const_iterator rbegin() const { return this->get().rbegin(); }
  reverse_const_iterator rend() const { return this->get().rend(); }

  size_t size() const { return this->get().size(); }
  bool empty() const { return this->get().empty(); }

  const_reference front() const { return this->get().front(); }
  const_reference back() const { return this->get().back(); }

  const char* c_str() const { return this->get().c_str(); }
  flyweight_string_no_track substr(size_t pos = 0, size_t len = std::string::npos) const {
    return flyweight_string_no_track(get().substr(pos, len));
  }

  /* Operators */
  flyweight_string_no_track& operator=(const flyweight_string_no_track& f) {
    value_ = f.value_;
    return *this;
  }
  flyweight_string_no_track& operator=(flyweight_string_no_track&& f) {
    value_ = f.value_;
    return *this;
  }
  bool operator==(const flyweight_string_no_track& f) const { return value_ == f.value_; }
  bool operator!=(const flyweight_string_no_track& f) const { return value_ != f.value_; }
  bool operator<(const flyweight_string_no_track& f) const { return value_ < f.value_; }

 private:
  flyweight_type value_;
};

template <class Tag>
std::ostream& operator<<(std::ostream& ostream, const flyweight_string_no_track<Tag>& flyweight) {
  return (ostream << flyweight.get());
}
}  // namespace util
}  // namespace aseq

namespace std {

template <class T>
struct hash<aseq::util::flyweight_string_no_track<T> > {
  std::size_t operator()(const aseq::util::flyweight_string_no_track<T>& k) const {
    std::hash<const void*> hasher;
    return hasher(&k.get());
  }
};
}