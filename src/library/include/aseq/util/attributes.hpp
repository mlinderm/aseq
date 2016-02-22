//
// Created by Michael Linderman on 12/22/15.
//

#pragma once

#include <iostream>
#include <vector>
#include <map>

#include "aseq/util/flyweight.hpp"
#include "aseq/util/any.hpp"

namespace aseq {
namespace util {

class Attributes {
  // We expect to have a small set of attributes keys overall, so use flyweight to optimize
  // storage
  struct attributes_tag {};
  typedef std::map<flyweight_string_no_track<attributes_tag>, any> AttrMap;

 public:
  typedef AttrMap::key_type key_type;
  typedef AttrMap::mapped_type mapped_type;
  typedef AttrMap::value_type value_type;
  typedef AttrMap::iterator iterator;
  typedef AttrMap::const_iterator const_iterator;

  // Common attribute types (for convenience)
  typedef int32_t Integer;
  typedef std::vector<Integer> Integers;
  typedef float Float;
  typedef std::vector<Float> Floats;
  typedef char Character;
  typedef std::vector<Character> Characters;
  typedef std::string String;
  typedef std::vector<String> Strings;

  const_iterator find(const key_type& k) const { return attrs_.find(k); }

  const_iterator begin() const { return attrs_.begin(); }
  const_iterator end() const { return attrs_.end(); }

  template <typename T>
  const T& at(const key_type& k) const {
    return any_cast<const T&>(attrs_.at(k));
  }

  template <typename T>
  T& at(const key_type& k) {
    return any_cast<T&>(attrs_.at(k));
  }

  template <typename T>
  const T& at_or(const key_type& k, const T& v) const {
    auto i = find(k);
    return i != end() ? any_cast<const T&>(i->second) : v;
  }

  mapped_type& operator[](const key_type& k) { return attrs_[k]; }

  size_t erase(const key_type& k) { return attrs_.erase(k); }

  template <class T>
  std::pair<iterator, bool> emplace(const key_type& k, const T& v) {
    return attrs_.emplace(k, mapped_type(v));
  }

  template <class T>
  std::pair<iterator, bool> emplace(const key_type& k, T&& v) {
    return attrs_.emplace(k, mapped_type(std::forward<T>(v)));
  }

 private:
  AttrMap attrs_;
};

template <>
inline const Attributes::mapped_type& Attributes::at<Attributes::mapped_type>(
    const Attributes::key_type& k) const {
  return attrs_.at(k);
}

template <>
inline std::pair<Attributes::iterator, bool> Attributes::emplace<Attributes::mapped_type>(
    const key_type& k, const mapped_type& v) {
  return attrs_.emplace(k, v);
}

template <>
inline std::pair<Attributes::iterator, bool> Attributes::emplace<Attributes::mapped_type>(
    const key_type& k, mapped_type&& v) {
  return attrs_.emplace(k, std::forward<mapped_type>(v));
}

class HasAttributes {
 public:
  HasAttributes() {}
  HasAttributes(Attributes&& attrs) : attrs_(std::move(attrs)) {}

  const Attributes& attributes() const { return attrs_; }

  bool HasAttribute(const Attributes::key_type& key) const {
    return attrs_.find(key) != attrs_.end();
  }

  template <typename T>
  const T& GetAttribute(const Attributes::key_type& key) const {
    return attrs_.at<T>(key);
  }

  template <typename T>
  T& GetAttribute(const Attributes::key_type& key) {
    return attrs_.at<T>(key);
  }

  template <typename T>
  T& GetOrAddAttribute(const Attributes::key_type& key, const T& val) {
    Attributes::iterator iter;
    std::tie(iter, std::ignore) = attrs_.emplace(key, val);
    return any_cast<T&>(iter->second);
  }

  template <typename T>
  const T& SetAttribute(const Attributes::key_type& key, T&& v) {
    return any_cast<const T&>(attrs_[key] = std::forward<T>(v));
  }

  void EraseAttribute(const Attributes::key_type& key) { attrs_.erase(key); }

 protected:
  Attributes attrs_;
};

}  // namespace util
}  // namespace aseq