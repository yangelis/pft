// Copyright 2020 Ioannis Angelis <john_agelis@hotmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#ifndef __PFT_H_
#define __PFT_H_

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring> // for memset
#include <iostream>
#include <limits>
#include <vector>
#include <tuple>

#ifdef PFT_USE_ROOT
#include <TLorentzVector.h>
#include <TMath.h>
#endif

namespace pft {
//////////////////////////////////////////////////
// Particle struct
//////////////////////////////////////////////////
struct Particle {
  int pdg_id;
#ifdef PFT_USE_ROOT
  TLorentzVector vec4;
  Particle() : pdg_id(0) { vec4.SetPxPyPzE(0, 0, 0, 0); } // default constructor

  Particle(int _id, double _px, double _py, double _pz, double _e)
      : pdg_id(_id) {
    vec4.SetPxPyPzE(_px, _py, _pz, _e);
  }

  Particle(const int &_id, const TLorentzVector &p) : pdg_id(_id), vec4(p) {}

  Particle(const Particle &x) : pdg_id(x.pdg_id), vec4(x.vec4) {}

  double DeltaPhi(const Particle &p) const { return vec4.DeltaPhi(p.vec4); }
  double DeltaPhi(const TLorentzVector &v) const { return vec4.DeltaPhi(v); }

  double DeltaR(const Particle &p) const { return vec4.DeltaR(p.vec4); }

  void clear() {
    pdg_id = 0;
    vec4 = {0, 0, 0, 0};
  }

  Particle operator+(const Particle &p) const {
    return Particle(0, vec4 + p.vec4);
  }

  Particle operator+(const TLorentzVector &v) const {
    return Particle(0, vec4 + v);
  }

  bool operator==(const Particle &p) const {
    return (pdg_id == p.pdg_id && vec4 == p.vec4);
  }

  bool operator!=(const Particle &p) const {
    return (pdg_id != p.pdg_id || vec4 != p.vec4);
  }

  ~Particle() {
    pdg_id = 0;
    vec4.Clear();
  }
#else
  // px, py, pz momenta and energy
  double px, py, pz, e;

  Particle()
      : pdg_id(0), px(0.0), py(0.0), pz(0.0), e(0.0) {} // default constructor

  Particle(int _id, double _px, double _py, double _pz, double _e)
      : pdg_id(_id), px(_px), py(_py), pz(_pz), e(_e) {}
#endif
};

//////////////////////////////////////////////////
// Maybe
//////////////////////////////////////////////////
template <typename T>
struct Maybe {
  bool has_value;
  T unwrap;

  bool operator!=(const Maybe<T> &that) const { return !(*this == that); }

  bool operator==(const Maybe<T> &that) const {
    if (this->has_value && that.has_value) {
      return this->unwrap == that.unwrap;
    }

    return !this->has_value && !that.has_value;
  }
};

//////////////////////////////////////////////////
// StringView
//////////////////////////////////////////////////
struct StringView {
  size_t count{0};
  const char *data{nullptr};

  void chop(size_t n) {
    if (n > count) {
      data += count;
      count = 0;
    } else {
      data += n;
      count -= n;
    }
  }

  StringView chop_by_delim(char delim) {
    assert(data);

    size_t i = 0;
    while (i < count && data[i] != delim)
      i++;
    StringView result = {i, data};
    chop(i + 1);

    return result;
  }
};

StringView operator""_sv(const char *data, size_t count) {
  return {count, data};
}

std::vector<std::string> split_by(StringView &sv, char delim) {
  std::vector<std::string> vec;
  StringView aug = {};
  while (0 < sv.count) {
    aug = sv.chop_by_delim(delim);
    vec.emplace_back(aug.data, aug.count);
  }
  return vec;
}

StringView cstr_as_sv(const char *cstr) { return {strlen(cstr), cstr}; }

StringView string_as_sv(const std::string &s) { return {s.length(), s.data()}; }

Maybe<StringView> read_file_as_string_view(const char *filename) {
  FILE *f = fopen(filename, "rb");
  if (!f)
    return {};

  int err = fseek(f, 0, SEEK_END);
  if (err < 0)
    return {};

  long size = ftell(f);
  if (size < 0)
    return {};

  err = fseek(f, 0, SEEK_SET);
  if (err < 0)
    return {};

  auto data = malloc(size);
  if (!data)
    return {};

  size_t read_size = fread(data, 1, size, f);
  if (read_size != (size_t)size && ferror(f))
    return {};

  fclose(f);
  return {true, {static_cast<size_t>(size), static_cast<const char *>(data)}};
}

void ignore_header_lines(std::vector<std::string> &vec, int lines) {
  vec.erase(vec.begin(), vec.begin() + lines);
}

std::vector<float> as_float(const std::vector<std::string> &vec) {
  std::vector<float> buffer(vec.size());

  for (size_t i = 0; i < vec.size(); ++i) {
    buffer[i] = std::stof(vec[i]);
  }
  return buffer;
}

void print1(FILE *stream, StringView view) {
  fwrite(view.data, 1, view.count, stream);
}

//////////////////////////////////////////////////
// Matrix
//////////////////////////////////////////////////

template <typename T>
struct Matrix1v {
  size_t rows, cols;
  std::vector<T> data;
  Matrix1v() : rows(0), cols(0), data(rows * cols, 0) {}
  Matrix1v(size_t r, size_t c) : rows(r), cols(c), data(rows * cols, 0) {}

  T &operator()(size_t i, size_t j) { return data.at(j + i * cols); }

  const T &operator()(size_t i, size_t j) const {
    return data.at(j + i * cols);
  }

  T trace() {
    T tr = 0;
    if (cols == rows) {
      for (int i = 0; i < cols; i++) {
        tr += ((*this)(i, i));
      }
      return tr;
    } else {
      return std::numeric_limits<T>::quiet_NaN();
    }
  }

  void transpose() {
    std::vector<T> tempv;
    for (int i = 0; i < cols; i++) {
      for (int j = 0; j < rows; j++) {
        tempv.push_back((*this)(j, i));
      }
    }
    std::swap(cols, rows);
    data.clear();
    data = tempv;
  }
};

template <typename T>
struct Matrix2v {
  size_t rows, cols;
  T *data;
  Matrix2v(size_t r, size_t c) : rows(c), cols(c) {
    data = (T *)malloc(rows * cols * sizeof(T));
    memset(data, 0, (rows * cols) * sizeof(*data));
  }

  T &operator()(size_t i, size_t j) { return data[j + i * cols]; }

  const T &operator()(size_t i, size_t j) const { return data[j + i * cols]; }
};


//////////////////////////////////////////////////
// Enumeration
//////////////////////////////////////////////////
// http://reedbeta.com/blog/python-like-enumerate-in-cpp17/
template <typename T, typename TIter = decltype(std::begin(std::declval<T>())),
          typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T &&iterable) {
  struct iterator {
    size_t i;
    TIter iter;
    bool operator!=(const iterator &other) const { return iter != other.iter; }
    void operator++() {
      ++i;
      ++iter;
    }
    auto operator*() const { return std::tie(i, *iter); }
  };
  struct iterable_wrapper {
    T iterable;
    auto begin() { return iterator{0, std::begin(iterable)}; }
    auto end() { return iterator{0, std::end(iterable)}; }
  };
  return iterable_wrapper{std::forward<T>(iterable)};
}

} // namespace pft

#endif // __PFT_H_
