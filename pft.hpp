// Copyright 2021 Ioannis Angelis <john_agelis@hotmail.com>
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
// ============================================================
//
// ChangeLog:
//   0.0.9    print1<complex>, to_f64, contains,
//            max_id_and_val, move_if, adjacent_transformN,
//            adjacent_findN, min_id_and_val
//            remove AParse
//   0.0.8    pad_right_until
//   0.0.7    linspace, pad_left, pad_right
//            zip_with, zip_to_pair
//            remove zip,
//            GetVectorsSize, map
//            takeFromIdx, Slice, chunks
//            Some
//   0.0.6    mod, findfirst, findall,
//            panic, unwrap_or_panic
//   0.0.5    map, take, filter, var,
//            arange, normalize, abs,
//            zip, pad<vector>
//            operators *,<, > for vector
//            where, argsort
//   0.0.4    remove Particle,
//            println for std::vector<T>
//   0.0.3    vec_from_range
//   0.0.2    pop, drop, readlines
//   0.0.1    Maybe<T>, StringView,
//            Particles_t, Matrix
//
// ============================================================

#ifndef PFT_H_
#define PFT_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

using c8 = char;
// Unsigned
using u8  = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

// Signed
using i8  = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

// Floating point
using f32 = float;
using f64 = double;

template <typename T>
using Vec = std::vector<T>;
//////////////////////////////////////////////////
// Operators for vector<T>
//////////////////////////////////////////////////
template <typename T>
Vec<T> operator*(const T& elem, const Vec<T>& rhs)
{
  Vec<T> res(rhs.size(), 0);
  for (std::size_t i = 0; i < rhs.size(); ++i) {
    res[i] = elem * rhs[i];
  }
  return res;
}

template <typename T>
Vec<T> operator/(const Vec<T>& rhs, const T& elem)
{
  Vec<T> res(rhs.size(), 0);
  for (std::size_t i = 0; i < rhs.size(); ++i) {
    res[i] = rhs[i] / elem;
  }
  return res;
}

template <typename T>
Vec<i32> operator<(const Vec<T>& lhs, const T& elem)
{
  Vec<i32> res(lhs.size());

  for (std::size_t i = 0; i < lhs.size(); ++i) {
    if (lhs[i] < elem) {
      res[i] = 1;
    } else {
      res[i] = 0;
    }
  }
  return res;
}

template <typename T>
Vec<i32> operator>(const T& elem, const Vec<T>& rhs)
{
  return rhs < elem;
}

template <typename T>
Vec<i32> operator>(const Vec<T>& lhs, const T& elem)
{
  Vec<i32> res(lhs.size());

  for (std::size_t i = 0; i < lhs.size(); ++i) {
    if (lhs[i] > elem) {
      res[i] = 1;
    } else {
      res[i] = 0;
    }
  }
  return res;
}

template <typename T>
Vec<i32> operator<(const T& elem, const Vec<T>& rhs)
{
  return rhs > elem;
}

template <typename T>
Vec<i32> operator==(const Vec<T>& lhs, const T& elem)
{
  Vec<i32> ret(lhs.size(), 1);

  for (std::size_t i = 0; i < lhs.size(); ++i) {
    if (lhs[i] != elem) {
      ret[i] = 0;
    }
  }
  return ret;
}

namespace pft {

//////////////////////////////////////////////////
// Slice
//////////////////////////////////////////////////
struct Slice {
  i64 start;
  i64 stop;
  constexpr auto size() const -> std::size_t
  {
    return static_cast<std::size_t>(stop - start);
  }
};

// Maybe and StringView are based on https://github.com/rexim/aids
//////////////////////////////////////////////////
// Maybe
//////////////////////////////////////////////////
template <typename T>
struct Maybe {
  bool has_value;
  T unwrap;

  constexpr Maybe() : has_value(false) {}
  constexpr Maybe(bool f, T& val) : has_value(f), unwrap(val) {}
  constexpr Maybe(bool f, T&& val) noexcept
      : has_value(f), unwrap(std::move(val))
  {
  }

  constexpr bool operator!=(const Maybe<T>& that) const
  {
    return !(*this == that);
  }

  constexpr bool operator==(const Maybe<T>& that) const
  {
    if (this->has_value && that.has_value) {
      return this->unwrap == that.unwrap;
    }

    return !this->has_value && !that.has_value;
  }

  template <typename F>
  constexpr auto bind(F&& f) const
  {
    if (has_value) {
      return Maybe(true, f(unwrap));
    }
    return Maybe();
  }
};

template <typename T>
constexpr auto Some(T m) -> Maybe<T>
{
  return {true, m};
}

//////////////////////////////////////////////////
// StringView
//////////////////////////////////////////////////
struct StringView : public std::string_view {

  StringView() : std::string_view("") {}
  StringView(const std::string& s) : std::string_view(s) {}
  StringView(const std::string_view& s) : std::string_view(s) {}
  StringView(std::string&& s) : std::string_view(std::move(s)) {}
  StringView(std::string_view&& s) : std::string_view(s) {}
  StringView(const c8* s) : std::string_view(s) {}
  StringView(const c8* s, std::size_t l) : std::string_view(s, l) {}

  auto chop(std::size_t n) -> StringView
  {
    if (this->size() > n) {
      return this->substr(n);
    } else {
      return {};
    }
  }

  auto chop_by_delim(c8 delim) -> StringView
  {
    std::size_t i = 0;
    while (i < this->size() && this->data()[i] != delim) {
      ++i;
    }
    StringView result = {this->data(), i};
    this->remove_prefix(i + 1);
    return result;
  }

  auto as_string() const -> std::string { return std::string(*this); }
};

//////////////////////////////////////////////////
// Particles struct usefull for Geant4
//////////////////////////////////////////////////
struct Particles_t {
  Vec<i32> det_id, parent_id, trid, n_secondaries;
  Vec<f64> times, edep, energy, posX, posY, posZ;
  Vec<f64> theta, phi, trlen;

  void Reserve(std::size_t nparticles)
  {
    det_id.reserve(nparticles);
    parent_id.reserve(nparticles);
    trid.reserve(nparticles);
    times.reserve(nparticles);
    edep.reserve(nparticles);
    energy.reserve(nparticles);
    posX.reserve(nparticles);
    posY.reserve(nparticles);
    posZ.reserve(nparticles);
    theta.reserve(nparticles);
    phi.reserve(nparticles);
    trlen.reserve(nparticles);
    n_secondaries.reserve(nparticles);
  }

  void ClearVecs()
  {
    det_id.clear();
    parent_id.clear();
    trid.clear();
    times.clear();
    edep.clear();
    energy.clear();
    posX.clear();
    posY.clear();
    posZ.clear();
    theta.clear();
    phi.clear();
    trlen.clear();
    n_secondaries.clear();
  }
};

// StringView utilities
static inline StringView operator""_sv(const c8* data, std::size_t count)
{
  return {data, count};
}

static inline auto trimr(StringView& s) -> StringView&
{
  auto i = s.find_last_not_of(" \t\n\r\f\v");
  if (i != std::string_view::npos) {
    s = s.substr(0, i + 1);
  }

  return s;
}

static inline auto triml(StringView& s) -> StringView&
{
  auto i = s.find_first_not_of(" \t\n\r\f\v");
  if (i != std::string_view::npos) {
    s.remove_prefix(i);
  }

  return s;
}

static inline auto split_by(StringView view, c8 delim) -> Vec<StringView>
{
  view = trimr(view);
  Vec<StringView> ret;
  while (!view.empty()) {
    auto len = view.find(delim);
    if (len == std::string_view::npos) {
      ret.emplace_back(view);
      break;
    }
    ret.emplace_back(view.substr(0, len));
    view = view.substr(len + 1);
  }
  return ret;
}

static inline auto split_by(std::string& str, c8 delim) -> Vec<StringView>
{
  Vec<StringView> vec;
  StringView temp{str};
  StringView aug = {};
  while (!temp.empty()) {
    aug = temp.chop_by_delim(delim);
    vec.emplace_back(aug.data(), aug.size());
  }
  return vec;
}

static inline auto cstr_as_sv(const c8* cstr) -> StringView
{
  return {cstr, strlen(cstr)};
}

static inline auto string_as_sv(const std::string& s) -> StringView
{
  return {s.data(), s.length()};
}

static inline auto to_i32(StringView s) -> i32
{
  return std::stoi(std::string(s));
}

static inline auto to_f64(pft::StringView s) -> f64
{
  return std::stod(std::string(s));
}

static inline auto read_file_as_string_view(const c8* filename)
    -> Maybe<StringView>
{
  FILE* f = fopen(filename, "rb");
  if (f == nullptr) {
    return {};
  }

  i32 err = fseek(f, 0, SEEK_END);
  if (err < 0) {
    return {};
  }

  long size = ftell(f);
  if (size < 0) {
    return {};
  }

  err = fseek(f, 0, SEEK_SET);
  if (err < 0) {
    return {};
  }

  auto data = malloc(size);
  if (data == nullptr) {
    return {};
  }

  std::size_t read_size = fread(data, 1, size, f);
  if (read_size != static_cast<std::size_t>(size) && ferror(f) != 0) {
    return {};
  }

  fclose(f);
  return Some(
      StringView{static_cast<const c8*>(data), static_cast<std::size_t>(size)});
}

static inline auto readlines(const c8* filename, const c8 delim = '\n')
    -> Vec<StringView>
{

  auto file   = read_file_as_string_view(filename);
  auto result = split_by(file.unwrap, delim);
  return result;
}

static inline void ignore_header_lines(Vec<StringView>& vec, i32 lines)
{
  vec.erase(std::begin(vec), std::begin(vec) + lines);
}

static inline auto as_floats(const Vec<StringView>& vec) -> Vec<f32>
{
  Vec<f32> buffer(vec.size());

  for (std::size_t i = 0; i < vec.size(); ++i) {
    buffer[i] = std::stof(std::string(vec[i]));
  }
  return buffer;
}

//////////////////////////////////////////////////
// Matrix
//////////////////////////////////////////////////
template <typename T>
struct Matrix {
  std::size_t rows, cols;
  T* data{nullptr};

  constexpr Matrix() : rows(0), cols(0), data(nullptr) {}
  constexpr Matrix(std::size_t r, std::size_t c) : rows(r), cols(c)
  {
    data = new T[r * c]();
  }
  constexpr Matrix(const Matrix<T>& m) : rows(m.rows), cols(m.cols)
  {
    data = new T[rows * cols]();
    std::memcpy(data, m.data, sizeof(T) * rows * cols);
  }

  template <std::size_t r, std::size_t c>
  constexpr Matrix(const T (&m)[r][c]) : Matrix(r, c)
  {
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        (*this)(i, j) = m[i][j];
      }
    }
  }

  ~Matrix() { delete[] this->data; }

  constexpr T& operator()(std::size_t i, std::size_t j)
  {
    return data[j + i * cols];
  }
  constexpr T operator()(std::size_t i, std::size_t j) const
  {
    return data[j + i * cols];
  }

  constexpr Matrix<T>& operator=(const Matrix<T>& a)
  {
    if (this == &a) {
      return *this;
    }
    if (data != nullptr) {
      delete[] data;
    }
    rows       = a.rows;
    cols       = a.cols;
    this->data = new T[rows * cols]();
    std::memcpy(data, a.data, sizeof(T) * rows * cols);
    return *this;
  }

  constexpr void diagonal(T d = static_cast<T>(1.0))
  {
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        if (i == j) {
          (*this)(i, j) = d;
        }
      }
    }
  }

  constexpr auto trace() const -> T
  {
    T tr = 0;
    if (cols == rows) {
      for (std::size_t i = 0; i < cols; ++i) {
        tr += ((*this)(i, i));
      }
      return tr;
    } else {
      return std::numeric_limits<T>::quiet_NaN();
    }
  }

  constexpr auto transpose() const -> Matrix<T>
  {
    Matrix<T> result(cols, rows);
    for (std::size_t i = 0; i < cols; ++i) {
      for (std::size_t j = 0; j < rows; ++j) {
        result(i, j) = ((*this)(j, i));
      }
    }
    return result;
  }

  constexpr auto mult(const Matrix<T>& b) -> Matrix<T>
  {
    if (cols != b.rows) {
      exit(1);
    }
    Matrix<T> ret(rows, b.cols);
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < b.cols; ++j) {
        for (std::size_t k = 0; k < cols; ++k) {
          ret(i, j) += ((*this)(i, k)) * b(k, j);
        }
      }
    }
    return ret;
  }

  constexpr auto getColumn(std::size_t c) const -> Vec<T>
  {
    Vec<T> ret(rows);
    for (std::size_t i = 0; i < rows; ++i) {
      ret[i] = (*this)(i, c);
    }
    return ret;
  }

  constexpr auto GetMinor(std::size_t d) const -> Matrix<T>
  {
    Matrix<T> ret(rows, cols);
    for (std::size_t i = 0; i < d; ++i) {
      ret(i, i) = 1;
    }

    for (std::size_t i = d; i < rows; ++i) {
      for (std::size_t j = d; j < cols; ++j) {
        ret(i, j) = (*this)(i, j);
      }
    }

    return ret;
  }
};

//////////////////////////////////////////////////
// Printers
//////////////////////////////////////////////////
static inline void print1(FILE* stream, c8 x) { fputc(x, stream); }
static inline void print1(FILE* stream, c8* x)
{
  fwrite(x, 1, std::strlen(x), stream);
}
static inline void print1(FILE* stream, const c8* x)
{
  fwrite(x, 1, std::strlen(x), stream);
}

static inline void print1(FILE* stream, u32 x) { fprintf(stream, "%u", x); }
static inline void print1(FILE* stream, u64 x) { fprintf(stream, "%lu", x); }
// static inline void print1(FILE* stream, std::size_t x) { fprintf(stream, "%zu", x); }

static inline void print1(FILE* stream, i16 x) { fprintf(stream, "%hd", x); }
static inline void print1(FILE* stream, i32 x) { fprintf(stream, "%d", x); }
static inline void print1(FILE* stream, i64 x) { fprintf(stream, "%ld", x); }

static inline void print1(FILE* stream, f32 f) { fprintf(stream, "%8.4f", f); }
static inline void print1(FILE* stream, f64 f) { fprintf(stream, "%8.4f", f); }

static inline void print1(FILE* stream, StringView view)
{
  fwrite(view.data(), 1, view.size(), stream);
}

static inline void print1(FILE* stream, std::complex<f64> f)
{
  fprintf(stream, "{%8.4f, %8.4f}", f.real(), f.imag());
}
// NOTE: forward declare to use Maybe<Vec>
template <typename T>
static inline void print1(FILE* stream, const Vec<T>& v);

template <typename T>
static inline void print1(FILE* stream, const Maybe<T>& m)
{
  print1(stream, "Maybe{ ");
  print1(stream, m.has_value);
  print1(stream, ", ");
  print1(stream, m.unwrap);
  print1(stream, " }");
}

template <typename T, typename U>
static inline void print1(FILE* stream, const std::pair<T, U>& pr)
{
  print1(stream, "(");
  print1(stream, pr.first);
  print1(stream, ", ");
  print1(stream, pr.second);
  print1(stream, ")");
}

template <typename T>
static inline void print1(FILE* stream, const Vec<T>& v)
{
  if (v.empty()) {
    return;
  }
  const std::size_t n = v.size();
  fprintf(stream, "vector(size=%zu) ", v.size());
  print1(stream, "{");
  for (std::size_t i = 0; i < n - 1; ++i) {
    print1(stream, v[i]);
    print1(stream, ", ");
  }
  print1(stream, v[n - 1]);
  print1(stream, "}");
}

// TODO: make pretty printer for matrices
template <typename T>
static inline void print1(FILE* stream, const Matrix<T>& mat)
{
  for (std::size_t i = 0; i < mat.rows; ++i) {
    for (std::size_t j = 0; j < mat.cols; ++j) {
      print1(stream, mat(i, j));
      print1(stream, "  ");
    }
    fprintf(stream, "\n");
  }
}

template <typename... Types>
static inline void print(FILE* stream, Types... args)
{
  (print1(stream, args), ...);
}

static inline void println(FILE* stream) { print1(stream, '\n'); }

template <typename... Types>
static inline void println(FILE* stream, Types... args)
{
  print(stream, args...);
  print1(stream, '\n');
}

// stolen from https://github.com/rexim/aids/blob/master/aids.hpp
template <typename... Args>
[[noreturn]] void panic(Args... args)
{
  println(stderr, args...);
  exit(1);
}

// stolen from https://github.com/rexim/aids/blob/master/aids.hpp
template <typename T, typename... Args>
static inline auto unwrap_or_panic(Maybe<T> maybe, Args... args) -> T
{
  if (!maybe.has_value) {
    panic(args...);
  }

  return maybe.unwrap;
}

//////////////////////////////////////////////////
// Utils
//////////////////////////////////////////////////
template <class InputIt, class OutputIt, class UnaryPred,
          typename = typename std::enable_if_t<
              std::is_convertible_v<
                  typename std::iterator_traits<InputIt>::iterator_category,
                  std::input_iterator_tag> &&
              std::is_same_v<
                  typename std::iterator_traits<OutputIt>::iterator_category,
                  std::output_iterator_tag>>>
constexpr static inline auto move_if(InputIt first, InputIt last,
                                     OutputIt result, UnaryPred pred)
    -> OutputIt
{
  // auto left = last;
  for (; first != last; ++first) {
    if (pred(*first)) {
      *result = std::move(*first);
      // std::swap(*first, *(--left));
      ++result;
    }
  }
  return result;
}

template <class Iter, std::size_t N, std::size_t... Is>
constexpr static inline auto unpack_iters_impl(std::array<Iter, N>& arr,
                                               std::index_sequence<Is...>)
    -> decltype(std::tie(*arr[Is]...))
{
  return std::tie(*arr[Is]...);
}

template <std::size_t N, class ForwardIt, class Predicate,
          typename = typename std::enable_if_t<std::is_convertible_v<
              typename std::iterator_traits<ForwardIt>::iterator_category,
              std::forward_iterator_tag>>>
constexpr static inline auto adjacent_findN(ForwardIt first, ForwardIt last,
                                            Predicate pred) -> ForwardIt
{
  if (first != last) {
    std::array<ForwardIt, N> iters;
    for (std::size_t i = 0; i < N; ++i) {
      iters[i] = std::next(first, i);
    }

    for (; iters.back() != last;
         std::for_each(std::rbegin(iters), std::rend(iters),
                       [](auto& x) { std::advance(x, 1); })) {
      if (std::apply(std::forward<Predicate>(pred),
                     unpack_iters_impl(iters, std::make_index_sequence<N>{}))) {
        return iters[0];
      }
    }
  }
  return last;
}

template <std::size_t N, class InputIt, class OutputIt, class Op,
          typename = typename std::enable_if_t<
              std::is_convertible_v<
                  typename std::iterator_traits<InputIt>::iterator_category,
                  std::input_iterator_tag> &&
              std::is_same_v<
                  typename std::iterator_traits<OutputIt>::iterator_category,
                  std::output_iterator_tag>>>
constexpr static inline auto adjacent_transformN(InputIt first, InputIt last,
                                                 OutputIt result, Op op)
    -> OutputIt
{
  if (first != last) {
    using val_t = typename std::iterator_traits<InputIt>::value_type;
    std::array<InputIt, N> iters;
    for (std::size_t i = 0; i < N; ++i) {
      iters[i] = std::next(first, i + 1);
    }
    val_t acc(*first);
    *result = acc;

    for ((void)++result; iters[N - 2] != last;
         std::for_each(std::rbegin(iters), std::rend(iters),
                       [](auto& x) { std::advance(x, 1); }),
         (void)++result) {
      const auto params = std::tuple_cat(
          unpack_iters_impl(iters, std::make_index_sequence<N - 1>{}),
          std::tie(acc));
      val_t val(*iters[0]);

      *result = std::apply(op, params);

      acc = std::move(val);
    }
  }
  return result;
}

// http://reedbeta.com/blog/python-like-enumerate-in-cpp17/
template <typename T, typename TIter = decltype(std::begin(std::declval<T>())),
          typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T&& iterable)
{
  struct iterator {
    std::size_t i;
    TIter iter;
    bool operator!=(const iterator& other) const { return iter != other.iter; }
    void operator++()
    {
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

template <typename F, typename T, typename U,
          typename R = Vec<decltype(std::declval<F>()(std::declval<T>(),
                                                      std::declval<U>()))>>
static inline auto zip_with(F&& fn, const Vec<T>& a, const Vec<U>& b) -> R
{
  const std::size_t n = std::min(a.size(), b.size());
  R ret;
  ret.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    ret.emplace_back(fn(a[i], b[i]));
  }
  return ret;
}

template <typename T, typename U>
static inline auto zip_to_pair(const Vec<T>& a, const Vec<U>& b)
{
  auto MakePair = [](const T& x, const U& y) { return std::make_pair(x, y); };

  return zip_with(MakePair, a, b);
}

template <typename T, typename R, typename FoldOp>
static inline auto foldl(const R& i, const Vec<T>& xs, FoldOp&& fn) -> R
{
  auto ret = i;
  for (const auto& x : xs) {
    ret = fn(ret, x);
  }
  return ret;
}

template <typename T>
static inline auto sum(const Vec<T>& xs) -> T
{
  return std::accumulate(std::cbegin(xs), std::cend(xs), T());
  // return foldl(T(), xs, [](const T& a, const T& b) -> T { return a + b; });
}

template <typename T>
static inline auto mean(const Vec<T>& xs) -> f64
{
  return f64(sum(xs) / (xs.size()));
}

// using the two pass formula
template <typename T>
static inline auto var(const Vec<T>& xs) -> f64
{
  const std::size_t n = xs.size();
  if (n < std::size_t(2)) {
    fprintf(stderr, "Need atleast 2 elements for variance\n");
    return 0;
  }
  f64 sum_squares = 0.0, squared_sum = 0.0;
  const auto xs_mean = mean(xs);
  auto pred          = [&sum_squares, &squared_sum, &xs_mean](const T& x) {
    squared_sum += x - xs_mean;
    sum_squares += (x - xs_mean) * (x - xs_mean);
  };
  std::for_each(std::cbegin(xs), std::cend(xs), pred);
  const auto N = static_cast<f64>(n);
  return (sum_squares - squared_sum * squared_sum / N) / (N - 1);
}

template <typename T>
static inline auto stdev(const Vec<T>& xs) -> T
{
  return std::sqrt(var(xs));
}

template <typename... Types>
static inline auto GetVectorsSize(const Vec<Types>&... vs) -> std::size_t
{
  constexpr const auto nArgs = sizeof...(Types);
  const std::size_t sizes[]  = {vs.size()...};
  if (nArgs > 1) {
    for (std::size_t i = 1; i < nArgs; ++i) {
      if (sizes[0] == sizes[i]) {
        continue;
      }
      pft::panic("Vectors have different lengths");
    }
  }
  return sizes[0];
}

template <typename F, typename T,
          typename R = decltype(std::declval<F>()(std::declval<T>()))>
static inline auto map(F&& fn, const Vec<T>& input) -> Vec<R>
{
  const auto size = input.size();
  Vec<R> ret;
  ret.reserve(size);
  std::transform(std::cbegin(input), std::cend(input), std::back_inserter(ret),
                 fn);
  return ret;
}

template <typename F, typename T,
          typename R = decltype(std::declval<F>()(std::declval<T>()))>
static inline auto map(F&& fn, Vec<T>&& input) -> Vec<R>
{
  Vec<R> ret;
  ret.reserve(input.size());
  std::transform(std::make_move_iterator(std::begin(input)),
                 std::make_move_iterator(std::end(input)),
                 std::back_inserter(ret), std::forward<F>(fn));
  return ret;
}

template <typename F, typename... Types,
          typename R = decltype(std::declval<F>()(std::declval<Types>()...))>
static inline auto map(F&& fn, const Vec<Types>&... input) -> Vec<R>
{
  const auto size = GetVectorsSize(input...);
  Vec<R> ret;
  ret.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    ret.emplace_back(fn(input[i]...));
  }
  return ret;
}

template <typename F, typename T>
static inline auto filter(F&& fn, const Vec<T>& v) -> Vec<T>
{
  const auto n = v.size();
  Vec<T> ret;
  ret.reserve(n);
  for (auto&& val : v) {
    if (fn(val)) {
      ret.emplace_back(val);
    }
  }
  return ret;
}

template <typename F, typename T>
static inline auto filter(F&& fn, Vec<T>&& v) -> Vec<T>
{
  Vec<T> ret;
  ret.reserve(v.size());
  move_if(std::begin(v), std::end(v), std::back_inserter(ret),
          std::forward<F>(fn));
  (void)std::exchange(v, {}); // "destroy" old vec
  return ret;
}

template <typename T>
static inline void pop(Vec<T>& vec, std::size_t elements)
{
  vec.erase(std::begin(vec), std::begin(vec) + elements);
}

template <typename T>
static inline void drop(Vec<T>& vec, std::size_t elements)
{
  for (std::size_t i = 0; i < elements; ++i) {
    vec.pop_back();
  }
}

template <typename T>
static inline auto take(const Vec<T>& vec, std::size_t elements) -> Vec<T>
{
  return {std::cbegin(vec), std::cbegin(vec) + elements};
}

template <typename T>
static inline auto take(const Vec<T>& vec, const Slice& slice) -> Vec<T>
{
  return Vec<T>{std::cbegin(vec) + slice.start,
                std::cbegin(vec) +
                    std::min(static_cast<size_t>(slice.stop), vec.size())};
}

template <typename T>
static inline auto take(Vec<T>&& vec, const pft::Slice& slice) -> Vec<T>
{
  Vec<T> ret;
  ret.reserve(slice.size());
  std::move(std::begin(vec) + slice.start,
            std::begin(vec) +
                std::min(static_cast<size_t>(slice.stop), vec.size()),
            std::back_inserter(ret));
  return ret;
}

template <typename T>
static inline auto take(const Vec<T>& vec, const Vec<i32>& indices) -> Vec<T>
{
  const std::size_t n = indices.size();
  Vec<T> ret;
  ret.reserve(n);

  auto l = [&vec](const i32& i) { return vec[i]; };
  std::transform(std::cbegin(indices), std::cend(indices),
                 std::back_inserter(ret), l);
  return ret;
}

template <typename T>
static inline auto takeFromIdx(const Vec<T>& vec, const Vec<i32>& keep_indices)
    -> Vec<T>
{
  const std::size_t n = keep_indices.size();
  Vec<T> ret;
  ret.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    if (keep_indices[i]) {
      ret.emplace_back(vec[i]);
    }
  }
  return ret;
}

template <typename T>
static inline auto vec_from_range(i64 low, i64 high) -> Vec<T>
{
  std::size_t elements = 0;
  if (low < 0) {
    elements = high - low + 1;
  } else {
    elements = high + low + 1;
  }
  Vec<T> v(elements);
  i64 i = low;
  for (auto& x : v) {
    x = static_cast<T>(i);
    ++i;
  }
  return v;
}

template <typename T>
static inline auto arange(T start, T stop, T step = 1) -> Vec<T>
{
  // TODO: this is bad if T=unsigned type
  const std::size_t n = std::abs(stop - start);
  Vec<T> values;
  values.reserve(n);

  if (step > 0) {
    for (T value = start; value < stop; value += step) {
      values.emplace_back(value);
    }
  } else {
    for (T value = start; value > stop; value += step) {
      values.emplace_back(value);
    }
  }

  return values;
}

template <typename T>
static inline auto linspace(T start, T stop, std::size_t num = 50,
                            bool endpoint = true) -> Vec<T>
{

  if (num == 0) {
    return Vec<T>(1, 0);
  }

  if (num == 1) {
    Vec<T> values = {start};
    return values;
  }

  if (endpoint) {
    if (num == 2) {
      Vec<T> values = {start, stop};
      return values;
    }

    Vec<T> values(num, 0);
    values[0]       = start;
    values[num - 1] = stop;
    const T step    = (stop - start) / static_cast<T>(num - 1);

    for (std::size_t i = 1; i < num - 1; ++i) {
      values[i] = start + static_cast<T>(i) * step;
    }
    return values;
  }

  if (num == 2) {
    const T step  = (stop - start) / num;
    Vec<T> values = {start, start + step};
    return values;
  }
  Vec<T> values(num, 0);
  values[0] = start;

  const T step = (stop - start) / static_cast<T>(num);

  for (std::size_t i = 1; i < num; ++i) {
    values[i] = start + static_cast<T>(i) * step;
  }

  return values;
}

template <typename ContainerIn, typename ContainerOut = Vec<ContainerIn>>
static inline auto chunks(std::size_t n, const ContainerIn& xs) -> ContainerOut
{
  const auto input_size = xs.size();
  auto ids =
      pft::arange<i64>(0, static_cast<i64>(input_size), static_cast<i64>(n));
  auto N = input_size / n;

  N += input_size % n != 0;
  ids.push_back(static_cast<i64>(input_size));

  ContainerOut ret(N);
  std::transform(std::begin(ret), std::end(ret), std::begin(ret),
                 [i = 0, &xs, &ids](auto) mutable {
                   return take(xs, pft::Slice{ids[i], ids[i++ + 1]});
                 });
  return ret;
}

template <typename T>
static inline auto norm(const Vec<T>& xs)
{
  auto sum2 =
      foldl(T(), xs, [](const T& a, const T& b) -> T { return a + b * b; });
  return std::sqrt(sum2);
}

template <typename T>
static inline auto normalise(const Vec<T>& xs)
{
  const auto d = norm(xs);
  return map([&d](const T& a) -> T { return a / d; }, xs);
}

template <typename T>
static inline auto abs(Vec<T>& vec)
{
  auto abs_lambda = [](auto x) { return std::abs(x); };
  return map(abs_lambda, vec);
}

template <typename T>
static inline auto pad_right(const Vec<T>& in, std::size_t padwidth,
                             T pad_value = 0) -> Vec<T>
{
  const auto out_size = padwidth + in.size();
  Vec<T> out(out_size, pad_value);
  std::copy(std::cbegin(in), std::cend(in), std::begin(out));
  return out;
}

template <typename T>
static inline auto pad_right(const Vec<T>& in, const Slice& slice) -> Vec<T>
{
  const auto out_size = in.size() + slice.size();
  Vec<T> out;
  out.reserve(out_size);
  std::copy(std::cbegin(in), std::cend(in), std::back_inserter(out));
  std::copy(std::cbegin(in) + slice.start, std::cbegin(in) + slice.stop,
            std::back_inserter(out));
  return out;
}

template <typename T>
static inline auto pad_left(const Vec<T>& in, std::size_t padwidth,
                            T pad_value = 0) -> Vec<T>
{
  auto out = pad_right(in, padwidth, pad_value);
  std::rotate(std::begin(out), std::end(out) - padwidth, std::end(out));
  return out;
}

template <typename T>
static inline auto pad(const Vec<T>& in, std::size_t padwidth, T pad_value = 0)
    -> Vec<T>
{
  const auto out_size = 2 * padwidth + in.size();
  Vec<T> out(out_size, pad_value);
  std::copy(std::cbegin(in), std::cend(in), std::begin(out) + padwidth);
  return out;
}

/// Pad a vector from the right with a slice of its elements till that vector
/// has length of max_len
template <typename T>
static inline auto pad_right_until(const Vec<T>& input, const pft::Slice& slice,
                                   std::size_t max_len) -> Vec<T>
{
  if (max_len < input.size()) {
    pft::panic("Padded vector cant have smaller size!");
  }
  Vec<T> out;
  out.reserve(max_len);
  std::copy(std::cbegin(input), std::cend(input), std::back_inserter(out));
  const std::size_t len_missing  = max_len - input.size();
  const std::size_t n_iterations = len_missing / slice.size();
  const std::size_t elements_left_out =
      len_missing - n_iterations * slice.size();
  for (std::size_t i = 0; i < n_iterations; ++i) {
    std::copy(std::cbegin(input) + slice.start, std::cbegin(input) + slice.stop,
              std::back_inserter(out));
  }

  if (elements_left_out != 0) {
    std::copy(std::cbegin(input) + slice.start,
              std::cbegin(input) + slice.start + elements_left_out,
              std::back_inserter(out));
  }

  return out;
}

template <typename T>
static inline auto where(const Vec<i32>& c, const Vec<T>& v1, const Vec<T>& v2)
    -> Vec<T>
{
  const std::size_t n = c.size();
  Vec<T> ret(n);
  for (std::size_t i = 0; i < n; ++i) {
    ret[i] = (c[i] != 0 ? v1[i] : v2[i]);
  }

  return ret;
}

template <typename T>
static inline auto where(const Vec<i32>& c, const Vec<T>& v1, T v2) -> Vec<T>
{
  const std::size_t n = c.size();
  Vec<T> ret(n);
  for (std::size_t i = 0; i < n; ++i) {
    ret[i] = (c[i] != 0 ? v1[i] : v2);
  }

  return ret;
}

template <typename T>
static inline auto where(const Vec<i32>& c, T v1, const Vec<T>& v2) -> Vec<T>
{
  const std::size_t n = c.size();
  Vec<T> ret(n);
  for (std::size_t i = 0; i < n; ++i) {
    ret[i] = (c[i] != 0 ? v1 : v2[i]);
  }

  return ret;
}

template <typename T>
static inline auto where(const Vec<i32>& c, T v1, T v2) -> Vec<T>
{
  const std::size_t n = c.size();
  Vec<T> ret(n);
  for (std::size_t i = 0; i < n; ++i) {
    ret[i] = (c[i] != 0 ? v1 : v2);
  }

  return ret;
}

template <typename T>
static inline auto argsort(const Vec<T>& xs) -> Vec<i32>
{
  auto idx     = arange<i32>(0, static_cast<i32>(xs.size()));
  const auto f = [&xs](std::size_t i1, std::size_t i2) {
    return xs[i1] < xs[i2];
  };
  std::sort(std::begin(idx), std::end(idx), f);

  return idx;
}

// sane mod function
template <typename T>
constexpr static inline auto mod(T a, T b) -> T
{
  return (a % b + b) % b;
}

// return the index
template <typename T, typename Op>
static inline auto findfirst(Op&& fn, const Vec<T>& h) -> Maybe<i64>
{
  auto res = std::find(std::cbegin(h), std::cend(h), fn);
  if (res != h.end()) {
    return {1, std::distance(std::cbegin(h), res)};
  } else {
    return {0, -1};
  }
}

// return the indices
template <typename T, typename Op>
static inline auto findall(Op&& fn, const Vec<T>& h) -> Vec<i64>
{
  const auto n = h.size();
  Vec<i64> ret;
  ret.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (fn(h[i])) {
      ret.emplace_back(i);
    }
  }
  return ret;
}

template <
    class ForwardIt,
    typename Id_t  = typename std::iterator_traits<ForwardIt>::difference_type,
    typename Val_t = typename std::iterator_traits<ForwardIt>::value_type>
constexpr static inline auto max_id_and_val(ForwardIt first, ForwardIt last)
    -> std::pair<Id_t, Val_t>
{
  const auto max_el = std::max_element(first, last);
  const auto i_max  = std::distance(first, max_el);
  return std::make_pair(i_max, *max_el);
}

template <
    class ForwardIt,
    typename Id_t  = typename std::iterator_traits<ForwardIt>::difference_type,
    typename Val_t = typename std::iterator_traits<ForwardIt>::value_type>
constexpr static inline auto min_id_and_val(ForwardIt first, ForwardIt last)
    -> std::pair<Id_t, Val_t>
{
  const auto min_el = std::min_element(first, last);
  const auto i_min  = std::distance(first, min_el);
  return std::make_pair(i_min, *min_el);
}

template <typename T>
static inline bool contains(const Vec<T>& v, T e)
{
  return std::find(std::cbegin(v), std::cend(v), e) != std::cend(v);
}
} // namespace pft

#endif // PFT_H_
