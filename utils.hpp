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
//   0.0.2    NextPowerOf2, deconvln_use_r2c,deconvln_use_c2c, power_spectrum
//   0.0.1    abs_complex, peak_prominences, Peak, peak_widths, local_maxima,
//            select_by_property, unpack_condition_args,
//            select_peaks_by_distance, find_peaks, prod, factorial, vandermonde,
//            SG, vmadd, compute_householder_factor, vnorm, vdiv, QRDecomposition,
//            LUdecomposition, savgol_coeffs, FFTW_R2C_1D, FFTW_C2R_1D,
//            convln_use_r2c, convln_use_c2c
// =============================================================================
#ifndef UTILS_H_
#define UTILS_H_
#include "pft.hpp"
#include <complex>

#ifdef UTILS_USE_FFTW
#include <fftw3.h>
#endif

using pft::Maybe;

namespace utils {
constexpr static inline auto NextPowerOf2(u64 x) -> u64
{
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >> 16);
  x |= (x >> 32);
  return x + 1;
}

template <typename T>
static inline auto abs_complex(const Vec<std::complex<T>>& vec)
{
  auto abs_lambda = [](const auto& x) { return std::abs(x); };
  return pft::map(abs_lambda, vec);
}

static inline auto peak_prominences(const Vec<f64>& xs, const Vec<i32>& peaks,
                                    i32 wlen)
    -> std::tuple<Vec<f64>, Vec<i32>, Vec<i32>>
{
  const std::size_t peaks_size = peaks.size();

  Vec<f64> prominences(peaks_size);
  Vec<i32> left_bases(peaks_size);
  Vec<i32> right_bases(peaks_size);

  f64 left_min, right_min;
  i32 peak, i_min, i_max, i;
  for (std::size_t pi = 0; pi < peaks_size; ++pi) {
    peak  = peaks[pi];
    i_min = 0;
    i_max = static_cast<i32>(xs.size()) - 1;

    if (!(i_min <= peak && peak <= i_max)) {
      fprintf(stderr, "Oopsie\n");
      exit(1);
    }

    if (2 <= wlen) {
      i_min = std::max(peak - wlen / 2, i_min);
      i_max = std::min(peak + wlen / 2, i_max);
    }

    i = left_bases[pi] = peak;
    left_min           = xs[peak];
    while (i_min <= i && xs[i] <= xs[peak]) {
      if (xs[i] < left_min) {
        left_min       = xs[i];
        left_bases[pi] = i;
      }
      --i;
    }

    i = right_bases[pi] = peak;
    right_min           = xs[peak];
    while (i <= i_max && xs[i] <= xs[peak]) {
      if (xs[i] < right_min) {
        right_min       = xs[i];
        right_bases[pi] = i;
      }
      ++i;
    }
    prominences[pi] = xs[peak] - std::max(left_min, right_min);
  }
  return std::make_tuple(prominences, left_bases, right_bases);
}

struct Peak {
  i32 index;
  f64 height;
  f64 width, width_height;
  f64 left_p, right_p;
};

/// Find the width of each peak at a relative height
static inline auto peak_widths(const Vec<f64>& signal,
                               const Vec<i32>& peaks_indices, f64 rel_height)
    -> Vec<Peak>
{
  const std::size_t n_of_peaks = peaks_indices.size();
  const i32 wlen               = -1;
  assert(rel_height > 0);
  const auto [prominences, left_bases, right_bases] =
      peak_prominences(signal, peaks_indices, wlen);
  Vec<f64> widths(n_of_peaks);
  Vec<f64> width_heights(n_of_peaks);
  Vec<f64> left_ips(n_of_peaks);
  Vec<f64> right_ips(n_of_peaks);

  f64 height;
  i32 i_min, i_max, peak_index;

  for (std::size_t i = 0; i < n_of_peaks; ++i) {
    i_min      = left_bases[i];
    i_max      = right_bases[i];
    peak_index = peaks_indices[i];
    if (!(0 <= i_min && i_min <= peak_index && peak_index <= i_max &&
          i_max < static_cast<i32>(signal.size()))) {
      fprintf(stderr, "Oopsie\n");
      exit(1);
    }
    height           = signal[peak_index] - prominences[i] * rel_height;
    width_heights[i] = height;

    // intersection point on left side
    auto p = peak_index;
    while (i_min < p && height < signal[p]) {
      --p;
    }
    auto left_ip = static_cast<f64>(p);

    if (signal[p] < height) {
      left_ip += (height - signal[p]) / (signal[p + 1] - signal[p]);
    }

    // intersection point on right side
    p = peak_index;
    while (p < i_max && height < signal[p]) {
      ++p;
    }

    auto right_ip = (f64)p;
    if (signal[p] < height) {
      right_ip -= (height - signal[p]) / (signal[p - 1] - signal[p]);
    }
    widths[i]    = right_ip - left_ip;
    left_ips[i]  = left_ip;
    right_ips[i] = right_ip;
  }

  Vec<Peak> widths_vec(n_of_peaks);
  for (std::size_t i = 0; i < n_of_peaks; ++i) {
    widths_vec[i] = {peaks_indices[i], signal[peaks_indices[i]],
                     widths[i],        width_heights[i],
                     left_ips[i],      right_ips[i]};
  }
  return widths_vec;
}

template <typename T>
static inline auto local_maxima(Vec<T> x)
    -> std::tuple<Vec<i32>, Vec<i32>, Vec<i32>>
{
  const i32 n = static_cast<i32>(x.size());
  Vec<i32> midpoints(n / 2, 0);
  Vec<i32> left_edges(n / 2, 0);
  Vec<i32> right_edges(n / 2, 0);
  std::size_t m = 0; // index pointer to the end

  i32 i           = 1;
  const i32 i_max = n - 1;

  while (i < i_max) {
    if (x[i - 1] < x[i]) {
      auto i_ahead = i + 1;

      while (i_ahead < i_max && x[i_ahead] == x[i]) {
        ++i_ahead;
      }

      if (x[i_ahead] < x[i]) {
        left_edges[m]  = i;
        right_edges[m] = i_ahead - 1;
        midpoints[m]   = (left_edges[m] + right_edges[m]) / 2;
        ++m;
        i = i_ahead;
      }
    }
    ++i;
  }
  midpoints.resize(m);
  left_edges.resize(m);
  right_edges.resize(m);

  return {midpoints, left_edges, right_edges};
}

template <typename T>
static inline auto select_by_property(Vec<T> p, T pmin, T pmax) -> Vec<i32>
{
  Vec<i32> keep(p.size(), 1);

  // if (pmin.has_value) {
  for (std::size_t i = 0; i < p.size(); ++i) {
    keep[i] &= pmin <= p[i];
  }
  // }

  // if (pmax.has_value) {
  for (std::size_t i = 0; i < p.size(); ++i) {
    keep[i] &= p[i] <= pmax;
  }
  // }
  return keep;
}

template <typename T>
static inline auto unpack_condition_args(const std::pair<T, T>& interval,
                                         const Vec<f64>& xs,
                                         const Vec<i32>& peaks)
    -> std::pair<f64, f64>
{

  (void)xs;
  (void)peaks;
  // TODO: implement unpacking for when T is a container
  auto [imin, imax] = interval;

  return std::make_pair(imin, imax);
}

static inline auto select_peaks_by_distance(const Vec<i32>& peaks,
                                            const Vec<f64>& priority,
                                            f64 distance) -> Vec<i32>
{

  const i32 peaks_size = static_cast<i32>(peaks.size());

  distance = ceil(distance);
  Vec<i32> keep(peaks_size, true);

  auto priority_to_position = pft::argsort(priority);
  i32 j                     = 0;
  for (i32 i = peaks_size - 1; i > -1; --i) {
    j = priority_to_position[i];
    if (!keep[j]) {
      continue;
    }

    auto k = j - 1;
    while (0 <= k && (peaks[j] - peaks[k]) < distance) {
      keep[k] = false;
      --k;
    }
    k = j + 1;
    while (k < peaks_size && (peaks[k] - peaks[j]) < distance) {
      keep[k] = false;
      ++k;
    }
  }

  return keep;
}

/// Return the indices of the peaks
static inline auto find_peaks(const Vec<f64>& xs,
                              Maybe<std::pair<f64, f64>> height,
                              f64 distance = 0.0) -> Vec<i32>
{
  auto [peaks_indices, l_edges, r_edges] = local_maxima(xs);

  if (height.has_value) {
    auto peak_heights = pft::take(xs, peaks_indices);
    auto [hmin, hmax] = unpack_condition_args(height.unwrap, xs, peaks_indices);
    auto keep         = select_by_property(peak_heights, hmin, hmax);

    peaks_indices = pft::takeFromIdx(peaks_indices, keep);
  }
  // distance between peaks
  if (distance != 0.0) {
    auto keep = select_peaks_by_distance(
        peaks_indices, pft::take(xs, peaks_indices), distance);
    // TODO: properly test that
    peaks_indices = pft::takeFromIdx(peaks_indices, keep);
  }

  return peaks_indices;
}

template <typename T>
static inline auto prod(const Vec<T>& x) -> T
{
  auto ret =
      pft::foldl(T(1), x, [](const T& a, const T& b) -> T { return a * b; });
  return ret;
}

static inline auto factorial(i64 n) -> f64
{
  static f64 table[171];
  static bool init = true;
  if (init) {
    init     = false;
    table[0] = 1;
    for (std::size_t i = 1; i < 171; ++i) {
      table[i] = static_cast<f64>(i) * table[i - 1];
    }
  }
  return table[n];
}

template <typename T>
static inline auto vandermonde(i64 halfWindow, i64 polyDeg) -> pft::Matrix<T>
{
  assert(halfWindow >= 0);
  assert(polyDeg >= 0);
  // arange return an exclusive range
  auto x = pft::arange<T>(-halfWindow, halfWindow + 1);

  std::size_t n = polyDeg + 1;
  auto m        = x.size();

  pft::Matrix<T> V(m, n);

  for (std::size_t i = 0; i < m; ++i) {
    V(i, 0) = T(1);
  }

  for (std::size_t j = 1; j < n; ++j) {
    for (std::size_t i = 0; i < m; ++i) {
      V(i, j) = x[i] * V(i, j - 1);
    }
  }

  return V;
}

template <typename T>
static inline auto SG(i64 halfWindow, i64 polyDeg) -> pft::Matrix<T>
{
  assert(2 * halfWindow > polyDeg);

  auto V = vandermonde<T>(halfWindow, polyDeg);

  // QR decomposition
  auto [Q, R] = qr(V);

  auto SG = R / Q.transpose();

  auto n = SG.rows;
  auto m = SG.cols;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      SG(i, j) *= factorial(static_cast<i64>(i));
    }
  }

  return SG.transpose();
}

/// Calculate Vi = Ai + c*Bi
template <typename T>
static inline auto vmadd(const Vec<T>& a, const Vec<T>& b, T s) -> Vec<T>
{
  const auto n = a.size();
  Vec<T> ret(n);
  for (std::size_t i = 0; i < n; ++i) {
    ret[i] = a[i] + s * b[i];
  }

  return ret;
}

template <typename T>
static inline auto compute_householder_factor(const Vec<T>& v) -> pft::Matrix<T>
{
  const auto n = v.size();
  pft::Matrix<T> ret(n, n);

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      ret(i, j) = -2 * v[i] * v[j];
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    ret(i, i) += 1;
  }
  return ret;
}

/// Return the norm of a vector
template <typename T>
static inline auto vnorm(const Vec<T>& xs) -> f64
{
  return std::sqrt(std::transform_reduce(std::cbegin(xs), std::cend(xs), T(),
                                         std::plus{},
                                         [](const auto& x) { return x * x; }));
}

/// Divide a vector with a value
template <typename T>
static inline auto vdiv(const Vec<T>& x, T d) -> Vec<T>
{
  auto ret = pft::map([&d](const T& a) -> T { return a / d; }, x);
  return ret;
}

template <typename T>
static inline auto QRDecomposition(const pft::Matrix<T>& mat)
    -> std::pair<pft::Matrix<T>, pft::Matrix<T>>
{
  pft::Matrix<T> Q;
  pft::Matrix<T> R;

  const std::size_t m = mat.rows;
  const std::size_t n = mat.cols;

  // vector if factors Q1. Q2.... Qm
  Vec<pft::Matrix<T>> qv(m);

  pft::Matrix<T> z = mat;
  pft::Matrix<T> z1;
  for (std::size_t k = 0; k < n && k < m - 1; ++k) {
    Vec<T> e(m, 0), x(m, 0);
    f64 a;

    z1 = z.GetMinor(k);
    z  = z1;

    x = z.getColumn(k);
    a = pft::norm(x);

    if (mat(k, k) > 0) {
      a = -a;
    }

    for (std::size_t i = 0; i < m; ++i) {
      e[i] = (i == k) ? 1 : 0;
    }

    e     = vmadd(x, e, a);
    e     = vdiv(e, pft::norm(e));
    qv[k] = compute_householder_factor(e);

    z1 = qv[k].mult(z);
    z  = z1;
  }
  Q = qv[0];

  for (std::size_t i = 1; i < n && i < m - 1; ++i) {
    z1 = qv[i].mult(Q);
    Q  = z1;
  }
  R = Q.mult(mat);

  return {Q.transpose(), R};
}

struct LUdecomposition {
  using Mat_t = pft::Matrix<f64>;
  size_t n;
  Mat_t lu;
  Vec<i32> indices;
  f64 d;

  LUdecomposition(const Mat_t& a) : n(a.rows), lu(a), indices(n)
  {
    constexpr f64 tiny = 1.0e-40;
    size_t i, imax, j, k;
    f64 big, temp;
    Vec<f64> vv(n, 0.0);
    d = 1.0;
    for (i = 0; i < n; ++i) {
      big = 0.0;
      for (j = 0; j < n; ++j) {
        temp = std::abs(lu(i, j));
        if (temp > big) {
          big = temp;
        }
      }
      vv[i] = 1.0 / big;
    }

    for (k = 0; k < n; ++k) {
      big = 0.0;
      for (i = k; i < n; ++i) {
        temp = vv[i] * std::abs(lu(i, k));
        if (temp > big) {
          big  = temp;
          imax = i;
        }
      }
      if (k != imax) {
        for (j = 0; j < n; ++j) {
          temp        = lu(imax, j);
          lu(imax, j) = lu(k, j);
          lu(k, j)    = temp;
        }
        d        = -d;
        vv[imax] = vv[k];
      }
      indices[k] = static_cast<i32>(imax);
      if (lu(k, k) == 0.0) {
        lu(k, k) = tiny;
      }
      for (i = k + 1; i < n; ++i) {
        lu(i, k) /= lu(k, k);
        temp = lu(i, k);
        for (j = k + 1; j < n; ++j) {
          lu(i, j) -= temp * lu(k, j);
        }
      }
    }
  }
  ~LUdecomposition() = default;

  auto solve(const Vec<f64>& b) -> Vec<f64>
  {
    Vec<f64> x(b);

    std::size_t ii = 0, ip, j;
    f64 sum;
    if (b.size() != n || x.size() != n) {
      fprintf(stderr, "Bad sizes\n");
      exit(1);
    }

    for (std::size_t i = 0; i < n; ++i) {
      ip    = indices[i];
      sum   = x[ip];
      x[ip] = x[i];
      if (ii != 0) {
        for (j = ii - 1; j < i; ++j) {
          sum -= lu(i, j) * x[j];
        }
      } else if (sum != 0.0) {
        ii = i + 1;
      }
      x[i] = sum;
    }

    for (i32 i = static_cast<i32>(n) - 1; i >= 0; --i) {
      sum = x[i];
      for (j = i + 1; j < n; ++j) {
        sum -= lu(i, j) * x[j];
      }
      x[i] = sum / lu(i, i);
    }

    return x;
  }

  void solve(Mat_t& b, Mat_t& x)
  {
    std::size_t i, m = b.cols;
    Vec<f64> xx(n);
    for (std::size_t j = 0; j < m; ++j) {
      for (i = 0; i < n; ++i) {
        xx[i] = b(i, j);
      }
      xx = solve(xx);
      for (i = 0; i < n; ++i) {
        x(i, j) = xx[i];
      }
    }
  }

  void inverse(Mat_t& ainv)
  {
    ainv = Mat_t(n, n);
    ainv.diagonal();
    solve(ainv, ainv);
  }
  // f64 det();
};

/// Savitzky-Golay Coeffs for 1-D filter
/// np: n points, must be an odd number
/// nl, nr: n points to left and right
/// ld: order of the derivative, default = 0 (no derivation)
/// m: order of the polynomial
static inline auto savgol_coeffs(const i32 np, const i32 nl, const i32 nr,
                                 const i32 ld, const i32 m) -> Vec<f64>
{

  if (np < nl + nr + 1 || nr < 0 || nr < 0 || ld > m || nl + nr < m) {
    fprintf(stderr, "Bad arguments\n");
    exit(1);
  }

  i32 k, mm, imj;
  // i32 kk;
  f64 sum, fac;
  Vec<i32> indx(m + 1);
  pft::Matrix<f64> a(m + 1, m + 1);

  for (i32 ipj = 0; ipj <= (m << 1); ++ipj) {
    sum = (ipj ? 0.0 : 1.0);
    for (k = 1; k <= nr; ++k) {
      sum += std::pow((f64)k, (f64)ipj);
    }
    for (k = 1; k <= nr; ++k) {
      sum += std::pow((f64)(-k), (f64)ipj);
    }

    mm = std::min(ipj, 2 * m - ipj);
    for (imj = -mm; imj <= mm; imj += 2) {
      a((ipj + imj) / 2, (ipj - imj) / 2) = sum;
    }
  }
  LUdecomposition alud(a);
  Vec<f64> b(m + 1, 0.0);
  b[ld] = 1.0;
  b     = alud.solve(b);

  std::size_t i = 0;
  Vec<f64> c(np, 0.0);
  for (k = -nl; k <= nr; ++k) {
    sum = b[0];
    fac = 1.0;
    for (mm = 1; mm <= m; ++mm) {
      fac *= k;
      sum += b[mm] * fac;
    }
    // kk    = (np - k) % np; // beware of c++ modulo
    // c[kk] = sum;
    c[i++] = sum;
  }

  return c;
}

#ifdef UTILS_USE_FFTW
struct FFTW_R2C_1D {
  std::size_t input_size{0};
  std::size_t output_size{0};
  f64* input_buffer;
  fftw_complex* output_buffer = nullptr;

  FFTW_R2C_1D(std::size_t n) : input_size(n), output_size(n / 2 + 1)
  {
    input_buffer  = fftw_alloc_real(input_size);
    output_buffer = fftw_alloc_complex(output_size);
    p = fftw_plan_dft_r2c_1d(n, input_buffer, output_buffer, FFTW_ESTIMATE);
  }

  ~FFTW_R2C_1D()
  {
    fftw_destroy_plan(p);
    fftw_free(input_buffer);
    fftw_free(output_buffer);
  }

  void set_input(const Vec<f64>& vec)
  {
    assert(vec.size() == input_size);
    memcpy(input_buffer, vec.data(), sizeof(f64) * vec.size());
  }

  void set_input_zeropadded(const f64* buffer, std::size_t size)
  {
    assert(size <= input_size);
    memcpy(input_buffer, buffer, sizeof(f64) * size);
    memset(&input_buffer[size], 0, sizeof(f64) * (input_size - size));
  }

  void set_input_zeropadded(const Vec<f64>& vec)
  {
    set_input_zeropadded(&vec[0], vec.size());
  }

  void set_input_pad_with_itself(const Vec<f64>& vec, pft::Slice slice)
  {
    const auto buffer = pft::pad_right(vec, slice);
    assert(buffer.size() == input_size);
    memcpy(input_buffer, buffer.data(), sizeof(f64) * buffer.size());
  }

  void execute() { fftw_execute(p); }

  auto get_output_as_array() -> fftw_complex* { return output_buffer; }

  auto get_output_as_vec() const -> Vec<std::complex<f64>>
  {
    Vec<std::complex<f64>> ret(output_size);

    for (std::size_t i = 0; i < output_size; ++i) {
      ret[i] = {output_buffer[i][0], output_buffer[i][1]};
    }
    return ret;
  }

private:
  fftw_plan p;
};

struct FFTW_C2R_1D {
  std::size_t input_size{0};
  std::size_t output_size{0};
  fftw_complex* input_buffer;
  f64* output_buffer;
  i32 total_size;

  FFTW_C2R_1D(std::size_t n_real_samples)
      : input_size(n_real_samples / 2 + 1), output_size(n_real_samples)
  {
    input_buffer  = fftw_alloc_complex(input_size);
    output_buffer = fftw_alloc_real(output_size);
    p = fftw_plan_dft_c2r_1d(n_real_samples, input_buffer, output_buffer,
                             FFTW_ESTIMATE);
  }

  ~FFTW_C2R_1D()
  {
    fftw_destroy_plan(p);
    fftw_free(input_buffer);
    fftw_free(output_buffer);
  }

  void set_input(const fftw_complex* buffer)
  {
    memcpy(input_buffer, buffer, sizeof(fftw_complex) * input_size);
  }

  void set_points(const fftw_complex* buffer)
  {
    const i32 sizein = i32(f64(total_size) * (input_size / 2 + 1) / input_size);
    for (i32 i = 0; i < sizein; ++i) {
      input_buffer[i][0] = buffer[i][0];
      input_buffer[i][1] = buffer[i][1];
    }
  }

  void set_input_zeropadded(const fftw_complex* buffer, std::size_t size)
  {
    assert(size <= input_size);
    memcpy(input_buffer, buffer, sizeof(fftw_complex) * size);
    memset(&input_buffer[size], 0, sizeof(fftw_complex) * (input_size - size));
  }

  void execute() { fftw_execute(p); }

  auto get_output_as_vec() const -> Vec<f64>
  {
    Vec<f64> ret(output_buffer, output_buffer + output_size);
    return ret;
  }

  auto get_normalised_output_as_vec() const -> Vec<f64>
  {
    const auto ret = pft::map(
        [&](const auto& x) { return x / static_cast<f64>(output_size); },
        Vec<f64>(output_buffer, output_buffer + output_size));
    return ret;
  }

private:
  fftw_plan p;
};

static inline auto convln_use_r2c(const Vec<f64>& input, const Vec<f64>& kernel)
    -> Vec<f64>
{
  const auto n                = input.size();
  const auto m                = kernel.size();
  const std::size_t fft_shape = n + m - 1;

  const auto padded_input  = pft::pad_right(input, fft_shape - n);
  const auto padded_kernel = pft::pad_right(kernel, fft_shape - m);

  FFTW_R2C_1D input_fft(fft_shape);
  input_fft.set_input(padded_input);
  input_fft.execute();

  FFTW_R2C_1D kernel_fft(fft_shape);
  kernel_fft.set_input(padded_kernel);
  kernel_fft.execute();

  const auto input_fft_output  = input_fft.get_output_as_array();
  const auto kernel_fft_output = kernel_fft.get_output_as_array();

  // mutliply the ffts
  const auto temp_size        = input_fft.output_size;
  fftw_complex* fixed_product = fftw_alloc_complex(temp_size);

  const f64 mux = 1.0 / static_cast<f64>(fft_shape);
  for (std::size_t i = 0; i < temp_size; ++i) {
    fixed_product[i][0] = (input_fft_output[i][0] * kernel_fft_output[i][0] -
                           input_fft_output[i][1] * kernel_fft_output[i][1]) *
                          mux;
    fixed_product[i][1] = (input_fft_output[i][0] * kernel_fft_output[i][1] +
                           input_fft_output[i][1] * kernel_fft_output[i][0]) *
                          mux;
  }

  // iff of the product
  FFTW_C2R_1D product_fft(fft_shape);
  product_fft.set_input(fixed_product);
  product_fft.execute();

  const auto output = product_fft.get_output_as_vec();

  fftw_free(fixed_product);
  return output;
}

static inline auto convln_use_c2c(const Vec<f64>& input, const Vec<f64>& kernel)
    -> Vec<f64>
{

  const auto n                       = input.size();
  const auto m                       = kernel.size();
  const std::size_t convolution_size = n + m - 1;

  const auto padded_input  = pft::pad_right(input, convolution_size - n);
  const auto padded_kernel = pft::pad_right(kernel, convolution_size - m);

  fftw_complex* input_arr_in    = fftw_alloc_complex(convolution_size);
  fftw_complex* kernel_arr_in   = fftw_alloc_complex(convolution_size);
  fftw_complex* input_arr_out   = fftw_alloc_complex(convolution_size);
  fftw_complex* kernel_arr_out  = fftw_alloc_complex(convolution_size);
  fftw_complex* fixed_product   = fftw_alloc_complex(convolution_size);
  fftw_complex* convolution_out = fftw_alloc_complex(convolution_size);

  fftw_plan p1 = fftw_plan_dft_1d(convolution_size, input_arr_in, input_arr_out,
                                  FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan p2 = fftw_plan_dft_1d(convolution_size, kernel_arr_in,
                                  kernel_arr_out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan rev =
      fftw_plan_dft_1d(convolution_size, fixed_product, convolution_out,
                       FFTW_BACKWARD, FFTW_ESTIMATE);

  for (std::size_t i = 0; i < convolution_size; ++i) {
    input_arr_in[i][0]  = padded_input[i];
    input_arr_in[i][1]  = 0.0;
    kernel_arr_in[i][0] = padded_kernel[i];
    kernel_arr_in[i][1] = 0.0;
  }

  fftw_execute(p1);
  fftw_execute(p2);
  const f64 mux = 1.0 / static_cast<f64>(convolution_size);
  for (std::size_t i = 0; i < convolution_size; ++i) {
    fixed_product[i][0] = (input_arr_out[i][0] * kernel_arr_out[i][0] -
                           input_arr_out[i][1] * kernel_arr_out[i][1]) *
                          mux;
    fixed_product[i][1] = (input_arr_out[i][0] * kernel_arr_out[i][1] +
                           input_arr_out[i][1] * kernel_arr_out[i][0]) *
                          mux;
  }

  fftw_execute(rev);

  Vec<f64> ret(convolution_size);
  for (std::size_t i = 0; i < convolution_size; ++i) {
    ret[i] = convolution_out[i][0];
  }

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_destroy_plan(rev);
  fftw_free(input_arr_in);
  fftw_free(input_arr_out);
  fftw_free(kernel_arr_in);
  fftw_free(kernel_arr_out);
  fftw_free(fixed_product);
  fftw_free(convolution_out);
  return ret;
}

static inline auto deconvln_use_r2c(const Vec<f64>& a, const Vec<f64>& b)
    -> Vec<f64>
{
  const auto n           = a.size();
  const auto m           = b.size();
  const auto deconv_size = n - m + 1;
  const auto pad_size    = NextPowerOf2(std::max(n, m));

  const auto padded_a = pft::pad_right(a, pad_size - n);
  const auto padded_b = pft::pad_right(b, pad_size - m);

  utils::FFTW_R2C_1D a_fft(pad_size);
  a_fft.set_input(padded_a);
  a_fft.execute();

  utils::FFTW_R2C_1D b_fft(pad_size);
  b_fft.set_input(padded_b);
  b_fft.execute();

  const auto a_fft_out = a_fft.get_output_as_array();
  const auto b_fft_out = b_fft.get_output_as_array();

  const auto fft_size   = a_fft.output_size;
  fftw_complex* div_arr = fftw_alloc_complex(pad_size);

  const f64 mux = 1.0 / static_cast<f64>(pad_size);
  for (std::size_t i = 0; i < fft_size; ++i) {
    const auto mag2 =
        b_fft_out[i][0] * b_fft_out[i][0] + b_fft_out[i][1] * b_fft_out[i][1];
    div_arr[i][0] = (a_fft_out[i][0] * b_fft_out[i][0] +
                     a_fft_out[i][1] * b_fft_out[i][1]) /
                    mag2 * mux;
    div_arr[i][1] = (-a_fft_out[i][0] * b_fft_out[i][1] +
                     a_fft_out[i][1] * b_fft_out[i][0]) /
                    mag2 * mux;
  }

  // iff of the product
  utils::FFTW_C2R_1D div_fft(pad_size);
  div_fft.set_input(div_arr);
  div_fft.execute();

  const Vec<f64> res(div_fft.output_buffer,
                     div_fft.output_buffer + deconv_size);

  fftw_free(div_arr);
  return res;
}

static inline auto deconvln_use_c2c(const Vec<f64>& a, const Vec<f64>& b)
    -> Vec<f64>
{
  const auto n           = a.size();
  const auto m           = b.size();
  const auto deconv_size = n - m + 1;
  const auto pad_size    = NextPowerOf2(std::max(n, m));

  const auto padded_a = pft::pad_right(a, pad_size - n);
  const auto padded_b = pft::pad_right(b, pad_size - m);

  fftw_complex* a_cmp_arr = fftw_alloc_complex(pad_size);
  fftw_complex* b_cmp_arr = fftw_alloc_complex(pad_size);
  fftw_complex* a_fft_arr = fftw_alloc_complex(pad_size);
  fftw_complex* b_fft_arr = fftw_alloc_complex(pad_size);
  fftw_complex* div_arr   = fftw_alloc_complex(pad_size);
  f64* deconv_res         = fftw_alloc_real(pad_size);

  fftw_plan p1 = fftw_plan_dft_1d(pad_size, a_cmp_arr, a_fft_arr, FFTW_FORWARD,
                                  FFTW_ESTIMATE);
  fftw_plan p2 = fftw_plan_dft_1d(pad_size, b_cmp_arr, b_fft_arr, FFTW_FORWARD,
                                  FFTW_ESTIMATE);
  fftw_plan rev =
      fftw_plan_dft_c2r_1d(pad_size, div_arr, deconv_res, FFTW_ESTIMATE);

  for (std::size_t i = 0; i < pad_size; ++i) {
    a_cmp_arr[i][0] = padded_a[i];
    a_cmp_arr[i][1] = 0.0;
    b_cmp_arr[i][0] = padded_b[i];
    b_cmp_arr[i][1] = 0.0;
  }

  fftw_execute(p1);
  fftw_execute(p2);
  const f64 mux       = 1.0 / static_cast<f64>(pad_size);
  const auto fft_size = (pad_size >> 1);
  for (std::size_t i = 0; i < fft_size; ++i) {
    const auto mag2 =
        b_fft_arr[i][0] * b_fft_arr[i][0] + b_fft_arr[i][1] * b_fft_arr[i][1];
    div_arr[i][0] = (a_fft_arr[i][0] * b_fft_arr[i][0] +
                     a_fft_arr[i][1] * b_fft_arr[i][1]) *
                    mux / mag2;
    div_arr[i][1] = (-a_fft_arr[i][0] * b_fft_arr[i][1] +
                     a_fft_arr[i][1] * b_fft_arr[i][0]) *
                    mux / mag2;
  }

  fftw_execute(rev);

  Vec<f64> res(deconv_size);
  for (std::size_t i = 0; i < deconv_size; ++i) {
    res[i] = deconv_res[i] + deconv_res[pad_size - i - 1];
  }

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_destroy_plan(rev);
  fftw_free(a_cmp_arr);
  fftw_free(a_fft_arr);
  fftw_free(b_cmp_arr);
  fftw_free(b_fft_arr);
  fftw_free(div_arr);
  fftw_free(deconv_res);
  return res;
}

// TODO: test deconvolution corectly
static inline auto deconvln_use_r2c_wraparound_kernel(const Vec<f64>& input,
                                                      const Vec<f64>& kernel)
    -> Vec<f64>
{
  const auto n                = input.size();
  const auto m                = kernel.size();
  const std::size_t fft_shape = input.size() + kernel.size() - 1;

  const auto padded_input  = pft::pad_right(input, fft_shape - n);
  const auto padded_kernel = pft::pad_right(kernel, fft_shape - m);

  FFTW_R2C_1D input_fft(fft_shape);
  input_fft.set_input(padded_input);
  input_fft.execute();

  FFTW_R2C_1D kernel_fft(fft_shape);
  kernel_fft.set_input(padded_kernel);
  kernel_fft.execute();

  const auto input_fft_output  = input_fft.get_output_as_array();
  const auto kernel_fft_output = kernel_fft.get_output_as_array();

  fftw_complex* fixed_product = fftw_alloc_complex(fft_shape);

  for (std::size_t i = 0; i <= fft_shape / 2; ++i) {
    const auto mag2 = input_fft_output[i][0] * input_fft_output[i][0] +
                      input_fft_output[i][1] * input_fft_output[i][1];
    fixed_product[i][0] = (input_fft_output[i][0] * kernel_fft_output[i][0] -
                           input_fft_output[i][1] * kernel_fft_output[i][1]) /
                          mag2;
    fixed_product[i][1] = (input_fft_output[i][0] * kernel_fft_output[i][1] +
                           input_fft_output[i][1] * kernel_fft_output[i][0]) /
                          mag2;
  }

  for (std::size_t i = fft_shape / 2 + 1; i < fft_shape; ++i) {
    const auto mag2 =
        input_fft_output[fft_shape - i][0] *
            input_fft_output[fft_shape - i][0] +
        input_fft_output[fft_shape - i][1] * input_fft_output[fft_shape - i][1];
    fixed_product[i][0] = (input_fft_output[fft_shape - i][0] *
                               kernel_fft_output[fft_shape - i][0] -
                           input_fft_output[fft_shape - i][1] *
                               kernel_fft_output[fft_shape - i][1]) /
                          mag2;
    fixed_product[i][1] = (-(input_fft_output[fft_shape - i][0] *
                                 kernel_fft_output[fft_shape - i][1] +
                             input_fft_output[fft_shape - i][1] *
                                 kernel_fft_output[fft_shape - i][0])) /
                          mag2;
  }

  // iff of the product
  FFTW_C2R_1D product_fft(fft_shape);
  product_fft.set_points(fixed_product);
  product_fft.execute();

  const auto output = product_fft.get_normalised_output_as_vec();

  fftw_free(fixed_product);
  return output;
}

static inline auto power_spectrum(const Vec<f64>& timesteps,
                                  const Vec<f64>& signal)
    -> std::pair<Vec<f64>, Vec<f64>>
{
  const auto N        = signal.size();
  const auto timestep = timesteps[1] - timesteps[0];

  utils::FFTW_R2C_1D fft_signal(N);
  fft_signal.set_input(signal);
  fft_signal.execute();

  const auto fft_output = fft_signal.get_output_as_vec();

  const auto wf_power_spectrum =
      pft::map([&](const auto& x) { return x * x / (f64)N; },
               utils::abs_complex(fft_output));
  const auto freqs = pft::arange<f64>(0.0, timestep * 0.5, timestep / (f64)N);

  return std::make_pair(freqs, wf_power_spectrum);
}

#endif
} // namespace utils
#endif // UTILS_H_
