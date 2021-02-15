#ifndef UTILS_H_
#define UTILS_H_
#include "pft.hpp"
#include <complex>

using pft::Maybe;

namespace utils {
// so that we dont have to include complex in pft
template <typename T>
static inline auto abs_complex(const std::vector<std::complex<T>>& vec) {
  auto abs_lambda = [](const auto& x) { return std::abs(x); };
  return pft::map(abs_lambda, vec);
}

auto peak_prominences(const std::vector<f64>& xs, const std::vector<i32>& peaks,
                      i32 wlen)
    -> std::tuple<std::vector<f64>, std::vector<i32>, std::vector<i32>> {
  const std::size_t peaks_size = peaks.size();

  std::vector<f64> prominences(peaks_size);
  std::vector<i32> left_bases(peaks_size);
  std::vector<i32> right_bases(peaks_size);

  f64 left_min, right_min;
  i32 peak, i_min, i_max, i;
  for (std::size_t pi = 0; pi < peaks_size; ++pi) {
    peak  = peaks[pi];
    i_min = 0;
    i_max = xs.size() - 1;

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
auto peak_widths(const std::vector<f64>& signal,
                 const std::vector<i32>& peaks_indices, f64 rel_height)
    -> std::vector<Peak> {
  const std::size_t n_of_peaks = peaks_indices.size();
  const i32 wlen               = -1;
  assert(rel_height > 0);
  const auto [prominences, left_bases, right_bases] =
      peak_prominences(signal, peaks_indices, wlen);
  std::vector<f64> widths(n_of_peaks);
  std::vector<f64> width_heights(n_of_peaks);
  std::vector<f64> left_ips(n_of_peaks);
  std::vector<f64> right_ips(n_of_peaks);

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
    auto left_ip = (f64)p;

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

  std::vector<Peak> widths_vec(n_of_peaks);
  for (std::size_t i = 0; i < n_of_peaks; ++i) {
    widths_vec[i] = {peaks_indices[i], signal[peaks_indices[i]],
                     widths[i],        width_heights[i],
                     left_ips[i],      right_ips[i]};
  }
  return widths_vec;
}

template <typename T>
auto local_maxima(std::vector<T> x)
    -> std::tuple<std::vector<i32>, std::vector<i32>, std::vector<i32>> {
  const std::size_t n = x.size();
  std::vector<i32> midpoints(n / 2, 0);
  std::vector<i32> left_edges(n / 2, 0);
  std::vector<i32> right_edges(n / 2, 0);
  std::size_t m = 0; // index pointer to the end

  std::size_t i = 1;
  auto i_max    = n - 1;

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
auto select_by_property(std::vector<T> p, T pmin, T pmax) -> std::vector<i32> {
  std::vector<i32> keep(p.size(), 1);

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
auto unpack_condition_args(const std::pair<T, T>& interval,
                           const std::vector<f64>& xs,
                           const std::vector<i32>& peaks)
    -> std::pair<f64, f64> {

  (void)xs;
  (void)peaks;
  // TODO: implement unpacking for when T is a container
  auto [imin, imax] = interval;

  return std::make_pair(imin, imax);
}

auto select_peaks_by_distance(const std::vector<i32>& peaks,
                              const std::vector<f64>& priority, f64 distance)
    -> std::vector<i32> {

  i32 peaks_size = peaks.size();

  distance = ceil(distance);
  std::vector<i32> keep(peaks_size, true);

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
auto find_peaks(const std::vector<f64>& xs, Maybe<std::pair<f64, f64>> height,
                f64 distance = 0.0) -> std::vector<i32> {
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
auto prod(const std::vector<T>& x) -> T {
  auto ret =
      pft::foldl(T(1), x, [](const T& a, const T& b) -> T { return a * b; });
  return ret;
}

auto factorial(const i64& n) -> f64 {
  static f64 table[171];
  static bool init = true;
  if (init) {
    init     = false;
    table[0] = 1;
    for (std::size_t i = 1; i < 171; ++i) {
      table[i] = (f64)i * table[i - 1];
    }
  }
  return table[n];
}

template <typename T>
auto vandermonde(i64 halfWindow, i64 polyDeg) -> pft::Matrix<T> {
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
auto SG(i64 halfWindow, i64 polyDeg) -> pft::Matrix<T> {
  assert(2 * halfWindow > polyDeg);

  auto V = vandermonde<T>(halfWindow, polyDeg);

  // QR decomposition
  auto [Q, R] = qr(V);

  auto SG = R / Q.transpose();

  auto n = SG.rows;
  auto m = SG.cols;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      SG(i, j) *= factorial(i);
    }
  }

  return SG.transpose();
}

/// Calculate Vi = Ai + c*Bi
template <typename T>
auto vmadd(const std::vector<T>& a, const std::vector<T>& b, T s)
    -> std::vector<T> {
  const auto n = a.size();
  std::vector<T> ret(n);
  for (std::size_t i = 0; i < n; ++i) {
    ret[i] = a[i] + s * b[i];
  }

  return ret;
}

template <typename T>
auto compute_householder_factor(const std::vector<T>& v) -> pft::Matrix<T> {
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
auto vnorm(const std::vector<T>& xs) -> f64 {
  return std::sqrt(std::transform_reduce(std::cbegin(xs), std::cend(xs), T(),
                                         std::plus{},
                                         [](const auto& x) { return x * x; }));
}

/// Divide a vector with a value
template <typename T>
auto vdiv(const std::vector<T>& x, T d) -> std::vector<T> {
  auto ret = pft::map([&d](const T& a) -> T { return a / d; }, x);
  return ret;
}

template <typename T>
auto QRDecomposition(const pft::Matrix<T>& mat)
    -> std::pair<pft::Matrix<T>, pft::Matrix<T>> {
  pft::Matrix<T> Q;
  pft::Matrix<T> R;

  const std::size_t m = mat.rows;
  const std::size_t n = mat.cols;

  // vector if factors Q1. Q2.... Qm
  std::vector<pft::Matrix<T>> qv(m);

  pft::Matrix<T> z = mat;
  pft::Matrix<T> z1;
  for (std::size_t k = 0; k < n && k < m - 1; ++k) {
    std::vector<T> e(m, 0), x(m, 0);
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
  std::size_t n;
  Mat_t lu;
  std::vector<i32> indices;
  f64 d;
  LUdecomposition(const Mat_t& a);
  ~LUdecomposition() = default;
  auto solve(const std::vector<f64>& b) -> std::vector<f64>;
  void solve(Mat_t& b, Mat_t& x);
  void inverse(Mat_t& ainv);
  // f64 det();
};

LUdecomposition::LUdecomposition(const Mat_t& a)
    : n(a.rows), lu(a), indices(n) {

  constexpr f64 tiny = 1.0e-40;
  std::size_t i, imax, j, k;
  f64 big, temp;
  std::vector<f64> vv(n, 0.0);
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
    indices[k] = imax;
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

auto LUdecomposition::solve(const std::vector<f64>& b) -> std::vector<f64> {
  std::vector<f64> x(b);

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

  for (i32 i = n - 1; i >= 0; --i) {
    sum = x[i];
    for (j = i + 1; j < n; ++j) {
      sum -= lu(i, j) * x[j];
    }
    x[i] = sum / lu(i, i);
  }

  return x;
}

void LUdecomposition::solve(Mat_t& b, Mat_t& x) {
  std::size_t i, m = b.cols;
  std::vector<f64> xx(n);
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

void LUdecomposition::inverse(Mat_t& ainv) {
  ainv = Mat_t(n, n);
  ainv.diagonal();
  solve(ainv, ainv);
}

/// Savitzky-Golay Coeffs for 1-D filter
/// np: n points, must be an odd number
/// nl, nr: n points to left and right
/// ld: order of the derivative, default = 0 (no derivation)
/// m: order of the polynomial
auto savgol_coeffs(const i32 np, const i32 nl, const i32 nr, const i32 ld,
                   const i32 m) -> std::vector<f64> {

  if (np < nl + nr + 1 || nr < 0 || nr < 0 || ld > m || nl + nr < m) {
    fprintf(stderr, "Bad arguments\n");
    exit(1);
  }

  i32 k, mm, imj;
  // i32 kk;
  f64 sum, fac;
  std::vector<i32> indx(m + 1);
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
  std::vector<f64> b(m + 1, 0.0);
  b[ld] = 1.0;
  b     = alud.solve(b);

  std::size_t i = 0;
  std::vector<f64> c(np, 0.0);
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
#include <complex>
#include <fftw3.h>

struct FFTW_R2C_1D {
  std::size_t input_size{0};
  f64* const input_buffer;
  std::size_t output_size{0};
  fftw_complex* output_buffer = nullptr;

  FFTW_R2C_1D(std::size_t fft_size)
      : input_size(fft_size), input_buffer(fftw_alloc_real(fft_size)),
        output_size(input_size / 2 + 1),
        output_buffer(fftw_alloc_complex(output_size)) {
    p = fftw_plan_dft_r2c_1d(input_size, input_buffer, output_buffer,
                             FFTW_ESTIMATE);
  }

  ~FFTW_R2C_1D() {
    fftw_destroy_plan(p);
    fftw_free(input_buffer);
    fftw_free(output_buffer);
  }

  void set_input(const std::vector<f64>& vec) {
    memcpy(input_buffer, vec.data(), sizeof(f64) * vec.size());
  }

  void set_input_zeropadded(const f64* buffer, std::size_t size) {
    assert(size <= input_size);
    memcpy(input_buffer, buffer, sizeof(f64) * size);
    memset(&input_buffer[size], 0, sizeof(f64) * (input_size - size));
  }

  void set_input_zeropadded(const std::vector<f64>& vec) {
    set_input_zeropadded(&vec[0], vec.size());
  }

  void set_input_pad_with_itself(const std::vector<f64>& vec,
                                 pft::Slice slice) {
    const auto buffer = pft::pad_right(vec, slice);
    assert(buffer.size() == input_size);
    memcpy(input_buffer, buffer.data(), sizeof(f64) * buffer.size());
  }

  void execute() { fftw_execute(p); }

  auto get_output_as_array() -> fftw_complex* { return output_buffer; }

  auto get_output_as_vec() const -> std::vector<std::complex<f64>> {
    std::vector<std::complex<f64>> ret(output_size);

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
  fftw_complex* const input_buffer;
  std::size_t output_size{0};
  f64* output_buffer;

  FFTW_C2R_1D(std::size_t fft_size)
      : input_size(fft_size), input_buffer(fftw_alloc_complex(fft_size)),
        output_size(input_size), output_buffer(fftw_alloc_real(output_size)) {
    p = fftw_plan_dft_c2r_1d(input_size, input_buffer, output_buffer,
                             FFTW_ESTIMATE);
  }

  ~FFTW_C2R_1D() {
    fftw_destroy_plan(p);
    fftw_free(input_buffer);
    fftw_free(output_buffer);
  }

  void set_input_zeropadded(const fftw_complex* buffer, std::size_t size) {
    assert(size <= input_size);
    memcpy(input_buffer, buffer, sizeof(fftw_complex) * size);
    memset(&input_buffer[size], 0, sizeof(fftw_complex) * (input_size - size));
  }

  // void set_input_zeropadded(const std::vector<std::complec<f64>>& vec) {
  //   set_input_zeropadded(&vec[0], vec.size());
  // }

  void execute() { fftw_execute(p); }

  auto get_output_as_vec() const -> std::vector<f64> {
    std::vector<f64> ret(output_buffer, output_buffer + output_size);
    return ret;
  }

  auto get_normalised_output_as_vec() const -> std::vector<f64> {
    const auto ret =
        pft::map([&](const auto& x) { return x / output_size; },
                 std::vector<f64>(output_buffer, output_buffer + output_size));
    return ret;
  }

private:
  fftw_plan p;
};

auto convln(const std::vector<f64>& input, const std::vector<f64>& kernel)
    -> std::vector<f64> {
  fftw_plan p;
  const auto n = input.size();
  const auto m = kernel.size();

  const std::size_t padded_length = input.size() + kernel.size() - 1;

  // // zero pad input
  auto padded_input = pft::pad_right(input, padded_length - n);

  // fft of input
  fftw_complex* input_fft = nullptr;
  input_fft               = fftw_alloc_complex(padded_length);
  p = fftw_plan_dft_r2c_1d(padded_length, padded_input.data(), input_fft,
                           FFTW_ESTIMATE);
  fftw_execute(p);

  // zero pad kernel
  auto padded_kernel = pft::pad_right(kernel, padded_length - m);

  // fft of kernel
  fftw_complex* kernel_fft = nullptr;
  kernel_fft               = fftw_alloc_complex(padded_length);
  p = fftw_plan_dft_r2c_1d(padded_length, padded_kernel.data(), kernel_fft,
                           FFTW_ESTIMATE);
  fftw_execute(p);

  // mutliply the ffts
  fftw_complex* fixed_product = nullptr;
  fixed_product               = fftw_alloc_complex(padded_length);

  for (std::size_t i = 0; i < padded_length / 2; ++i) {
    fixed_product[i][0] =
        input_fft[i][0] * kernel_fft[i][0] - input_fft[i][1] * kernel_fft[i][1];
    fixed_product[i][1] =
        input_fft[i][0] * kernel_fft[i][1] + input_fft[i][1] * kernel_fft[i][0];
  }

  // iff of the prodyct
  f64* output = nullptr;
  output      = fftw_alloc_real(padded_length);
  p = fftw_plan_dft_c2r_1d(padded_length, fixed_product, output, FFTW_ESTIMATE);
  fftw_execute(p);

  // normalize output due to fftw scaling
  for (std::size_t i = 0; i < n; ++i) {
    output[i] /= padded_length;
  }

  std::vector<f64> ret(output, output + n);
  fftw_destroy_plan(p);
  fftw_free(input_fft);
  fftw_free(kernel_fft);
  fftw_free(fixed_product);
  fftw_free(output);
  return ret;
}
#endif
} // namespace utils
#endif // UTILS_H_
