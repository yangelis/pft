#ifndef UTILS_H_
#define UTILS_H_
#include "pft.hpp"
#include <cmath>

using pft::Maybe;

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
peak_prominences(const std::vector<double>& xs, const std::vector<int> peaks,
                 int wlen) {
  const size_t peaks_size = peaks.size();

  std::vector<double> prominences(peaks_size);
  std::vector<int> left_bases(peaks_size);
  std::vector<int> right_bases(peaks_size);

  double left_min, right_min;
  int peak, i_min, i_max, i;
  for (size_t pi = 0; pi < peaks_size; ++pi) {
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

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
           std::vector<double>>
peak_widths(const std::vector<double>& xs, const std::vector<int>& peaks,
            double rel_height) {
  const size_t peaks_size = peaks.size();
  const int wlen          = -1;
  assert(rel_height > 0);
  auto prom_data = peak_prominences(xs, peaks, wlen);
  auto [prominences, left_bases, right_bases] = prom_data;
  std::vector<double> widths(peaks_size);
  std::vector<double> width_heights(peaks_size);
  std::vector<double> left_ips(peaks_size);
  std::vector<double> right_ips(peaks_size);

  double height;
  int i_min, i_max, peak;

  for (size_t i = 0; i < peaks_size; ++i) {
    i_min = left_bases[i];
    i_max = right_bases[i];
    peak  = peaks[i];
    if (!(0 <= i_min && i_min <= peak && peak <= i_max && i_max < xs.size())) {
      fprintf(stderr, "Oopsie\n");
      exit(1);
    }
    height           = xs[peak] - prominences[i] * rel_height;
    width_heights[i] = height;

    // intersection point on left side
    auto p = peak;
    while (i_min < p && height < xs[p]) {
      --p;
    }
    auto left_ip = (double)p;

    if (xs[p] < height) {
      left_ip += (height - xs[p]) / (xs[p + 1] - xs[p]);
    }

    // intersection point on right side
    p = peak;
    while (p < i_max && height < xs[p]) {
      ++p;
    }

    auto right_ip = (double)p;
    if (xs[p] < height) {
      right_ip -= (height - xs[p]) / (xs[p - 1] - xs[p]);
    }
    widths[i]    = right_ip - left_ip;
    left_ips[i]  = left_ip;
    right_ips[i] = right_ip;
  }
  return std::make_tuple(widths, width_heights, left_ips, right_ips);
}

template <typename T>
std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
local_maxima(std::vector<T> x) {
  const size_t n = x.size();
  std::vector<int> midpoints(n / 2, 0);
  std::vector<int> left_edges(n / 2, 0);
  std::vector<int> right_edges(n / 2, 0);
  size_t m = 0; // index pointer to the end

  size_t i   = 1;
  auto i_max = n - 1;

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
std::deque<bool> select_by_property(std::vector<T> p, T pmin, T pmax) {
  std::deque<bool> keep(p.size(), 1);

  // if (pmin.has_value) {
  for (size_t i = 0; i < p.size(); ++i) {
    keep[i] &= pmin <= p[i];
  }
  // }

  // if (pmax.has_value) {
  for (size_t i = 0; i < p.size(); ++i) {
    keep[i] &= p[i] <= pmax;
  }
  // }
  return keep;
}

template <typename T>
std::pair<double, double> unpack_condition_args(std::pair<T, T> interval,
                                                std::vector<double> xs,
                                                std::vector<int> peaks) {

  // TODO: implement unpacking for when T is a container
  auto [imin, imax] = interval;

  return std::make_pair(imin, imax);
}

std::deque<bool> select_peaks_by_distance(std::vector<int> peaks,
                                          std::vector<double> priority,
                                          double distance) {

  int peaks_size = peaks.size();

  distance = ceil(distance);
  std::deque<bool> keep(peaks_size, 1);

  auto priority_to_position = pft::argsort(priority);
  int j                     = 0;
  for (int i = peaks_size - 1; i > -1; --i) {
    j = priority_to_position[i];
    if (keep[j] == 0) {
      continue;
    }

    auto k = j - 1;
    while (0 <= k && (peaks[j] - peaks[k]) < distance) {
      keep[k] = 0;
      --k;
    }
    k = j + 1;
    while (k < peaks_size && (peaks[k] - peaks[j]) < distance) {
      keep[k] = 0;
      ++k;
    }
  }

  return keep;
}

std::vector<int> find_peaks(const std::vector<double>& xs,
                            Maybe<std::pair<double, double>> height,
                            double distance = 0.0) {
  auto [peaks, l_edges, r_edges] = local_maxima(xs);

  if (height.has_value) {
    auto peak_heights = pft::take(xs, peaks);
    auto [hmin, hmax] = unpack_condition_args(height.unwrap, xs, peaks);
    auto keep         = select_by_property(peak_heights, hmin, hmax);

    peaks = pft::take(peaks, keep);
  }
  // distance between peaks
  if (distance != 0.0) {
    auto keep = select_peaks_by_distance(peaks, pft::take(xs, peaks), distance);
    // TODO: properly test that
    peaks = pft::take(peaks, keep);
  }

  return peaks;
}

template <typename T>
T prod(const std::vector<T>& x) {
  auto ret =
      pft::foldl(T(1), x, [](const T& a, const T& b) -> T { return a * b; });
  return ret;
}

double factorial(const i64& n) {
  static double table[171];
  static bool init = true;
  if (init) {
    init     = false;
    table[0] = 1;
    for (size_t i = 1; i < 171; ++i) {
      table[i] = (double)i * table[i - 1];
    }
  }
  return table[n];
}

template <typename T>
pft::Matrix<T> vandermonde(i64 halfWindow, i64 polyDeg) {
  assert(halfWindow >= 0);
  assert(polyDeg >= 0);
  // arange return an exclusive range
  auto x = pft::arange<T>(-halfWindow, halfWindow + 1);

  size_t n = polyDeg + 1;
  auto m   = x.size();

  pft::Matrix<T> V(m, n);

  for (size_t i = 0; i < m; ++i) {
    V(i, 0) = T(1);
  }

  for (size_t j = 1; j < n; ++j) {
    for (size_t i = 0; i < m; ++i) {
      V(i, j) = x[i] * V(i, j - 1);
    }
  }

  return V;
}

template <typename T>
pft::Matrix<T> SG(i64 halfWindow, i64 polyDeg) {
  assert(2 * halfWindow > polyDeg);

  auto V = vandermonde<T>(halfWindow, polyDeg);

  // QR decomposition
  auto [Q, R] = qr(V);

  auto SG = R / Q.transpose();

  auto n = SG.rows;
  auto m = SG.cols;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      SG(i, j) *= factorial(i);
    }
  }

  return SG.transpoze();
}

template <typename T>
std::vector<T> vmadd(const std::vector<T>& a, const std::vector<T>& b, T s) {
  const auto n = a.size();
  std::vector<T> ret;
  ret.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    ret.push_back(a[i] + s * b[i]);
  }

  return ret;
}

template <typename T>
pft::Matrix<T> compute_householder_factor(const std::vector<T>& v) {
  const auto n = v.size();
  pft::Matrix<T> ret(n, n);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      ret(i, j) = -2 * v[i] * v[j];
    }
  }
  for (size_t i = 0; i < n; ++i) {
    ret(i, i) += 1;
  }
  return ret;
}

template <typename T>
T vnorm(const std::vector<T>& x) {
  // T sum2 = 0;
  // for(size_t i =0;i<x.size();++i){
  //   sum2 += x[i]*x[i];
  // }

  // TODO: This is hackish
  auto sum2 =
      pft::foldl(T(), x, [](const T& a, const T& b) -> T { return a + b * b; });
  return std::sqrt(sum2);
}

template <typename T>
std::vector<T> vdiv(const std::vector<T>& x, T d) {
  auto ret = pft::map([&d](const T& a) -> T { return a / d; }, x);
  return ret;
}

template <typename T>
std::pair<pft::Matrix<T>, pft::Matrix<T>>
QRDecomposition(const pft::Matrix<T>& mat) {
  pft::Matrix<T> Q;
  pft::Matrix<T> R;

  const size_t m = mat.rows;
  const size_t n = mat.cols;

  // vector if factors Q1. Q2.... Qm
  std::vector<pft::Matrix<T>> qv(m);

  pft::Matrix<T> z = mat;
  pft::Matrix<T> z1;
  for (size_t k = 0; k < n && k < m - 1; ++k) {
    std::vector<T> e(m, 0), x(m, 0);
    double a;

    z1 = z.minor(k);
    z  = z1;

    x = z.getColumn(k);
    a = vnorm(x);

    if (mat(k, k) > 0) {
      a = -a;
    }

    for (size_t i = 0; i < m; ++i) {
      e[i] = (i == k) ? 1 : 0;
    }

    e     = vmadd(x, e, a);
    e     = vdiv(e, vnorm(e));
    qv[k] = compute_householder_factor(e);

    z1 = qv[k].mult(z);
    z  = z1;
  }
  Q = qv[0];

  for (size_t i = 1; i < n && i < m - 1; ++i) {
    z1 = qv[i].mult(Q);
    Q  = z1;
  }
  R = Q.mult(mat);

  return {Q.transpose(), R};
}

struct LUdcmp {
  using Mat_t = pft::Matrix<f64>;
  size_t n;
  Mat_t lu;
  std::vector<int> indx;
  f64 d;
  LUdcmp(const Mat_t& a);
  ~LUdcmp(){};
  std::vector<f64> solve(const std::vector<f64>& b);
  void solve(Mat_t& b, Mat_t& x);
  void inverse(Mat_t& ainv);
  // f64 det();
};

LUdcmp::LUdcmp(const Mat_t& a) : n(a.rows), lu(a), indx(n) {

  constexpr f64 tiny = 1.0e-40;
  size_t i, imax, j, k;
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
    indx[k] = imax;
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

std::vector<f64> LUdcmp::solve(const std::vector<f64>& b) {
  std::vector<f64> x(b);

  size_t i, ii = 0, ip, j;
  f64 sum;
  if (b.size() != n || x.size() != n) {
    fprintf(stderr, "Bad sizes\n");
    exit(1);
  }

  for (i = 0; i < n; ++i) {
    ip    = indx[i];
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

  for (int i = n - 1; i >= 0; --i) {
    sum = x[i];
    for (j = i + 1; j < n; ++j) {
      sum -= lu(i, j) * x[j];
    }
    x[i] = sum / lu(i, i);
  }

  return x;
}

void LUdcmp::solve(Mat_t& b, Mat_t& x) {
  size_t i, m = b.cols;
  std::vector<f64> xx(n);
  for (size_t j = 0; j < m; ++j) {
    for (i = 0; i < n; ++i) {
      xx[i] = b(i, j);
    }
    xx = solve(xx);
    for (i = 0; i < n; ++i) {
      x(i, j) = xx[i];
    }
  }
}

void LUdcmp::inverse(Mat_t& ainv) {
  ainv = Mat_t(n, n);
  ainv.diagonal();
  solve(ainv, ainv);
}

// Savitzky-Golay Coeffs for 1-D filter
// np: n points, must be an odd number
// nl, nr: n points to left and right
// ld: order of the derivative, default = 0 (no derivation)
// m: order of the polynomial
std::vector<f64> savgol_coeffs(const int np, const int nl, const int nr,
                               const int ld, const int m) {

  if (np < nl + nr + 1 || nr < 0 || nr < 0 || ld > m || nl + nr < m) {
    fprintf(stderr, "Bad arguments\n");
    exit(1);
  }

  int k, mm, imj, kk;
  f64 sum, fac;
  std::vector<int> indx(m + 1);
  pft::Matrix<f64> a(m + 1, m + 1);

  for (int ipj = 0; ipj <= (m << 1); ++ipj) {
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
  LUdcmp alud(a);
  std::vector<f64> b(m + 1, 0.0);
  b[ld] = 1.0;
  b     = alud.solve(b);

  size_t i = 0;
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

#endif // UTILS_H_
