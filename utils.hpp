#ifndef UTILS_H_
#define UTILS_H_
#include "pft.hpp"
#include <cmath>
#include <iostream>

using pft::Maybe;
using namespace std;

tuple<vector<double>, vector<int>, vector<int>>
peak_prominences(const vector<double>& xs, const vector<int> peaks, int wlen) {
  const size_t peaks_size = peaks.size();

  vector<double> prominences(peaks_size);
  vector<int> left_bases(peaks_size);
  vector<int> right_bases(peaks_size);

  double left_min, right_min;
  int peak, i_min, i_max, i;
  for (size_t pi = 0; pi < peaks_size; ++pi) {
    peak  = peaks[pi];
    i_min = 0;
    i_max = xs.size() - 1;

    if (!(i_min <= peak && peak <= i_max)) {
      cerr << "Oopsie\n";
      exit(1);
    }

    if (2 <= wlen) {
      i_min = max(peak - wlen / 2, i_min);
      i_max = min(peak + wlen / 2, i_max);
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
    prominences[pi] = xs[peak] - max(left_min, right_min);
  }
  return make_tuple(prominences, left_bases, right_bases);
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>>
peak_widths(const vector<double>& xs, const vector<int>& peaks,
            double rel_height) {
  const size_t peaks_size = peaks.size();
  const int wlen          = -1;
  assert(rel_height > 0);
  auto prom_data = peak_prominences(xs, peaks, wlen);
  auto [prominences, left_bases, right_bases] = prom_data;
  vector<double> widths(peaks_size);
  vector<double> width_heights(peaks_size);
  vector<double> left_ips(peaks_size);
  vector<double> right_ips(peaks_size);

  double height;
  int i_min, i_max, peak;

  for (size_t i = 0; i < peaks_size; ++i) {
    i_min = left_bases[i];
    i_max = right_bases[i];
    peak  = peaks[i];
    if (!(0 <= i_min && i_min <= peak && peak <= i_max && i_max < xs.size())) {
      cerr << "Oopsie\n";
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
  return make_tuple(widths, width_heights, left_ips, right_ips);
}

template <typename T>
tuple<vector<int>, vector<int>, vector<int>> local_maxima(vector<T> x) {
  const size_t n = x.size();
  vector<int> midpoints(n / 2, 0);
  vector<int> left_edges(n / 2, 0);
  vector<int> right_edges(n / 2, 0);
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
deque<bool> select_by_property(vector<T> p, T pmin, T pmax) {
  deque<bool> keep(p.size(), 1);

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
pair<double, double> unpack_condition_args(pair<T, T> interval,
                                           vector<double> xs,
                                           vector<int> peaks) {

  // TODO: implement unpacking for when T is a container
  auto [imin, imax] = interval;

  return make_pair(imin, imax);
}

deque<bool> select_peaks_by_distance(vector<int> peaks, vector<double> priority,
                                     double distance) {

  int peaks_size = peaks.size();

  distance = ceil(distance);
  deque<bool> keep(peaks_size, 1);

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

vector<int> find_peaks(const vector<double>& xs,
                       Maybe<pair<double, double>> height,
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
  pft::println(stdout, V);

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
  for(size_t i =0; i<n ; ++i){
    // SG(i,:) *= factorial(i);
  }

  return SG.transpoze();
}

template <typename T>
T prod(const std::vector<T>& x) {
  auto ret =
      pft::foldl(T(1), x, [](const T& a, const T& b) -> T { return a * b; });
  return ret;
}
template <typename T>
T factorial(const T& n) {
  auto ret = prod(pft::arange(T(1), n + T(1)));
  return ret;
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

#endif // UTILS_H_
