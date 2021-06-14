#define UTILS_USE_FFTW
#include "../utils.hpp"
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <deque>

using namespace std;

constexpr double freq      = 3.0;
constexpr double amplitude = 100.0;

vector<double> fftme(const vector<double>& in) {
  utils::FFTW_R2C_1D input_fft(in.size());
  input_fft.set_input_zeropadded(in);
  input_fft.execute();

  auto out = input_fft.get_output_as_vec();

  vector<double> vec(out.size(), 0);

  for (size_t i = 0; i < out.size(); ++i) {
    vec[i] = out[i].real();
  }

  // vec = pft::normalize(vec);
  auto ret = utils::abs_complex(out);

  return ret;
}

int main() {
  TApplication app("", nullptr, nullptr);
  auto canvas = new TCanvas();
  canvas->Divide(2, 1);

  auto wave = [&](double x) {
    return amplitude * TMath::Sin(2.0 * TMath::Pi() * freq * x);
  };
  auto xs    = pft::linspace(0.0, 2.0 * TMath::Pi(), 1000, true);
  auto ys    = pft::map(wave, xs);
  auto graph = new TGraph(xs.size(), xs.data(), ys.data());
  graph->SetDrawOption("alp");

  auto vec = fftme(ys);

  const auto dt = xs[1] - xs[0];
  auto tpCount  = ys.size();
  const auto N  = ys.size() / 2 + 1;
  auto vals     = pft::arange(0, (int)tpCount);
  // auto timePeriod  = tpCount / samplingFrequency;
  auto fa          = 1.0 / dt;
  auto frequencies = pft::linspace(0.0, fa / 2.0, N, true);
  auto graph2      = new TGraph();

  pft::println(stdout, "dt=", dt, ", fa=", fa);
  for (size_t i = 0; i < N; ++i) {
    // graph2->SetPoint(i, frequencies[i], out[i][0]);
    graph2->SetPoint(i, frequencies[i], 2.0 * vec[i] / N);
  }
  graph2->SetLineColor(kRed);

  // auto xx    = pft::arange(0, 10);
  // auto paded = pft::pad(xx, 69, 2);
  // auto bu    = xx < 5;

  // auto xxx  = pft::where(bu, xx, 10 * xx);
  // auto xxxx = pft::where(xx < 5, xx, 10 * xx);

  auto peaks  = utils::find_peaks(ys, {});
  auto graph3 = new TGraph(peaks.size(), pft::take(xs, peaks).data(),
                           pft::take(ys, peaks).data());
  // for (size_t i = 0; i < peaks.size(); ++i) {
  //   // pft::println(stdout, make_pair(xs[peaks[i]], ys[peaks[i]]));
  //   graph3->SetPoint(i, xs[peaks[i]], ys[peaks[i]]);
  // }
  graph3->SetMarkerColor(kRed);
  graph3->SetMarkerStyle(23);
  graph3->SetDrawOption("ap");

  canvas->cd(1);
  auto tm1 = new TMultiGraph();
  tm1->Add(graph, "PL");
  tm1->Add(graph3, "p");
  tm1->Draw("A");
  canvas->Draw();

  canvas->cd(2);
  graph2->Draw("al");
  canvas->Draw();

  auto yyy = pft::take(ys, peaks);
  auto idy = pft::argsort<int>({69, 2930, 123, 44, 420});

  cout << "DONE\n";
  app.Run();
  return 0;
}
