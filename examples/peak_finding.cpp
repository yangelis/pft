#define UTILS_USE_FFTW
#include "../utils.hpp"
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <iostream>

using namespace std;
using pft::Maybe;

int main() {
  TApplication app("", nullptr, nullptr);

  auto func = new TF1("func", "sin(x)/x", -1.0, 10);
  // func->Draw();

  auto xs      = pft::arange(-0.5, 10.0, 0.01);
  const auto n = xs.size();
  vector<double> ys(n);
  for (size_t i = 0; i < n; ++i) {
    ys[i] = func->Eval(xs[i]);
  }

  auto graph = new TGraph(ys.size(), xs.data(), ys.data());

  Maybe<pair<double, double>> heights = {true, {0.02, 1.0}};
  auto peaks                          = utils::find_peaks(ys, heights);
  auto peaks_properties               = utils::peak_widths(ys, peaks, 0.90);
  // auto [widths, widths_heights, left_p, right_p] = peak_widths(ys, peaks, 0.90);

  auto graph_peaks = new TGraph();
  for (size_t i = 0; i < peaks.size(); ++i) {
    graph_peaks->SetPoint(i, xs[peaks[i]], ys[peaks[i]]);
  }
  graph_peaks->SetMarkerColor(kRed);
  graph_peaks->SetMarkerStyle(23);
  graph_peaks->SetDrawOption("ap");

  auto graph_peaks_width = new TGraph();
  for (size_t i = 0; i < peaks.size(); ++i) {
    graph_peaks_width->SetPoint(i, xs[peaks_properties[i].left_p],
                                peaks_properties[i].width_height);
    graph_peaks_width->SetPoint(i + peaks.size(),
                                xs[peaks_properties[i].right_p],
                                peaks_properties[i].width_height);
  }
  graph_peaks_width->SetMarkerColor(kBlue);
  graph_peaks_width->SetMarkerStyle(24);
  graph_peaks_width->SetDrawOption("ap");

  auto tm1    = new TMultiGraph();
  auto canvas = new TCanvas();
  tm1->Add(graph, "PL");
  tm1->Add(graph_peaks, "p");
  tm1->Add(graph_peaks_width, "p");
  tm1->Draw("A");
  canvas->Draw();

  cout << "DONE\n";
  app.Run();
  return 0;
}
