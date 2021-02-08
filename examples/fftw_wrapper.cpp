#define UTILS_USE_FFTW
#include "../utils.hpp"
#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
#include <TMultiGraph.h>

using RDF = ROOT::RDataFrame;

int main(int argc, char* argv[]) {
  if (argc == 1) {
    pft::panic("NO FILE GIVEN");
  }
  TApplication app("", nullptr, nullptr);

  RDF df("t", argv[1]);

  auto df_aux = df.Define("xs",
                          [](TGraph& graph) {
                            auto points  = graph.GetX();
                            const auto n = graph.GetN();
                            return std::vector<double>(points, points + n);
                          },
                          {"Waveform"})
                    .Define("ys",
                            [](TGraph& graph) {
                              auto points  = graph.GetY();
                              const auto n = graph.GetN();
                              return std::vector<double>(points, points + n);
                            },
                            {"Waveform"});

  const auto data_xs = df_aux.Take<std::vector<double>>("xs").GetValue();
  const auto data_ys = df_aux.Take<std::vector<double>>("ys").GetValue();

  const auto& xs = data_xs[0];
  const auto& ys = data_ys[0];

  auto canvas = new TCanvas();
  canvas->Divide(2, 1);

  auto graph = new TGraph(xs.size(), xs.data(), ys.data());
  graph->SetDrawOption("alp");

  utils::FFTW_R2C_1D r2c_fft(ys.size());
  r2c_fft.set_input_zeropadded(ys);
  r2c_fft.execute();

  auto ys_fft = r2c_fft.get_output_as_array();

  // utils::FFTW_C2R_1D c2r_fft(r2c_fft.output_size);
  utils::FFTW_C2R_1D c2r_fft(ys.size());
  c2r_fft.set_input_zeropadded(ys_fft, r2c_fft.output_size);
  c2r_fft.execute();
  // auto ys_ifft = c2r_fft.get_output_as_vec();
  auto ys_ifft = c2r_fft.get_normalised_output_as_vec();

  FILE* fd = fopen("fft_testing.txt", "w");
  for (size_t i = 0; i < xs.size(); ++i) {
    fprintf(fd, "%4.15f\n", ys[i] - ys_ifft[i]);
  }

  auto graph2 = new TGraph(ys_ifft.size(), xs.data(), ys_ifft.data());

  canvas->cd(1);
  graph->Draw();
  canvas->Draw();

  canvas->cd(2);
  graph2->Draw();
  canvas->Draw();

  pft::println(stdout, "DONE");
  app.Run();
  return 0;
}
