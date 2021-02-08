#define UTILS_USE_FFTW
#include "../utils.hpp"
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <complex>
#include <iostream>

using std::cout;
using std::vector;

int main(int argc, char* argv[]) {
  TApplication app("app", nullptr, nullptr);
  if (argc == 1) {
    cout << "NO FILE GIVEN\n";
    exit(1);
  }
  TGraph* graph_waveform = nullptr;

  auto input_file = new TFile(argv[1], "READ");
  auto input_tree = dynamic_cast<TTree*>(input_file->Get("t"));
  input_tree->SetBranchAddress("Waveform", &graph_waveform);
  for (size_t i = 0; i < 10; ++i) {
    input_tree->GetEvent(i % 2);
    auto gr_wf               = graph_waveform;
    const size_t n_of_points = gr_wf->GetN();
    auto xs                  = gr_wf->GetX();
    auto ys                  = gr_wf->GetY();
    std::vector<double> x_vec(xs, xs + n_of_points);
    std::vector<double> y_vec(ys, ys + n_of_points);
    auto savgol_coefs = utils::savgol_coeffs(33, 16, 16, 2, 3);
    auto test = utils::convln(y_vec, savgol_coefs);

    auto graph2 = new TGraph();
    for (size_t j = 0; j < test.size(); ++j) {
      graph2->SetPoint(j, x_vec[j], test[j]);
    }

    auto canvas = new TCanvas();
    canvas->Divide(2, 1);
    canvas->cd(1);
    gr_wf->DrawClone();
    canvas->cd(2);
    graph2->Draw();
    canvas->Draw();
  }
  cout << "Close the program\n";
  app.Run();

  return 0;
}
