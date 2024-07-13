void LoadStyle();
void SetLegend(TLegend *);

void performance_plot() {
    LoadStyle();
    gStyle->SetOptStat(0);

    const double resolution = 1.;

    TFile *fIn = new TFile("merged_outputexp.root");

    TH1D *histMassSEPM = (TH1D*) fIn->Get("invMassDist_massp");
    TH1D *histV2SEPM = (TH1D*) fIn->Get("v22DistAll_v2n");
    histV2SEPM->Scale(1. / resolution);

    TH1D *histV4SEPM = (TH1D*) fIn->Get("v24Dist_v24p");
    histV4SEPM->Scale(1. / resolution);

    int lastBin = histMassSEPM->GetNbinsX();
    histMassSEPM->SetBinContent(lastBin + 1, 0); // Set overflow bin to zero

    std::cout << histV2SEPM->GetBinContent(histV2SEPM->FindBin(3.)) << std::endl;

    histMassSEPM->SetTitle("");
    histV2SEPM->SetTitle("");
    histV4SEPM->SetTitle("");

    histMassSEPM->GetXaxis()->SetRangeUser(2.5, 4.0);
    histMassSEPM->GetYaxis()->SetRangeUser(1000, 60000);
    histMassSEPM->GetYaxis()->SetLabelSize(0.08);
    histMassSEPM->GetYaxis()->SetTitle("Counts / (20 MeV/c^{2})");
    histMassSEPM->GetYaxis()->SetTitleSize(0.1);
    histMassSEPM->GetYaxis()->SetTitleOffset(0.75);

    histV2SEPM->GetXaxis()->SetRangeUser(2.5, 4.0);
    histV2SEPM->GetXaxis()->SetLabelSize(0.08);
    histV2SEPM->GetXaxis()->SetTitle("");
    histV2SEPM->GetXaxis()->SetTitleSize(0.09);
    histV2SEPM->GetYaxis()->SetRangeUser(histV2SEPM->GetMinimum()* 1.4, histV2SEPM->GetMaximum() * 1.6);
    histV2SEPM->GetYaxis()->SetLabelSize(0.07);
    histV2SEPM->GetYaxis()->CenterTitle(true);
    histV2SEPM->GetYaxis()->SetTitle("#it{v}_{2}^{#mu#mu}{2}  a.u.");
    histV2SEPM->GetYaxis()->SetTitleOffset(0.75);
    histV2SEPM->GetYaxis()->SetTitleSize(0.1);

    histV4SEPM->GetXaxis()->SetRangeUser(2.5, 4.0);
    histV4SEPM->GetXaxis()->SetLabelSize(0.08);
    histV4SEPM->GetXaxis()->SetTitle("");
    histV4SEPM->GetXaxis()->SetTitleSize(0.1);
    histV4SEPM->GetYaxis()->SetRangeUser(histV4SEPM->GetMinimum()* 1.4, histV4SEPM->GetMaximum() * 1.6);
    histV4SEPM->GetYaxis()->SetLabelSize(0.07);
    histV4SEPM->GetYaxis()->CenterTitle(true);
    histV4SEPM->GetYaxis()->SetTitle("#it{v}_{2}^{#mu#mu}{4}  a.u.");
    histV4SEPM->GetYaxis()->SetTitleOffset(0.75);
    histV4SEPM->GetYaxis()->SetTitleSize(0.1);

    histMassSEPM->SetLineColor(kBlack);
    histMassSEPM->SetMarkerStyle(20);
    histMassSEPM->SetMarkerSize(0.5);
    histV2SEPM->SetLineColor(kBlack);
    histV2SEPM->SetMarkerStyle(20);
    histV2SEPM->SetMarkerSize(0.5);
    histV4SEPM->SetMarkerSize(0.5);

    TFile *fFuncIn = new TFile("merged_outputexp.root");
    TF1 *funcMassfit = (TF1*) fFuncIn->Get("fitDistrib_massr");
    funcMassfit->SetLineColor(kRed);
    TF1 *funcMassSigBkg = (TF1*) fFuncIn->Get("fitSignal_massjpsi");
    funcMassSigBkg->SetLineColor(kGreen+3);
    TF1 *funcMassSigBkg2 = (TF1*) fFuncIn->Get("fitSignal_masspsi2s");
    funcMassSigBkg2->SetLineColor(kGreen + 1);
    TF1 *funcMassBkg = (TF1*) fFuncIn->Get("fitBG_massb");
    funcMassBkg->SetLineColor(kBlue + 1);
    TF1 *funcV2SigBkg = (TF1*) fFuncIn->Get("fitV2_v22r");
    funcV2SigBkg->SetLineColor(kRed + 1);
    TF1 *funcV2Bkg = (TF1*) fFuncIn->Get("fitBG_v22b");
    funcV2Bkg->SetLineColor(kBlue + 1);

    TF1 *funcV4SigBkg = (TF1*) fFuncIn->Get("fitV2_v24r");
    funcV4SigBkg->SetLineColor(kRed + 1);
    TF1 *funcV4Bkg = (TF1*) fFuncIn->Get("fitBG_v24b");
    funcV4Bkg->SetLineColor(kBlue + 1);

    TCanvas *canvasMassV2Pt = new TCanvas("canvasMassV2Pt", "", 800, 2400);

    canvasMassV2Pt->cd();
    gStyle->SetOptStat(0);

    TPad *pad1 = new TPad("pad1", "pad1", 0.05, 0.65, 0.95, 0.95);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    gPad->SetLogy(true);
    histMassSEPM->Draw("EP");
    funcMassBkg->Draw("SAME");
    funcMassfit->Draw("SAME");
    funcMassSigBkg->Draw("SAME");
    //funcMassSigBkg2->Draw("SAME");

    TLegend *legendHist = new TLegend(0.71, 0.4, 0.9, 0.91, " ", "brNDC");
    // TLegend *legendHist = new TLegend(0.7, 0.41, 0.8, 0.84, " ", "brNDC");
    SetLegend(legendHist);
    legendHist->SetTextSize(0.07);
    legendHist->SetHeader("Opposite-sign pairs");
    legendHist->AddEntry(histMassSEPM, "Data", "PLE");
     legendHist->AddEntry(funcMassSigBkg, "J/#psi signal", "L");
    legendHist->AddEntry(funcMassfit, "Total fit", "L");
    legendHist->AddEntry(funcMassBkg, "Background", "L");

    //legendHist->AddEntry(funcMassSigBkg2, "#psi(2S) signal", "L");
    legendHist->Draw();

    TLatex *latexTitle = new TLatex();
    latexTitle->SetTextSize(0.07);
    latexTitle->SetNDC();
    latexTitle->SetTextFont(42);
    latexTitle->DrawLatex(0.18, 0.87, "ALICE Performance, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    // latexTitle->DrawLatex(0.19, 0.73, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    // latexTitle->DrawLatex(0.19, 0.59, "0#minus50 %, 0 < #it{p}_{T} < 100 GeV/#it{c}");

    canvasMassV2Pt->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.05, 0.35, 0.95, 0.65);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0);
    pad2->Draw();
    pad2->cd();
    histV2SEPM->Draw("EP");
    funcV2Bkg->Draw("SAME");
    funcV2SigBkg->Draw("SAME");


    /////////////////////////////////

    TLatex *latexTitlev2 = new TLatex();
    latexTitlev2->SetTextSize(0.07);
    latexTitlev2->SetNDC();
    latexTitlev2->SetTextFont(42);
    // latexTitlev2->DrawLatex(0.19, 0.87, "ALICE Performance, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitlev2->DrawLatex(0.2, 0.89, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    latexTitlev2->DrawLatex(0.2, 0.76, "0#minus50 %, 0 < #it{p}_{T} < 100 GeV/#it{c}");

    //////////////////////////////////////

    canvasMassV2Pt->cd();
    TPad *pad3 = new TPad("pad3", "pad3", 0.05, 0.05, 0.95, 0.35);
    pad3->SetTopMargin(0);
    pad3->Draw();
    pad3->cd();
    histV4SEPM->Draw("EP");
    funcV4Bkg->Draw("SAME");
    funcV4SigBkg->Draw("SAME");

    canvasMassV2Pt->cd();
    TLatex *latexAxis = new TLatex();
    latexAxis->SetTextSize(0.035);
    latexAxis->SetNDC();
    latexAxis->SetTextFont(42);
    latexAxis->DrawLatex(0.70, 0.030, "#it{M}_{#mu#mu} (GeV/#it{c}^{2})");

    canvasMassV2Pt->Update();
    canvasMassV2Pt->SaveAs("flow_performance_exp.pdf");


}

void LoadStyle() {
    int font = 42;
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetStatBorderSize(1);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetDrawBorder(0);
    gStyle->SetTextFont(font);
    gStyle->SetStatFontSize(0.05);
    gStyle->SetStatX(0.97);
    gStyle->SetStatY(0.98);
    gStyle->SetStatH(0.03);
    gStyle->SetStatW(0.3);
    gStyle->SetTickLength(0.02, "y");
    gStyle->SetEndErrorSize(3);
    gStyle->SetLabelSize(0.05, "xyz");
    gStyle->SetLabelFont(font, "xyz");
    gStyle->SetLabelOffset(0.01, "xyz");
    gStyle->SetTitleFont(font, "xyz");
    gStyle->SetTitleOffset(0.9, "x");
    gStyle->SetTitleOffset(1.02, "y");
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetMarkerSize(1.3);
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetHatchesSpacing(0.0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(0.0);
    gStyle->SetTitleSize(0.05, "X");
    gStyle->SetTitleSize(0.045, "Y");
    gStyle->SetLabelSize(0.045, "X");
    gStyle->SetLabelSize(0.045, "Y");
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.35, "Y");
}

void SetLegend(TLegend *legend) {
    legend->SetBorderSize(0);
    legend->SetFillColor(10);
    legend->SetFillStyle(1);
    legend->SetLineStyle(0);
    legend->SetLineColor(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
}
