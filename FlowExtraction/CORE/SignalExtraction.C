// Author: Victor Valencia
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <TStyle.h>
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TFitResult.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include <fstream>

#include "InitializeFunctions.C"
#include "SetRangeAndNameTest.C"

float VDN2Value(float cor2d, float cor2)
{
    if (cor2 < 0)
        return -2;
    return cor2d / TMath::Sqrt(cor2);
};
float VDN2Error(float cor2d, float cor2de, float cor2, float cor2e)
{
    float sqrtv = cor2de * cor2de / cor2 + 0.25 * cor2d * cor2d * cor2e * cor2e / (cor2 * cor2 * cor2);
    if (sqrtv < 0)
        return 0;
    return TMath::Sqrt(sqrtv);
};

std::string GetFunctionName(Efunction func)
{
    switch (func)
    {
    case kVWGQuadratic:
        return "VWG Quadratic";
    case kPol2Pol3:
        return "Pol2 Pol3";
    case kDoubleExp:
        return "Double Exp";
    case kExp:
        return "Exp";
    case kVWG:
        return "VWG";
    case kCBExtended:
        return "CB Extended";
    case kNA60:
        return "NA60";
    case kv2Pol2vsMassCB2VWG:
        return "Mass_VWG_CB2-V2_POL2";
    case kv2Pol2vsMassNA60QVWG:
        return "Mass_QVWG_NA60-V2_POL2";
    case kv2Pol2vsMassNA60VWG:
        return "Mass_VWG_NA60-V2_POL2";
    case kv2Pol2vsMassNA60Pol2Pol3:
        return "Mass_POL2POL3_NA60-V2_POL2";
    case kv2Pol2vsMassCB2Pol2Pol3:
        return "Mass_POL2POL3_CB2-V2_POL2";
    case kv2Pol2vsMassCB2QVWG:
        return "Mass_QVWG_CB2-V2_POL2";
    //
    case kv2Pol3vsMassCB2VWG:
        return "Mass_VWG_CB2-V2_POL3";
    case kv2Pol3vsMassNA60QVWG:
        return "Mass_QVWG_NA60-V2_POL3";
    case kv2Pol3vsMassNA60VWG:
        return "Mass_VWG_NA60-V2_POL3";
    case kv2Pol3vsMassNA60Pol2Pol3:
        return "Mass_POL2POL3_NA60-V2_POL3";
    case kv2Pol3vsMassCB2Pol2Pol3:
        return "Mass_POL2POL3_CB2-V2_POL3";
    case kv2Pol3vsMassCB2QVWG:
        return "Mass_QVWG_CB2-V2_POL3";
    //
    case kv2PolExpvsMassCB2VWG:
        return "Mass_VWG_CB2-V2_PolExp";
    case kv2PolExpvsMassNA60QVWG:
        return "Mass_QVWG_NA60-V2_PolExp";
    case kv2PolExpvsMassNA60VWG:
        return "Mass_VWG_NA60-V2_PolExp";
    case kv2PolExpvsMassNA60Pol2Pol3:
        return "Mass_Pol2Pol3Exp_NA60-V2_PolExp";
    case kv2PolExpvsMassCB2Pol2Pol3:
        return "Mass_Pol2Pol3Exp_CB2-V2_PolExp";
    case kv2PolExpvsMassCB2QVWG:
        return "Mass_QVWG_CB2-V2_PolExp";
    }
}

// TH1* GetInvMassHisto(TFile* file, Double_t minMas, Double_t maxMass, Double_t minPt, Double_t maxPt, Double_t minY, Double_t maxY,  Int_t minCent, Int_t maxCent, ECentralityEstimator estimator);
// TH1* GetInvMassHisto(TFile* file, Double_t minMas, Double_t maxMass, Double_t minPt, Double_t maxPt,  Int_t minCent, Int_t maxCent);
TH1 *GetInvMassHisto2(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent);
TH1 *GetInvMassHisto3(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent);
TH1 *GetInvMassHisto4(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent);
TH1 *Getv24Histo(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent);
TH1 *GetvepHisto(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent);
TH1 *Getv22Histo(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent);

bool isTailsParameter(TString nameTails);

void SaveFitResults(const std::string &title, const std::string &filename, Double_t v2jpsi, Double_t Errorv2jpsi);

void LoadData(TChain *fChain, const char *TreeName, const char *FileName);

void CreateBins(double *axis, double min, double max, int Nbins);

void writeToFile(const char *filename, TH1F *hist);

bool getTChain(TChain *fChain, string FileName, string TreeName);

#define JPSI_MASS 3.096916
#define PSI2S_MASS 3.686109

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
std::vector<Double_t> SignalExtraction(const char *file, Float_t minY, Float_t maxY, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent, Float_t minMass, Float_t maxMass, Efunction fBackGround, Efunction fSignal, Efunction fBackGroundv2, Efunction fSignalv2, Etails pTails, Double_t minFit, Double_t maxFit, Bool_t isDisplayed, Bool_t isSaved)
{

    ifstream file2("input_extraction_txt/input.txt");
    string line2;

    vector<pair<string, vector<double>>> data;

    while (getline(file2, line2))
    {
        string str = line2;
        getline(file2, line2);
        stringstream ss(line2);
        vector<double> params;
        string p;
        while (getline(ss, p, ','))
        {
            params.push_back(stod(p));
        }
        data.push_back(make_pair(str, params));
    }
    for (const auto &pair : data)
    {
        addV2MassModel(V2MassModelDictionary[pair.first], pair.second);
    }

    std::vector<Double_t> results;
    Int_t nb_bg = GetNPar(fBackGround);
    Int_t nb_sig = GetNPar(fSignal);

    TString nameTest = SetNameTest(fBackGround, fSignal, pTails, minFit, maxFit);
    TString rangeName = SetRangeName(minY, maxY, minPt, maxPt, minCent, maxCent);
    TString tailsName = SetTailsName(pTails, fSignal, fBackGround, minY, maxY, minPt, maxPt, minFit, maxFit);
    TString rangeNameFull;
    rangeNameFull.Form("rapidity%.1f-%.1f_pT%.2f-%.2f_centrality%1.f-%1.f", minY, maxY, minPt, maxPt, minCent, maxCent);
    if (!isTailsParameter(tailsName))
        return results;

    TH1 *invMassDist;

    invMassDist = GetInvMassHisto2(file, minMass, maxMass, minPt, maxPt, minCent, maxCent);

    // Fit Background
    TF1 *BGFunction = BackGroundFunction(fBackGround, 0., 10.);

    reject = kTRUE;
    Double_t n_BG = invMassDist->GetBinContent(invMassDist->GetXaxis()->FindBin(minFit));
    BGFunction->SetParameter(0, n_BG);

    int secur1 = 0;
    TFitResultPtr fitStatus1;
    Double_t chi21 = 3;
    int covMatrixStatus1 = 2;

    do
    {
        if (fBackGround == kDoubleExp)
            fitStatus1 = invMassDist->Fit("fitBG", "SN", "", 2.2, maxFit); // eventmixing : pas likelyhood
        else
            fitStatus1 = invMassDist->Fit("fitBG", "LSN", "", 2.2, maxFit); // Direct fit
        chi21 = BGFunction->GetChisquare() / BGFunction->GetNDF();
        covMatrixStatus1 = fitStatus1->CovMatrixStatus();

        secur1++;
        if (secur1 > 3)
            break;
        cout << "FitStatus = " << fitStatus1 << "     chi2/NDF = " << chi21 << "     cov matrix status = " << covMatrixStatus1 << endl;
    } while (fitStatus1 != 0 || chi21 > 2.5 || covMatrixStatus1 != 3);

    reject = kFALSE;

    // Fit Invariant mass distribution
    TF1 *fitFunction = DistributionFunction(fBackGround, fSignal, tailsName, 0., 10.);
    Double_t n_jpsi = invMassDist->GetBinContent(invMassDist->GetXaxis()->FindBin(JPSI_MASS));
    n_jpsi -= BGFunction->Eval(JPSI_MASS);

    for (int i = 0; i < nb_bg; i++)
    {
        fitFunction->SetParameter(i, BGFunction->GetParameter(i));
    }
    fitFunction->SetParameter(nb_bg, n_jpsi);

    int secur2 = 0;
    TFitResultPtr fitStatus;
    Double_t chi22 = 3;
    int covMatrixStatus2 = 2;

    do
    {
        if (fBackGround == kDoubleExp)
            fitStatus = invMassDist->Fit("fitDistrib", "SN", "", minFit, maxFit); // eventmixing : pas likelyhood
        else
            fitStatus = invMassDist->Fit("fitDistrib", "LSN", "", minFit, maxFit); // Direct fit
        chi22 = fitFunction->GetChisquare() / fitFunction->GetNDF();
        covMatrixStatus2 = fitStatus->CovMatrixStatus();
        cout << "FitStatus = " << fitStatus << "     chi2/NDF = " << chi22 << "     cov matrix status = " << covMatrixStatus2 << endl;

        secur2++;
        if (secur2 > 30)
        {
            cout << "______________________________" << endl;
            cout << " The fit has not converged " << endl;
            cout << "______________________________" << endl;
            break;
        }
    } while (fitStatus != 0 || chi22 > 2.5 || covMatrixStatus2 != 3);

    // Draw Function
    TF1 *JPsiFunction = SignalFunction(fSignal, tailsName, kJPsi, 0., 10.);
    TF1 *Psi2SFunction = SignalFunction(fSignal, tailsName, kPsi2S, 0., 10.);
    Double_t para[nb_bg + nb_sig + 1];

    fitFunction->GetParameters(para);
    BGFunction->SetParameters(para);
    JPsiFunction->SetParameters(&(para[nb_bg]));
    Psi2SFunction->SetParameters(&(para[nb_bg]));
    Psi2SFunction->SetParameter(0, para[nb_bg + nb_sig]);
    Psi2SFunction->SetParameter(1, para[nb_bg + 1] + 3.686109 - 3.096916);
    Psi2SFunction->SetParameter(2, para[nb_bg + 2] * 1.05);

    Double_t paraJPsi[nb_sig];
    JPsiFunction->GetParameters(paraJPsi);
    Double_t paraPsi2S[nb_sig];
    Psi2SFunction->GetParameters(paraPsi2S);

    // Compute Number of JPsi
    double dx = invMassDist->GetBinWidth(1);
    Double_t N_JPsi = JPsiFunction->Integral(0, 10) / dx;
    Double_t mean_JPsi = fitFunction->GetParameter(nb_bg + 1);
    Double_t Errmean_JPsi = fitFunction->GetParError(nb_bg + 1);
    Double_t sigma_JPsi = fitFunction->GetParameter(nb_bg + 2);
    Double_t Errsigma_JPsi = fitFunction->GetParError(nb_bg + 2);

    Double_t covmat[nb_sig][nb_sig];
    for (int k = 0; k < nb_sig; k++)
    {
        for (int t = 0; t < nb_sig; t++)
        {
            covmat[k][t] = (fitStatus->GetCovarianceMatrix())(k + nb_bg, t + nb_bg);
        }
    }

    Double_t Err_JPsi = JPsiFunction->IntegralError(0, 10, &paraJPsi[0], &covmat[0][0], 0.05) / (dx);
    cout << "Number of J/psi " << N_JPsi << " with associated error " << Err_JPsi << endl;

    // Compute Number of Psi2S
    Double_t N_Psi2S = Psi2SFunction->Integral(0, 9) / dx;

    Double_t covmatPsi2S[nb_sig][nb_sig];
    int m = 0;
    int n = 0;
    for (int k = 0; k < nb_sig; k++)
    {
        for (int t = 0; t < nb_sig; t++)
        {
            if (k == 0)
                m = nb_sig;
            else
                m = k;
            if (t == 0)
                n = nb_sig;
            else
                n = t;

            covmatPsi2S[k][t] = (fitStatus->GetCovarianceMatrix())(m + nb_bg, n + nb_bg);
        }
    }

    Double_t Err_Psi2S = Psi2SFunction->IntegralError(0, 9., &paraPsi2S[0], &covmatPsi2S[0][0], 0.05) / (dx);
    cout << "Number of Psi2S " << N_Psi2S << " with associated error " << Err_Psi2S << endl;

    // Set Results
    results.push_back(N_JPsi);
    results.push_back(Err_JPsi);
    results.push_back(N_Psi2S);
    results.push_back(Err_Psi2S);
    results.push_back(mean_JPsi);
    results.push_back(Errmean_JPsi);
    results.push_back(sigma_JPsi);
    results.push_back(Errsigma_JPsi);
    results.push_back(chi22);
    results.push_back(0.);
    results.push_back(fitStatus);
    results.push_back(covMatrixStatus2);

    // Display
    if (isDisplayed == kTRUE)
    {

        std::string massTitle = Form("Mass_%s[Pt%.1f-%.1f]",
                                     nameTest.Data(),
                                     minPt, maxPt); // Adds minPt, maxPt to the title

        // TCanvas *cInvMassDist = new TCanvas(Form("%s_%s", rangeNameFull.Data(), nameTest.Data()), Form("Range %s for test %s", rangeNameFull.Data(), nameTest.Data()));
        TCanvas *cInvMassDist = new TCanvas("cImvMassDist", "c", 700, 700);
        cInvMassDist->SetRightMargin(0.03);
        cInvMassDist->SetTopMargin(0.08);

        gStyle->SetOptFit(1111);
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
        gStyle->SetHistLineColor(1);

        invMassDist->SetMarkerStyle(kFullCircle);
        invMassDist->SetMarkerSize(0.5);
        invMassDist->GetYaxis()->SetRangeUser(10,1000000);
        gPad->SetLogy();
        invMassDist->Draw("P");
        invMassDist->GetXaxis()->SetRangeUser(2, 5);
        fitFunction->SetLineColor(kRed);
        fitFunction->SetRange(minFit, maxFit);
        fitFunction->SetNpx(10000);
        fitFunction->Draw("SAME");

        BGFunction->SetLineColor(kBlue);
        BGFunction->SetLineStyle(7);
        BGFunction->SetRange(minFit, maxFit);
        BGFunction->Draw("SAME");
        JPsiFunction->SetLineColor(kRed);
        JPsiFunction->SetNpx(10000);
        JPsiFunction->Draw("SAME");
        Psi2SFunction->SetLineColor(kGreen + 2);
        Psi2SFunction->SetNpx(10000);
        Psi2SFunction->Draw("SAME");

        // TLatex *text = new TLatex();
        // text->SetNDC();
        // text->SetTextAlign(12);
        // text->SetTextFont(43);
        // text->SetTextSize(19);
        // text->DrawLatex(0.7, 0.85, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
        // TLatex *text2 = new TLatex();
        // text2->SetNDC();
        // text2->SetTextAlign(12);
        // text2->SetTextFont(43);
        // text2->SetTextSize(17);
        // text2->DrawLatex(0.7, 0.81, Form("Centrality %.1f-%.1f%%", minCent, maxCent));
        // text2->DrawLatex(0.7, 0.78, "2.5 < #it{y} < 4");
        // text2->DrawLatex(0.7, 0.74, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", minPt, maxPt));
        // // text2->DrawLatex(0.7, 0.74, Form("#it{p}_{T} < %.1f GeV/#it{c}", maxPt));
        // TLatex *text3 = new TLatex();
        // text3->SetNDC();
        // text3->SetTextAlign(12);
        // text3->SetTextFont(43);
        // text3->SetTextSize(15);
        // text3->DrawLatex(0.7, 0.68, Form("N_{J/#psi} = %i #pm %i", TMath::Nint(N_JPsi), TMath::Nint(Err_JPsi)));
        // for (int np = 1; np < 3; np++)
        // {
        //     text3->DrawLatex(0.7, 0.68 - np * 0.04, Form("%s = %.3f #pm %.3f", fitFunction->GetParName(np + nb_bg), fitFunction->GetParameter(np + nb_bg), fitFunction->GetParError(np + nb_bg)));
        // }
        // text3->DrawLatex(0.7, 0.68 - 3 * 0.04, Form("#chi^{2}/NDF = %.3f", chi22));

        // Create a TLatex object using a pointer
        // TLatex *latexMASS = new TLatex();
        // latexMASS->SetNDC();
        // latexMASS->SetTextAlign(12);
        // latexMASS->SetTextFont(43);
        // latexMASS->SetTextSize(17);                                                             // Size in pixel height
        // latexMASS->DrawLatex(0.23, 0.55, "ALICE Performance, Pb-Pb #sqrt{s_{NN}} = 5.36 TeV "); // Draw the text at the specified position

        // // Create a TLatex object using a pointer
        // TLatex *latexJPSI = new TLatex();
        // latexJPSI->SetNDC();
        // latexJPSI->SetTextAlign(12);
        // latexJPSI->SetTextFont(43);
        // latexJPSI->SetTextSize(17);                                                       // Size in pixel height
        // latexJPSI->DrawLatex(0.63, 0.55, "Inclusive J/#psi #rightarrow #mu^{+}#mu^{-} "); // Draw the text at the specified position
        // latexJPSI->DrawLatex(0.63, 0.50, " MCH + MID ");                                  // Draw the text at the specified position
        // latexJPSI->DrawLatex(0.63, 0.45, " 2.5 < |#it{y}^{#mu#mu}| < 4 ");                // Draw the text at the specified position
        // latexJPSI->DrawLatex(0.63, 0.40, "Centrality 0-50%");                             // Draw the text at the specified position
        // latexJPSI->DrawLatex(0.63, 0.35, Form("%.1f < p^{#mu#mu}_{T} < %.1f GeV/c", minPt, maxPt));

        ///////////////////////////////

        TLatex *textMass = new TLatex();
        textMass->SetNDC();
        textMass->SetTextAlign(12);
        textMass->SetTextFont(43);
        textMass->SetTextSize(17); // Size in pixel height
        textMass->DrawLatex(0.2, 0.86, "ALICE Performance, Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
        textMass->DrawLatex(0.2, 0.83, "Inclusive J/#psi #rightarrow #mu^{+}#mu^{-} "); // Draw the text at the specified position
        textMass->DrawLatex(0.2, 0.80, "MCH + MID");
        textMass->DrawLatex(0.2, 0.77, "Centrality 0-50%");
        // textMass->DrawLatex(0.2, 0.76, Form("Centrality %.1f-%.1f%%", minCent, maxCent));
        textMass->DrawLatex(0.2, 0.74, "2.5 < |#it{y}^{#mu#mu}| < 4");
        // textMass->DrawLatex(0.7, 0.745, Form("#it{p}_{T} < %.1f GeV/#it{c}", maxPt));
        textMass->DrawLatex(0.2, 0.71, Form("%.1f < #it{p}^{#mu#mu}_{T} < %.1f GeV/c", minPt, maxPt));

        /////////////////////////////////

        TLegend *legend_lines = new TLegend(0.62, 0.6, 0.8, 0.75); // Define the legend position and size
        legend_lines->SetBorderSize(0);
        legend_lines->SetFillColor(0);
        // legend_linesV2->SetTextSize(0.03);
        legend_lines->SetTextAlign(12);
        legend_lines->SetTextFont(43);
        legend_lines->SetTextSize(17);
        legend_lines->AddEntry(invMassDist, "Data", "ep");
        legend_lines->AddEntry(fitFunction, "Total fit", "l");
        legend_lines->AddEntry(BGFunction, "Background", "l");
        legend_lines->Draw("same");

        //////////////////////////////////////////////////////////////////////////////////////////////////////////// BEGINNING FITTING V2:

        TH1 *v24Dist;
        TH1 *vepDist = new TH1F();

        // v24Dist = Getv24Histo(file, minMass, maxMass, minPt, maxPt, minCent, maxCent);
        v24Dist = Getv22Histo(file, minMass, maxMass, minPt, maxPt, minCent, maxCent);
        //  vepDist = GetvepHisto(file, minMass, maxMass, minPt, maxPt, minCent, maxCent);

        v24Dist->SetMarkerStyle(kFullCircle);
        v24Dist->SetMarkerSize(0.5);
        v24Dist->SetTitle("");

        Double_t *params = fitFunction->GetParameters();

        TF1 *fitV2 = DistributionFunction(fBackGround, fSignal, fBackGroundv2, fSignalv2, minMass, maxMass);

        // Convert enum values to string names
        std::string bgName = GetFunctionName(fBackGroundv2);
        std::string signalName = GetFunctionName(fSignalv2);

        std::string title = Form("%sm[%.1f-%.1f][Pt%.1f-%.1f]",
                                 signalName.c_str(),
                                 minFit, maxFit,
                                 minPt, maxPt); // Adds minPt, maxPt to the title

        // v24Dist->SetTitle(title.c_str());

        v24Dist->Sumw2(); // SUMW2

        double maxYv2 = v24Dist->GetMaximum();
        double minYv2 = v24Dist->GetMaximum();

       // invMassDist->SetTitle(Form("Mass_%s", nameTest.Data()));
        invMassDist->GetYaxis()->SetRangeUser(1, 500000);
        invMassDist->SetTitle("");

        TCanvas *cv24Dist = new TCanvas("cv24Dist","cv24Distß", 700, 700);
        cv24Dist->SetRightMargin(0.03);
        cv24Dist->SetTopMargin(0.08);
        fitV2->GetParameters();

        ////////////////////////////////////// FITTING V2

        int securv2 = 0;
        Double_t chi2v2 = 3;

        do
        {
            v24Dist->Fit("fitV2", "R");
            chi2v2 = fitV2->GetChisquare() / fitV2->GetNDF();
            cout << "     chi2v2/NDF = " << chi2v2 << endl;

            securv2++;
            if (securv2 > 30)
            {
                cout << "______________________________" << endl;
                cout << " The fit has not converged " << endl;
                cout << "______________________________" << endl;
                break;
            }
        } while (chi2v2 > 2.5);

        ///////////////////////////////////////

        TLatex *textv2 = new TLatex();
        textv2->SetNDC();
        textv2->SetTextAlign(12);
        textv2->SetTextFont(43);
        textv2->SetTextSize(17); // Size in pixel height
        textv2->DrawLatex(0.2, 0.86, "ALICE Performance, Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
        textv2->DrawLatex(0.2, 0.83, "Inclusive J/#psi #rightarrow #mu^{+}#mu^{-} "); // Draw the text at the specified position
        textv2->DrawLatex(0.2, 0.80, "MCH + MID");
        textv2->DrawLatex(0.2, 0.77, "Centrality 0-50%");
        // textv2->DrawLatex(0.2, 0.76, Form("Centrality %.1f-%.1f%%", minCent, maxCent));
        textv2->DrawLatex(0.2, 0.74, "2.5 < |#it{y}^{#mu#mu}| < 4");
        // textv2->DrawLatex(0.7, 0.745, Form("#it{p}_{T} < %.1f GeV/#it{c}", maxPt));
        textv2->DrawLatex(0.2, 0.71, Form("%.1f < #it{p}^{#mu#mu}_{T} < %.1f GeV/c", minPt, maxPt));

        Double_t v2jpsi = fitV2->GetParameter(GetNPar(fSignalv2) - 1);
        Double_t Errorv2jpsi = fitV2->GetParError(GetNPar(fSignalv2) - 1);

        // textv2->DrawLatex(0.2, 0.61, Form("#it{v}_{2}^{J/#psi} = %.3f #pm %.3f", v2jpsi, Errorv2jpsi));
        // textv2->DrawLatex(0.2, 0.60, Form("#chi^{2}/NDF = %.3f ", chi2v2));

        fitV2->SetLineColor(kRed);


        ///////////////////////////////////////////////////////////////////////////////////////****************








        //////////////////////////////////////////////////////////////////////////////////////

        fitV2->Draw("same");

        //////////////

        TF1 *bkgV4 = BackGroundFunction(fBackGroundv2, minMass, maxMass);

        Double_t *paramsV4 = fitV2->GetParameters();

        if (paramsV4)
        {
            int nParams = fitV2->GetNpar();
            int nParamsSignal = GetNPar(fSignalv2); // Retrieve number of parameters once
            int nParamsBkgv2 = GetNPar(fBackGroundv2);
            int nParamsBkg = GetNPar(fBackGround);
            for (int i = 0; i < nParamsBkgv2; ++i)
            {
                int paramIndex = nParamsSignal - 1; // Assuming V2jpsi parameter is always one.
                if (paramIndex < nParams)
                {
                    bkgV4->SetParameter(i, paramsV4[paramIndex - nParamsBkgv2 + i]);
                    cout << "Parameters finalll" << i << " = " << paramsV4[paramIndex - nParamsBkgv2 + i] << endl;
                }
                else
                {
                    cout << "Error: Parameter index out of bounds: " << paramIndex << endl;
                }
            }

            bkgV4->SetLineColor(kBlue); // Changed to green
            bkgV4->SetLineStyle(2);
            bkgV4->Draw("same");

            //////////////////////////////// Add Legend:

            TLegend *legend_linesV2 = new TLegend(0.62, 0.6, 0.8, 0.75); // Define the legend position and size
            legend_linesV2->SetBorderSize(0);
            legend_linesV2->SetFillColor(0);
            // legend_linesV2->SetTextSize(0.03);
            legend_linesV2->SetTextAlign(12);
            legend_linesV2->SetTextFont(43);
            legend_linesV2->SetTextSize(17);
            legend_linesV2->AddEntry(v24Dist, "Data", "ep");
            legend_linesV2->AddEntry(fitV2, "Total fit", "l");
            legend_linesV2->AddEntry(bkgV4, "Background", "l");
            legend_linesV2->Draw("same");

            // // Create a TLatex object using a pointer
            // TLatex *latex = new TLatex();
            // latex->SetNDC();
            // latex->SetTextAlign(12);
            // latex->SetTextFont(43);
            // latex->SetTextSize(17); // Size in pixel height
            // latex->DrawLatex(0.63, 0.88, "ALICE Performance"); // Draw the text at the specified position

            ///////////////////////////////////////////////

            gPad->Update(); // Ensure the drawing is updated
        }
        else
        {
            cout << "Error: Failed to get parameters from fitV2." << endl;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////END FITTING V2:

        if (isSaved)
        {

            std::string pdfPath = "";

            cInvMassDist->SaveAs(("output_extraction_pdf/" + pdfPath + massTitle + ".pdf").c_str());
            cv24Dist->SaveAs(("output_extraction_pdf/" + pdfPath + title + ".pdf").c_str());
            std::string filename = "output_extraction_txt/output_v2jpsi_sys.txt";

            SaveFitResults(title, filename, v2jpsi, Errorv2jpsi);
        }
    }
    // analysis->Close();
    return results;
}

//___________________________________________________________________________________________________________
//__________________________________________________________________________________________________________

// Function to get invariant mass histogram using a TChain
TH1 *GetInvMassHisto2(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent)
{

    // Create a TChain
    const char *treeName = "O2rtdimuonall";
    TChain *chain = new TChain(treeName);

    // Load data into the TChain
    cout << "GETINVHISTO ??? = " << filePattern << endl;
    // LoadData(chain, treeName, filePattern);
    string tree_str = treeName;
    string file_str = filePattern;
    getTChain(chain, file_str, tree_str);

    // Set up branches
    Float_t mass = 0;
    Float_t pt, pt1, pt2 = 0;
    Float_t cent, fEta = 0;
    int fSign = -999;
    chain->SetBranchAddress("fMass", &mass);
    chain->SetBranchAddress("fPt", &pt);
    chain->SetBranchAddress("fCentFT0C", &cent);
    chain->SetBranchAddress("fPt1", &pt1);
    chain->SetBranchAddress("fPt2", &pt2);
    chain->SetBranchAddress("fSign", &fSign);
    chain->SetBranchAddress("fEta", &fEta);

    // Create the histogram
    TH1D *invMassDist = new TH1D("invMassDist", "; M_{#mu#mu} (GeV/c^{2}); Counts / (20 MeV/c^{2})", 100, minMass, maxMass);

    // Fill the histogram
    Long64_t nentries = chain->GetEntries();
    for (Long64_t i = 0; i < nentries; i++)
    {

        chain->GetEntry(i);

        if (pt1 >= 1.0 && pt2 >= 1.0)
        {
            if (fSign == 0 && pt >= minPt && pt <= maxPt && cent >= minCent && cent <= maxCent && fEta > -4. && fEta < -2.5)
            {
                invMassDist->Fill(mass);
            }
        }
    }

    for (Long64_t i = 0; i < nentries; i++)
    {

        float massH = invMassDist->GetBinContent(i + 1);
        float emassH = invMassDist->GetBinError(i + 1);
        invMassDist->SetBinContent(i + 1, massH);
        invMassDist->SetBinError(i + 1, emassH);
    }

    delete chain;

    return invMassDist; // Return the filled histogram
}

TH1 *GetInvMassHisto3(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent)
{

    //*************************************************************************************************************
    filePattern = "/Users/valencia/Desktop/AnalysisFlow/Flow/download/runs_chi/AO2D_merged_30april.root"; // PATH
    //*************************************************************************************************************

    //*************************************************************************************************************
    // const char *treeName = "O2rtdimuonall";
    const char *treeName = "O2rtdileptonflow";
    //*************************************************************************************************************
    TChain *chain = new TChain(treeName);

    // Load data into the TChain
    cout << "GETINVHISTO ??? = " << filePattern << endl;
    // LoadData(chain, treeName, filePattern);
    string tree_str = treeName;
    string file_str = filePattern;
    getTChain(chain, file_str, tree_str);

    // Set up branches
    Float_t mass = 0;
    Float_t pt, pt1, pt2 = 0;
    Float_t cent, fEta = 0;
    int fSign = -999;
    chain->SetBranchAddress("fMass", &mass);
    chain->SetBranchAddress("fPt", &pt);
    chain->SetBranchAddress("fCentFT0C", &cent);
    chain->SetBranchAddress("fSign", &fSign);
    chain->SetBranchAddress("fEta", &fEta);

    // Create the histogram
    TH1D *invMassDist = new TH1D("invMassDist", "Invariant Mass Distribution; M_{#mu#mu} (GeV/#it{c}^{2}) ; Counts per bin", 100, minMass, maxMass);
    invMassDist->Sumw2(); // Enable calculation of the sum of squares of weights

    // Fill the histogram
    Long64_t nentries = chain->GetEntries();
    for (Long64_t i = 0; i < nentries; i++)
    {

        chain->GetEntry(i);

        if (fSign == 0 && pt >= minPt && pt <= maxPt && cent >= minCent && cent <= maxCent && fEta > -4. && fEta < -2.5)
        {
            invMassDist->Fill(mass);
        }
    }

    delete chain;

    return invMassDist; // Return the filled histogram
}

bool isTailsParameter(TString nameTails)
{
    std::ifstream file("input_extraction_txt/keys_list.txt"); // Open the text file
    std::string line;

    if (file.is_open())
    {
        while (getline(file, line))
        {
            if (line == nameTails)
            {
                file.close();
                return true; // Parameter was found
            }
        }
        file.close();
    }
    else
    {
        std::cerr << "Unable to open file";
    }
    return false; // Parameter was not found
}

void LoadData(TChain *fChain, const char *TreeName, const char *FileName)
{
    cout << "Filename????? = " << FileName << endl;
    TFile *fInput = TFile::Open(FileName);
    if (!fInput || fInput->IsZombie())
    {
        std::cerr << "Error opening file: " << FileName << std::endl;
        return;
    }

    TIter keyList(fInput->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)keyList()))
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl)
        {
            continue; // Skip keys that do not correspond to a recognized class
        }
        if (!cl->InheritsFrom("TDirectoryFile"))
            continue;
        const char *dir = key->GetName();
        cout << "DIR ????? = " << dir << endl;
        if (strcmp(dir, "parentFiles") == 0)
        {
            continue; // Skip the directory named 'parentFiles'
        }
        fChain->Add(Form("%s/%s/%s", FileName, dir, TreeName));
        cout << "ADDING TCHAIN ??? = " << Form("%s/%s/%s", FileName, dir, TreeName) << endl;
    }

    fInput->Close(); // Close the file after adding the directories to the chain
}

bool getTChain(TChain *fChain, string FileName, string TreeName)
{
    cout << "[INFO] Extrating TTree " << TreeName.c_str() << endl;

    TFile *inputFile = TFile::Open(FileName.c_str());
    cout << "opening file " << FileName.c_str() << endl;
    TIter keyList(inputFile->GetListOfKeys());
    TKey *key;

    while ((key = (TKey *)keyList()))
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TDirectoryFile"))
            continue;
        string dir = key->GetName();
        cout << "dir : " << dir.c_str() << endl;
        string nb = dir;
        nb.erase(0, 9);
        cout << "file name : " << dir.c_str() << endl;
        cout << "[INFO] Adding TFile " << FileName.c_str() << dir.c_str() << endl;
        fChain->Add(Form("%s/%s/%s", FileName.c_str(), dir.c_str(), TreeName.c_str()));
    }

    inputFile->Close();

    if (!fChain)
    {
        cout << "[ERROR] fChain was not created, some input files are missing" << endl;
        return false;
    }
    return true;
};

void CreateBins(double *axis, double min, double max, int Nbins)
{
    for (int i = 0; i < Nbins; i++)
    {
        axis[i] = min + i * (max - min) / Nbins;
    }
    axis[Nbins] = max;
}

void writeToFile(const char *filename, TH1F *hist)
{
    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write the header
    outFile << "bin_range\tbin_content\tbin_error\n";

    // Write the bin contents and errors to the file
    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        double binMin = hist->GetXaxis()->GetBinLowEdge(i);
        double binMax = hist->GetXaxis()->GetBinUpEdge(i);
        double binContent = hist->GetBinContent(i);
        double binError = hist->GetBinError(i);
        outFile << binMin << " " << binMax << "\t" << binContent << "\t" << binError << std::endl;
    }

    outFile.close();
}

// Function to get v24 histogram using a TChain
TH1 *Getv24Histo(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent)
{
    // Create a TChain

    const char *treeName = "O2rtdimuonall";
    TChain *chain = new TChain(treeName);

    // Load data into the TChain
    cout << "GETINVHISTO ??? = " << filePattern << endl;
    // LoadData(chain, treeName, filePattern);
    filePattern = "/Users/valencia/Desktop/AnalysisFlow/daiki/9nov/O2DQworkflows/O2DQworkflows/AO2D_test_ref_flow.root";

    // filePattern = "/Users/valencia/Desktop/Flow/download/544124/AOD/002/AO2D.root";



    string tree_str = treeName;
    string file_str = filePattern;
    getTChain(chain, file_str, tree_str);

    // Set up branches
    int nBins = 100;
    Float_t mass = 0;
    Float_t pt, pt1, pt2, V4 = 0;
    Float_t cent = 0;
    int fSign = -999;
    chain->SetBranchAddress("fMass", &mass);
    chain->SetBranchAddress("fPt", &pt);
    chain->SetBranchAddress("fCentFT0C", &cent);
    chain->SetBranchAddress("fPt1", &pt1);
    chain->SetBranchAddress("fPt2", &pt2);
    chain->SetBranchAddress("fV4", &V4);
    chain->SetBranchAddress("fSign", &fSign);

    // Create the histogram
    TH1D *v24Dist = new TH1D("v24Dist", "v24Dist; M_{#mu#mu} (GeV/#it{c}^{2}); #it{v}_{2}^{#mu#mu}{4}", 100, minMass, maxMass);
    // TH1F *v24Dist = new TH1F("v24Dist", "", nBins, minMass, maxMass);

    TH2F *histMassV4 = new TH2F("histMassV4", "", nBins, 2, 4, 200, -1, 1);

    // Fill the histogram
    Long64_t nentries = chain->GetEntries();
    cout << "??? nentries " << nentries << endl;
    for (Long64_t i = 0; i < nentries; i++)
    {

        chain->GetEntry(i);

        if (pt1 >= 1.0 && pt2 >= 1.0)
        {
            if (fSign == 0 && pt >= minPt && pt <= maxPt && cent >= minCent && cent <= maxCent)
            {
                histMassV4->Fill(mass, V4);
            }
        }
    }

    TProfile *histProjV4 = (TProfile *)histMassV4->ProfileX("histProjV4");
    v24Dist->GetYaxis()->SetRangeUser(-0.04, 0.06);
    v24Dist->GetYaxis()->SetTitleOffset(1.6);

    for (Long64_t i = 0; i < nentries; i++)
    {

        float v4proj = histProjV4->GetBinContent(i + 1);
        float ev4proj = histProjV4->GetBinError(i + 1);
        v24Dist->SetBinContent(i + 1, v4proj);
        v24Dist->SetBinError(i + 1, ev4proj);
    }

    delete chain;

    return v24Dist; // Return the filled histogram
}

//////////////

TH1 *Getv22Histo(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent)
{
    // Create a TChain

    const char *treeName = "O2rtdimuonall";
    TChain *chain = new TChain(treeName);

    // Load data into the TChain
    cout << "GETINVHISTO ??? = " << filePattern << endl;
    // LoadData(chain, treeName, filePattern);
    //  filePattern = "/Users/valencia/Desktop/Flow/download/544124/AOD/002/AO2D.root";
    filePattern = "/Users/valencia/Desktop/AnalysisFlow/daiki/9nov/O2DQworkflows/O2DQworkflows/AO2D_test_ref_flow.root";

    string tree_str = treeName;
    string file_str = filePattern;
    Long64_t nHisto = getTChain(chain, file_str, tree_str);

    // Set up branches
    int nBins = 100;
    Float_t mass = 0;
    Float_t pt, pt1, pt2, V4 = 0;
    Float_t cent, corr2ref, corr2poi = 0;
    int fSign = -999;
    chain->SetBranchAddress("fMass", &mass);
    chain->SetBranchAddress("fPt", &pt);
    chain->SetBranchAddress("fCentFT0C", &cent);
    chain->SetBranchAddress("fPt1", &pt1);
    chain->SetBranchAddress("fPt2", &pt2);
    chain->SetBranchAddress("fV4", &V4);
    chain->SetBranchAddress("fSign", &fSign);

    chain->SetBranchAddress("fCORR2REF", &corr2ref);
    chain->SetBranchAddress("fCORR2POI", &corr2poi);

    chain->SetBranchAddress("fV4", &V4);

    TH1D *v22DistAll = new TH1D;

    v22DistAll = new TH1D("v22DistAll", "v22Dist; M_{#mu#mu} (GeV/#it{c}^{2}); #it{v}_{2}^{#mu#mu}{2}", 100, minMass, maxMass);
    v22DistAll->GetYaxis()->SetTitleOffset(1.6);

    TH2F *histMassCorr2Poi = new TH2F("histMassCorr2Poi", "", nBins, minMass, maxMass, 200, -1, 1);
    TH2F *histMassCorr2Ref = new TH2F("histMassCorr2Ref", "", nBins, minMass, maxMass, 200, -1, 1);

    // Global histo
    TH2F *histMassV22 = new TH2F("histMassV22", "", nBins, 2, 4, 200, -1, 1);

    Long64_t nentries = chain->GetEntries();
    cout << "??? nentries " << nentries << endl;
    for (Long64_t i = 0; i < nentries; i++)
    {

        chain->GetEntry(i);

        if (pt1 >= 1.0 && pt2 >= 1.0)
        {
            if (fSign == 0 && pt >= minPt && pt <= maxPt && cent >= minCent && cent <= maxCent)
            {

                histMassV22->Fill(mass, corr2poi / TMath::Sqrt(corr2ref));

                histMassCorr2Poi->Fill(mass, corr2poi);
                histMassCorr2Ref->Fill(mass, corr2ref);
            }
        }
    }

    TProfile *histProjV22 = (TProfile *)histMassV22->ProfileX("histProjV2");
    TProfile *histProjCorr2Poi = (TProfile *)histMassCorr2Poi->ProfileX("histProjCorr2Poi");
    TProfile *histProjCorr2Ref = (TProfile *)histMassCorr2Ref->ProfileX("histProjCorr2Ref");

    // cout << "N entries histProjV22 All " << histProjV22->GetEntries() << endl;
    Long64_t nEntriesProj = histProjV22->GetEntries();
    // ??? for (Long64_t i = 0; i < nentries; i++)
    for (Long64_t i = 0; i < nEntriesProj; i++)
    {

        float v22proj = histProjV22->GetBinContent(i + 1);
        float ev22proj = histProjV22->GetBinError(i + 1);

        float corr2Poi = histProjCorr2Poi->GetBinContent(i + 1);
        float corr2Ref = histProjCorr2Ref->GetBinContent(i + 1);

        float errCorr2Poi = histProjCorr2Poi->GetBinError(i + 1);
        float errCorr2Ref = histProjCorr2Ref->GetBinError(i + 1);

        // v22DistAll->SetBinContent(i + 1, VDN2Value(corr2Poi, corr2Ref));
        v22DistAll->SetBinContent(i + 1, v22proj);
        v22DistAll->SetBinError(i + 1, ev22proj);
        // v22DistAll->SetBinError(i + 1, VDN2Error(corr2Poi, errCorr2Poi, corr2Ref, errCorr2Ref));
    }

    delete chain;

    return v22DistAll; // Return the filled histogram
}

//////////////////

TH1 *GetvepHisto(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent)
{
    // Initialize TChain
    std::array<float, 3> dimuonMassRange = {minMass, maxMass, 20}; // Change Binning as needed....q
    std::array<float, 3> dimuonPtRange = {minPt, maxPt, 3};
    std::array<float, 2> dimuonCentRange{0, 50};

    //*************************************************************************************************************
    filePattern = "/Users/valencia/Desktop/AnalysisFlow/Flow/download/runs_chi/AO2D_merged_30april.root"; // OPTIONAL
    //*************************************************************************************************************

    TChain *fChain_POI = new TChain();
    LoadData(fChain_POI, "O2rtdileptonflow", filePattern);

    int NBinsMass = static_cast<int>(dimuonMassRange[2]);
    int NBinsPt = static_cast<int>(dimuonPtRange[2]);
    int NBinsMult = 10;
    double Cent[11] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
    double *MassBins = new double[NBinsMass + 1];
    double *PtBins = new double[NBinsPt + 1];

    CreateBins(MassBins, dimuonMassRange[0], dimuonMassRange[1], NBinsMass);
    CreateBins(PtBins, dimuonPtRange[0], dimuonPtRange[1], NBinsPt);

    TProfile2D *Cos2DeltaPhi = new TProfile2D(
        "Cos2DeltaPhi", "Cos2DeltaPhi", NBinsMass, MassBins, NBinsPt, PtBins);
    TProfile2D *R2EP = new TProfile2D("R2EP", "Resolution EP", NBinsMass,
                                      MassBins, NBinsPt, PtBins);

    TProfile2D *U2Q2 =
        new TProfile2D("U2Q2", "U2Q2", NBinsMass, MassBins, NBinsPt, PtBins);
    TProfile2D *R2SP_AB = new TProfile2D("R2SP_AB", "Resolution SP AB", NBinsMass,
                                         MassBins, NBinsPt, PtBins);
    TProfile2D *R2SP_BC = new TProfile2D("R2SP_BC", "Resolution SP BC", NBinsMass,
                                         MassBins, NBinsPt, PtBins);
    TProfile2D *R2SP_AC = new TProfile2D("R2SP_AC", "Resolution SP AC", NBinsMass,
                                         MassBins, NBinsPt, PtBins);

    float CentFT0POI;
    float fMass;
    float fPt, fPt1, fPt2;
    float fEta, fEta1, fEta2;
    int fMultDimuons = 0.;
    float fU2Q2, fU3Q3, fCos2DeltaPhi, fCos3DeltaPhi;
    float fR2EP, fR2SP_AB, fR2SP_AC, fR2SP_BC;

    fChain_POI->SetBranchAddress("fCentFT0C", &CentFT0POI);
    fChain_POI->SetBranchAddress("fMass", &fMass);
    fChain_POI->SetBranchAddress("fPt", &fPt);
    fChain_POI->SetBranchAddress("fEta", &fEta);

    fChain_POI->SetBranchAddress("fU2Q2", &fU2Q2);
    fChain_POI->SetBranchAddress("fU3Q3", &fU3Q3);
    fChain_POI->SetBranchAddress("fCos2DeltaPhi", &fCos2DeltaPhi);
    fChain_POI->SetBranchAddress("fCos3DeltaPhi", &fCos3DeltaPhi);
    fChain_POI->SetBranchAddress("fR2EP_BC", &fR2EP);
    fChain_POI->SetBranchAddress("fR2SP_AB", &fR2SP_AB);
    fChain_POI->SetBranchAddress("fR2SP_BC", &fR2SP_BC);
    fChain_POI->SetBranchAddress("fR2SP_AC", &fR2SP_AC);

    TH1F *hist_v2EP = new TH1F("v2EP", "v2EP", NBinsMass, MassBins);
    hist_v2EP->GetXaxis()->SetTitle("M_{#mu#mu} (GeV/#it{c}^{2})");
    hist_v2EP->GetYaxis()->SetTitle("#it{v}_{2}^{#mu#mu}{EP}");

    TH1F *hist_v2SP = new TH1F("v2SP", "v2SP", NBinsMass, MassBins);
    hist_v2SP->GetXaxis()->SetTitle("M_{#mu #mu} (GeV/c^{2})");
    hist_v2SP->GetYaxis()->SetTitle("v^{#mu #mu}_{2}{SP}");

    // Loop over all dimuons
    for (int i = 0; i < fChain_POI->GetEntries(); i++)
    {
        fChain_POI->GetEntry(i);
        // Dimuons general selection
        // {Pt range, msass bin, centrality range}
        if (!(fPt > dimuonPtRange[0] && fPt <= dimuonPtRange[1] &&
              fMass > dimuonMassRange[0] && fMass <= dimuonMassRange[1] &&
              fEta > -4. && fEta < -2.5 && CentFT0POI > dimuonCentRange[0] &&
              CentFT0POI <= dimuonCentRange[1]))
        {
            continue;
        }

        // Fill (mass, pt, centrality) bins for SP and EP
        if (!(isnan(fR2EP) || isinf(fR2EP)))
        {
            U2Q2->Fill(fMass, fPt, fU2Q2);
            R2SP_AB->Fill(fMass, fPt, fR2SP_AB);
            R2SP_BC->Fill(fMass, fPt, fR2SP_BC);
            R2SP_AC->Fill(fMass, fPt, fR2SP_AC);
            Cos2DeltaPhi->Fill(fMass, fPt, fCos2DeltaPhi);
            R2EP->Fill(fMass, fPt, fR2EP);
        }
    }

    // Mass-dependent flow for POI
    TProfile *U2Q2Mass = U2Q2->ProfileX("u2q2mass", 1, NBinsPt);
    TProfile *R2SPABMass = R2SP_AB->ProfileX("r2spmassAB", 1, NBinsPt);
    TProfile *R2SPBCMass = R2SP_BC->ProfileX("r2spmassBC", 1, NBinsPt);
    TProfile *R2SPACMass = R2SP_AC->ProfileX("r2spmassAC", 1, NBinsPt);
    TProfile *Cos2DeltaPhiMass = Cos2DeltaPhi->ProfileX("cos2deltaphimass", 1, NBinsPt);
    TProfile *R2EPMass = R2EP->ProfileX("r2epmass", 1, NBinsPt);

    for (int i = 0; i < NBinsMass; i++)
    {

        // Scalar-Product & Event-Plane method

        ///////////////////////////////////////////////
        float r2spab = R2SPABMass->GetBinContent(i + 1);
        float r2spab_e = R2SPABMass->GetBinError(i + 1);

        float r2spbc = R2SPBCMass->GetBinContent(i + 1);
        float r2spbc_e = R2SPBCMass->GetBinError(i + 1);

        float r2spac = R2SPACMass->GetBinContent(i + 1);
        float r2spac_e = R2SPACMass->GetBinError(i + 1);

        // Calculation of r2sp
        float r2sp = r2spab * r2spac / r2spbc;

        // Propagation of errors for r2sp
        float r2spab_part = (r2spac / r2spbc) * r2spab_e;
        float r2spac_part = (r2spab / r2spbc) * r2spac_e;
        float r2spbc_part = -1.0 * (r2spab * r2spac / (r2spbc * r2spbc)) * r2spbc_e;

        float r2spe = sqrt(r2spab_part * r2spab_part + r2spac_part * r2spac_part + r2spbc_part * r2spbc_part);

        ////////////////////////////////////

        float u2q2 = U2Q2Mass->GetBinContent(i + 1);
        float u2q2e = U2Q2Mass->GetBinError(i + 1);

        // float r2sp = R2SPMass->GetBinContent(i + 1);
        // float r2spe = R2SPMass->GetBinError(i + 1);

        float v2sp = u2q2 / pow(r2sp, 1. / 2);
        float v2spe = pow(u2q2e * u2q2e / r2sp +
                              0.25 * pow(r2sp, -3) * u2q2 * u2q2 * r2spe * r2spe,
                          1. / 2);

        float cos2deltaphi = Cos2DeltaPhiMass->GetBinContent(i + 1);
        float cos2deltaphie = Cos2DeltaPhiMass->GetBinError(i + 1);
        float r2ep = R2EPMass->GetBinContent(i + 1);
        float r2epe = R2EPMass->GetBinError(i + 1);
        float v2ep = cos2deltaphi / pow(r2ep, 1. / 2);
        float v2epe = pow(cos2deltaphie * cos2deltaphie / r2ep +
                              0.25 * pow(r2ep, -3) * cos2deltaphi * cos2deltaphi *
                                  r2epe * r2epe,
                          1. / 2);

        hist_v2SP->SetBinContent(i + 1, v2sp);
        hist_v2EP->SetBinContent(i + 1, v2ep);
        hist_v2SP->SetBinError(i + 1, isnan(v2spe) ? 0. : v2spe);
        hist_v2EP->SetBinError(i + 1, isnan(v2epe) ? 0. : v2epe);
    }

    return hist_v2SP;
}

void SaveFitResults(const std::string &title, const std::string &filename, Double_t v2jpsi, Double_t Errorv2jpsi)
{
    // Open an ofstream for output
    std::ofstream outFile(filename, std::ios::app); // Open in append mode directly in the constructor

    if (!outFile.is_open())
    {
        std::cerr << "Failed to open " << filename << " for writing." << std::endl;
        return;
    }

    // Write the title and the values to the file
    outFile << title << std::endl;
    outFile << "v2jpsi: " << v2jpsi << ", Errorv2jpsi: " << Errorv2jpsi << std::endl;

    // Close the file
    outFile.close();
}