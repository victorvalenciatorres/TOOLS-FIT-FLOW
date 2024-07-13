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
#include <sstream>

#include "InitializeFunctions.C"
#include "SetRangeAndNameTest.C"
#include "getHistograms.C"

void SaveFitResults(const std::string &title, const std::string &filename, Double_t v2jpsi, Double_t Errorv2jpsi);

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

bool isTailsParameter(TString nameTails);

#define JPSI_MASS 3.096916
#define PSI2S_MASS 3.686109

//___________________________________________________________________________________________________________
// WITH Shuffling
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
            fitStatus1 = invMassDist->Fit("fitBG", "SN", "", 2.2, maxFit); // pas likelyhood
        else
            fitStatus1 = invMassDist->Fit("fitBG", "LSN", "", 2.2, maxFit); // Direct fit
        chi21 = BGFunction->GetChisquare() / BGFunction->GetNDF();
        covMatrixStatus1 = fitStatus1->CovMatrixStatus();

        secur1++;
        if (secur1 > 3)
            break;
        cout << "FitStatus = " << fitStatus1 << "     chi2/NDF = " << chi21 << "     cov matrix status = " << covMatrixStatus1 << endl;
    } while (fitStatus1 != 0 || chi21 > 2.5 || covMatrixStatus1 != 3);
    // while (fitStatus1 !=0 || covMatrixStatus1 != 3);

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
        // std::string massTitle = Form("Mass_ %s", nameTest.Data());

        std::string massTitle = Form("Mass_%s[Pt%.1f-%.1f]",
                                     nameTest.Data(),
                                     minPt, maxPt); // Adds minPt, maxPt to the title

        TCanvas *cInvMassDist = new TCanvas(Form("%s_%s", rangeNameFull.Data(), nameTest.Data()), Form("Range %s for test %s", rangeNameFull.Data(), nameTest.Data()));
        cInvMassDist->SetRightMargin(0.03);
        cInvMassDist->SetTopMargin(0.08);

        gStyle->SetOptFit(1111);
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
        gStyle->SetHistLineColor(1);
        // if(minCent<70)cInvMassDist->SetLogy();

        invMassDist->Draw("P");
        invMassDist->GetXaxis()->SetRangeUser(2, 5);
        fitFunction->SetLineColor(kBlue);
        fitFunction->SetRange(minFit, maxFit);
        fitFunction->Draw("SAME");

        BGFunction->SetLineColor(kBlue);
        BGFunction->SetLineStyle(7);
        BGFunction->SetRange(minFit, maxFit);
        BGFunction->Draw("SAME");
        JPsiFunction->SetLineColor(kRed);
        JPsiFunction->Draw("SAME");
        Psi2SFunction->SetLineColor(kGreen + 2);
        Psi2SFunction->Draw("SAME");

        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextAlign(12);
        text->SetTextFont(43);
        text->SetTextSize(19);
        text->DrawLatex(0.7, 0.85, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
        TLatex *text2 = new TLatex();
        text2->SetNDC();
        text2->SetTextAlign(12);
        text2->SetTextFont(43);
        text2->SetTextSize(17);
        text2->DrawLatex(0.7, 0.81, Form("Centrality %.1f-%.1f%%", minCent, maxCent));
        text2->DrawLatex(0.7, 0.78, "2.5 < #it{y} < 4");
        text2->DrawLatex(0.7, 0.74, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", minPt, maxPt));
        // text2->DrawLatex(0.7, 0.74, Form("#it{p}_{T} < %.1f GeV/#it{c}", maxPt));
        TLatex *text3 = new TLatex();
        text3->SetNDC();
        text3->SetTextAlign(12);
        text3->SetTextFont(43);
        text3->SetTextSize(15);
        text3->DrawLatex(0.7, 0.68, Form("N_{J/#psi} = %i #pm %i", TMath::Nint(N_JPsi), TMath::Nint(Err_JPsi)));
        for (int np = 1; np < 3; np++)
        {
            text3->DrawLatex(0.7, 0.68 - np * 0.04, Form("%s = %.3f #pm %.3f", fitFunction->GetParName(np + nb_bg), fitFunction->GetParameter(np + nb_bg), fitFunction->GetParError(np + nb_bg)));
        }
        text3->DrawLatex(0.7, 0.68 - 3 * 0.04, Form("#chi^{2}/NDF = %.3f", chi22));

        //////////////////////////////////////////////////////////////////////////////////////////////////////////// BEGINNING FITTING V2:

        int nSamples = 10; // SELECT NUMBER OF SAMPLES FOR SHUFFLING

        // Create the histogram
        TH1D **v24Dist;
        TH1D *v24DistAll;
        Double_t *errV24;
        TH1 *vepDist = new TH1F();

        // SELECT DESIRED HISTOGRAM TO FIT (Getv22Histo or Getv24Histo)

        std::pair<TH1D *, Double_t *> p = Getv22Histo(file, nSamples, minMass, maxMass, minPt, maxPt, minCent, maxCent);
        // std::pair<TH1D *, Double_t *> p = Getv22Histo(file, nSamples, minMass, maxMass, minPt, maxPt, minCent, maxCent);

        v24DistAll = p.first;
        errV24 = p.second;

        Double_t *params = fitFunction->GetParameters();

        TF1 *fitV2 = DistributionFunction(fBackGround, fSignal, fBackGroundv2, fSignalv2, minMass, maxMass);

        Double_t maxYv2 = -1.0;
        Double_t minYv2 = +1.0;

        invMassDist->SetTitle(Form("Mass_%s", nameTest.Data()));
        invMassDist->GetYaxis()->SetRangeUser(0, 50000);

        TCanvas *cv24Dist = new TCanvas();
        cv24Dist->SetRightMargin(0.03);
        cv24Dist->SetTopMargin(0.08);

        ////////////////////////////////////// FITTING V2

        fitV2->GetParameters();

        int securv2 = 0;
        Double_t chi2v2 = 3;

        Double_t v2jpsi_k[nSamples];
        Double_t v2jpsiCum = 0;
        Double_t v2jpsi = 0;
        Double_t Errorv2jpsi = 0;

        do
        {
            // cout << k << " sample " << v24Dist[k]->GetEntries() << endl;
            v24DistAll->Fit("fitV2", "R");
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

        v2jpsi = fitV2->GetParameter(GetNPar(fSignalv2) - 1);
        Errorv2jpsi = fitV2->GetParError(GetNPar(fSignalv2) - 1);

        Int_t nBinsAll = v24DistAll->GetNbinsX();

        // v24DistAll->Sumw2();   // SUMW2

        for (int i = 0; i < nBinsAll; ++i)
        {
            //cout << i << " binContent " << v24DistAll->GetBinContent(i) << endl;
            v24DistAll->SetBinError(i, errV24[i]);
        }

        TLatex *textv2 = new TLatex();
        textv2->SetNDC();
        textv2->SetTextAlign(12);
        textv2->SetTextFont(43);
        textv2->SetTextSize(17);
        textv2->DrawLatex(0.2, 0.85, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
        textv2->DrawLatex(0.2, 0.81, Form("Centrality %.1f-%.1f%%", minCent, maxCent));
        textv2->DrawLatex(0.2, 0.78, "2.5 < #it{y} < 4");
        // textv2->DrawLatex(0.7, 0.745, Form("#it{p}_{T} < %.1f GeV/#it{c}", maxPt));
        textv2->DrawLatex(0.2, 0.74, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", minPt, maxPt));

        // title
        textv2->DrawLatex(0.2, 0.7, Form("#it{v}_{2}^{J/#psi} = %.3f #pm %.3f", v2jpsi, Errorv2jpsi));
        textv2->DrawLatex(0.2, 0.65, Form("#chi^{2}/NDF = %.3f ", chi2v2));

        // Convert enum values to string names
        std::string bgName = GetFunctionName(fBackGroundv2);
        std::string signalName = GetFunctionName(fSignalv2);

        std::string title = Form("%sm[%.1f-%.1f][Pt%.1f-%.1f]",
                                 signalName.c_str(),
                                 minFit, maxFit,
                                 minPt, maxPt); // Adds minPt, maxPt to the title

        v24DistAll->SetTitle(title.c_str());

        fitV2->SetLineColor(kBlue);

        fitV2->Draw("same");

        delete[] errV24;

        ////////////// V2 BKG

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
                    // cout << "Parameters final" << i << " = " << paramsV4[paramIndex - nParamsBkgv2 + i] << endl;
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

            TLegend legend_linesV2(0.63, 0.77, 0.88, 0.87); // Define the legend position and size
            legend_linesV2.SetBorderSize(1);
            legend_linesV2.SetFillColor(0);
            legend_linesV2.SetTextSize(0.03);
            legend_linesV2.AddEntry(fitV2, "Total fit", "l");
            legend_linesV2.AddEntry(bkgV4, "Background", "l");
            legend_linesV2.Draw();

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

            cInvMassDist->SaveAs(("output_extraction_pdf/" + pdfPath + massTitle + "_NSAMPLES.pdf").c_str());
            cv24Dist->SaveAs(("output_extraction_pdf/" + pdfPath + title + "_NSAMPLES.pdf").c_str());
            std::string filename = "output_extraction_txt/output_v2jpsi_sys.txt";

            SaveFitResults(title, filename, v2jpsi, Errorv2jpsi);
        }
    }
    // analysis->Close();
    return results;
}

//___________________________________________________________________________________________________________
//__________________________________________________________________________________________________________

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