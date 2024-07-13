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

Long64_t getTChain(TChain *fChain, string FileName, string TreeName)
{
    cout << "[INFO] Extrating TTree " << TreeName.c_str() << endl;

    TFile *inputFile = TFile::Open(FileName.c_str());
    cout << "opening file " << FileName.c_str() << endl;
    TIter keyList(inputFile->GetListOfKeys());
    TKey *key;
    Long64_t nDF = 0;

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
        nDF += 1;
    }

    inputFile->Close();

    if (!fChain)
    {
        cout << "[ERROR] fChain was not created, some input files are missing" << endl;
        return -1;
    }
    return nDF;
};

//___________________________________________________________________________________________________________
//__________________________________________________________________________________________________________

// Function to get invariant mass histogram using a TChain
TH1 *GetInvMassHisto2(const char *filePattern, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent)
{

    const char *treeName = "O2rtdimuonall";
    TChain *chain = new TChain(treeName);

    cout << "GETINVHISTO ??? = " << filePattern << endl;
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
    TH1D *invMassDist = new TH1D("invMassDist", "Invariant Mass Distribution; m_{#mu#mu} (GeV/c^{2}); Counts per bin", 100, minMass, maxMass);

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
    filePattern = "/Users/valencia/Desktop/AnalysisFlow/Flow/download/runs_chi/AO2D_merged_30april.root"; // OPTIONAL
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
    TH1D *invMassDist = new TH1D("invMassDist", "Invariant Mass Distribution; m_{#mu#mu} (GeV/c^{2}); Counts per bin", 100, minMass, maxMass);
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
// std::pair<TH1D *, TH1D **> Getv24Histo(const char *filePattern, int nSamples, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent)
std::pair<TH1D *, Double_t *> Getv24Histo(const char *filePattern, int nSamples, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent)
{
    // Create a TChain

    const char *treeName = "O2rtdimuonall";
    TChain *chain = new TChain(treeName);

    // Load data into the TChain
    cout << "GETINVHISTO ??? = " << filePattern << endl;
    // LoadData(chain, treeName, filePattern);
    filePattern = "/Users/valencia/Desktop/AnalysisFlow/daiki/9nov/O2DQworkflows/O2DQworkflows/AO2D_test_ref_flow.root";

    string tree_str = treeName;
    string file_str = filePattern;
    Long64_t nHisto = getTChain(chain, file_str, tree_str);

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

    // // Create the histogram
    // TH1D *v24Dist = new TH1D("v24Dist", "v24Dist; m_{#mu#mu} (GeV/c^{2}); #it{v}_{2}^{#mu#mu}{4}", 100, minMass, maxMass)[nSamples];
    // // TH1F *v24Dist = new TH1F("v24Dist", "", nBins, minMass, maxMass);

    ///////////////////////////// added by me
    // TH1 **v24Dist = new *TH1 [nSamples];
    // for (int i = 0; i < nSamples; ++i)
    // {
    //     v24Dist[i] = new TH1D(Form("v24Dist_%d", i), "v24Dist; m_{#mu#mu} (GeV/c^{2}); #it{v}_{2}^{#mu#mu}", 100, minMass, maxMass);
    // }

    ////////////////////////////

    // (TH1*) v24Dist[nSamples]; // Declare an array of pointers to TH1D objects
    TH1D **v24Dist = new TH1D *[nSamples];
    TH1D *v24DistAll = new TH1D;
    ;
    // TH1 **v24Dist;
    v24DistAll = new TH1D("v24DistAll", "v24Dist; m_{#mu#mu} (GeV/c^{2}); #it{v}_{2}^{#mu#mu}", 100, minMass, maxMass);

    for (int i = 0; i < nSamples; ++i)
    {
        v24Dist[i] = new TH1D("v24Dist", "v24Dist; m_{#mu#mu} (GeV/c^{2}); #it{v}_{2}^{#mu#mu}", 100, minMass, maxMass);
    }

    // Fill the histogram

    // Global histo
    TH2F *histMassV4 = new TH2F("histMassV4", "", nBins, 2, 4, 200, -1, 1);
    // Sampling
    TH2F *histMassV4_k[nSamples];
    for (Int_t k = 0; k < nSamples; k++)
    {
        histMassV4_k[k] = new TH2F("histMassV4_k", "", nBins, 2, 4, 200, -1, 1);
    }
    // To get random between [0..nSamples-1]
    std::mt19937 mt(random_device{}());
    std::uniform_int_distribution<size_t> dist(0, nSamples - 1);
    //
    Long64_t nentries = chain->GetEntries();
    // cout << "??? nentries " << nentries << endl;
    for (Long64_t i = 0; i < nentries; i++)
    {

        chain->GetEntry(i);

        if (pt1 >= 1.0 && pt2 >= 1.0)
        {
            if (fSign == 0 && pt >= minPt && pt <= maxPt && cent >= minCent && cent <= maxCent)
            {
                histMassV4->Fill(mass, V4);
                // Draw random between [0..nSamples-1]
                Int_t k = dist(mt);
                histMassV4_k[k]->Fill(mass, V4);
            }
        }
    }

    TProfile *histProjV4 = (TProfile *)histMassV4->ProfileX("histProjV4");

    // cout << "N entries " << histProjV4->GetEntries() << endl;
    Long64_t nEntriesProj = histProjV4->GetEntries();

    for (Long64_t i = 0; i < nEntriesProj; i++)
    {

        float v4proj = histProjV4->GetBinContent(i + 1);
        float ev4proj = histProjV4->GetBinError(i + 1);
        v24DistAll->SetBinContent(i + 1, v4proj);
        v24DistAll->SetBinError(i + 1, ev4proj);
    }

    Double_t v24[nSamples][nBins];
    Double_t meanV24[nBins];
    Double_t *errV24 = new Double_t[nBins];
    for (Int_t i = 0; i < nBins; i++)
    {
        meanV24[i] = 0.0;
        errV24[i] = 0.0;
    }

    TCanvas *c1 = new TCanvas("c1", "TProfiles", 800, 600);
    c1->Divide(1, nSamples); // Divide the canvas into nSamples pads

    std::vector<TProfile *> profiles;
    TString xTitle = "m_{#mu#mu} (GeV/c^{2}";
    TString yTitle = "#it{v}_{2}^{#mu#mu}{4}";

    for (Int_t k = 0; k < nSamples; k++)
    {

        TString histName = TString::Format("histProjV4k_%d", k);
        TProfile *histProjV4k = (TProfile *)histMassV4_k[k]->ProfileX(histName);
        // TProfile *histProjV4k = (TProfile *)histMassV4_k[k]->ProfileX("histProjV4k");
        // cout << "entries histProjV24 " << histProjV4k->GetEntries() << endl;
        // cout << " X entries histProjV24 " << histProjV4k->GetNbinsX() << endl;
        Long64_t nBinContent = histProjV4k->GetNbinsX();
        for (Long64_t i = 0; i < nBinContent; i++)
        {
            v24[k][i] = histProjV4k->GetBinContent(i + 1);
            meanV24[i] += v24[k][i];
        }

        TString title = TString::Format("NSample %d, nEntries: %.1f, nEntriesTot: %.1f", k+1, histProjV4k->GetEntries(),histProjV4->GetEntries());

        histProjV4k->GetXaxis()->SetTitle(xTitle);
        histProjV4k->GetYaxis()->SetTitle(yTitle);
        histProjV4k->SetTitle(title);
        profiles.push_back(histProjV4k);
    }

    ////////////////////////////////////////////

    for (Int_t k = 0; k < nSamples; k++)
    {
        c1->cd(k + 1);
        profiles[k]->Draw();
    }

    c1->Update();

    ////////////////////////////////////////////////////
    for (Int_t i = 0; i < nBins; i++)
    {
        meanV24[i] = meanV24[i] / nSamples;
    }

    // Std Err
    for (Int_t k = 0; k < nSamples; k++)
    {
        for (Long64_t i = 0; i < nBins; i++)
        {
            Double_t x = v24[k][i] - meanV24[i];
            errV24[i] += (x * x);
        }
    }
    for (Int_t i = 0; i < nBins; i++)
    {
        errV24[i] = std::sqrt(errV24[i] / (nSamples -1));
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    delete chain;
    for (Int_t k = 0; k < nSamples; k++)
    {
        delete histMassV4_k[k];
    }

    // return v24Dist; // Return the filled histogram
    return std::make_pair(v24DistAll, errV24); // Return the filled histogram
    // ??? return static_cast<TH1 **>(v24Dist); // Return the filled histogram
}

////////////////////////////////////

std::pair<TH1D *, Double_t *> Getv22Histo(const char *filePattern, int nSamples, Float_t minMass, Float_t maxMass, Float_t minPt, Float_t maxPt, Float_t minCent, Float_t maxCent)
{
    // Create a TChain

    const char *treeName = "O2rtdimuonall";
    TChain *chain = new TChain(treeName);

    // Load data into the TChain
    cout << "GETINVHISTO ??? = " << filePattern << endl;
    // LoadData(chain, treeName, filePattern);
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

    // // Create the histogram
    // TH1D *v24Dist = new TH1D("v24Dist", "v24Dist; m_{#mu#mu} (GeV/c^{2}); #it{v}_{2}^{#mu#mu}{4}", 100, minMass, maxMass)[nSamples];
    // // TH1F *v24Dist = new TH1F("v24Dist", "", nBins, minMass, maxMass);

    ///////////////////////////// added by me
    // TH1 **v24Dist = new *TH1 [nSamples];
    // for (int i = 0; i < nSamples; ++i)
    // {
    //     v24Dist[i] = new TH1D(Form("v24Dist_%d", i), "v24Dist; m_{#mu#mu} (GeV/c^{2}); #it{v}_{2}^{#mu#mu}", 100, minMass, maxMass);
    // }

    ////////////////////////////

    // (TH1*) v24Dist[nSamples]; // Declare an array of pointers to TH1D objects
    TH1D **v22Dist = new TH1D *[nSamples];
    TH1D *v22DistAll = new TH1D;
    ;
    // TH1 **v24Dist;
    v22DistAll = new TH1D("v22DistAll", "v22Dist; m_{#mu#mu} (GeV/c^{2}); #it{v}_{2}^{#mu#mu}", 100, minMass, maxMass);

    for (int i = 0; i < nSamples; ++i)
    {
        v22Dist[i] = new TH1D("v22Dist", "v22Dist; m_{#mu#mu} (GeV/c^{2}); #it{v}_{2}^{#mu#mu}", 100, minMass, maxMass);
    }

    TH2F *histMassCorr2Poi = new TH2F("histMassCorr2Poi", "", nBins, minMass, maxMass, 200, -1, 1);
    TH2F *histMassCorr2Ref = new TH2F("histMassCorr2Ref", "", nBins, minMass, maxMass, 200, -1, 1);

    // Global histo
    TH2F *histMassV22 = new TH2F("histMassV22", "", nBins, 2, 4, 200, -1, 1);
    // Sampling
    TH2F *histMassV22_k[nSamples];
    for (Int_t k = 0; k < nSamples; k++)
    {
        histMassV22_k[k] = new TH2F("histMassV22_k", "", nBins, 2, 4, 200, -1, 1);
    }
    // To get random between [0..nSamples-1]
    std::mt19937 mt(random_device{}());
    std::uniform_int_distribution<size_t> dist(0, nSamples - 1);
    //
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

                // Draw random between [0..nSamples-1]
                Int_t k = dist(mt);
                histMassV22_k[k]->Fill(mass, corr2poi / TMath::Sqrt(corr2ref));


            }
        }
    }

    TProfile *histProjV22 = (TProfile *)histMassV22->ProfileX("histProjV2");
    TProfile *histProjCorr2Poi = (TProfile *)histMassCorr2Poi->ProfileX("histProjCorr2Poi");
    TProfile *histProjCorr2Ref = (TProfile *)histMassCorr2Ref->ProfileX("histProjCorr2Ref");

    cout << "N entries histProjV22 All " << histProjV22->GetEntries() << endl;
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

         v22DistAll->SetBinContent(i + 1, v22proj);
         v22DistAll->SetBinError(i + 1, ev22proj);

        // v22DistAll->SetBinContent(i + 1, VDN2Value(corr2Poi, corr2Ref));
        // v22DistAll->SetBinError(i + 1, VDN2Error(corr2Poi, errCorr2Poi, corr2Ref, errCorr2Ref));
    }



    Double_t v22[nSamples][nBins];

    Double_t meanV22[nBins];
    Double_t *errV22 = new Double_t[nBins];
    for (Int_t i = 0; i < nBins; i++)
    {
        meanV22[i] = 0.0;
        errV22[i] = 0.0;
    }

    TCanvas *c1 = new TCanvas("c1", "TProfiles", 800, 600);
    c1->Divide(1, nSamples); // Divide the canvas into nSamples pads

    std::vector<TProfile *> profiles;
    TString xTitle = "m_{#mu#mu} (GeV/c^{2}";
    TString yTitle = "#it{v}_{2}^{#mu#mu}{2}";

    Double_t meancorr2Ref;
    Double_t meancorr2Poi;

    for (Int_t k = 0; k < nSamples; k++)
    {

        TString histName = TString::Format("histProjV2k_%d", k);
        TProfile *histProjV22k = (TProfile *)histMassV22_k[k]->ProfileX(histName);


        // TProfile *histProjV4k = (TProfile *)histMassV4_k[k]->ProfileX("histProjV4k");
        // cout << "entries histProjV22 " << histProjV22k->GetEntries() << endl;
        // cout << " X entries histProjV22 " << histProjV22k->GetNbinsX() << endl;
        Long64_t nBinContent = histProjV22k->GetNbinsX();
        for (Long64_t i = 0; i < nBinContent; i++)
        {


            // corr2Ref[k][i] = histProjCorr2Ref_k->GetBinContent(i + 1);
            // corr2Poi[k][i] = histProjCorr2Poi_k->GetBinContent(i + 1);


            v22[k][i] = histProjV22k->GetBinContent(i + 1);
            // v22[k][i] =  VDN2Value(corr2poi[k][i], corr2ref[k][i]);
            meanV22[i] += v22[k][i];
        }


        histProjV22k->GetXaxis()->SetTitle(xTitle);
        histProjV22k->GetYaxis()->SetTitle(yTitle);
        profiles.push_back(histProjV22k);
    }

    ////////////////////////////////////////////

    for (Int_t k = 0; k < nSamples; k++)
    {
        c1->cd(k + 1);
        profiles[k]->Draw();
    }

    c1->Update();

    ////////////////////////////////////////////////////
    for (Int_t i = 0; i < nBins; i++)
    {
        meanV22[i] = meanV22[i] / nSamples;
    }

    // Std Err
    for (Int_t k = 0; k < nSamples; k++)
    {
        for (Long64_t i = 0; i < nBins; i++)
        {
            Double_t x = v22[k][i] - meanV22[i];
            errV22[i] += (x * x);
        }
    }
    for (Int_t i = 0; i < nBins; i++)
    {
        errV22[i] = std::sqrt(errV22[i] / (nSamples-1));
    }


    delete chain;
    for (Int_t k = 0; k < nSamples; k++)
    {
        delete histMassV22_k[k];
    }

    // return v24Dist; // Return the filled histogram
    return std::make_pair(v22DistAll, errV22); // Return the filled histogram
    // ??? return static_cast<TH1 **>(v24Dist); // Return the filled histogram
}



/////////////////////////////////////


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
    hist_v2EP->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
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

