// Author: Victor Valencia
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>
#include <TLatex.h>
#include <TDirectory.h>
// #include "SignalExtraction.C"
// #include "SignalExtraction.C"
#include "CORE/SignalExtractionNSAMPLES.C"

const char *fileName2 = "AO2D_test_ref_flow.root";
const char *fileLocation2 = "/Users/valencia/Desktop/AnalysisFlow/daiki/9nov/O2DQworkflows/O2DQworkflows";

// const char *fileName2 = "AO2D.root";
// const char *fileLocation2 = "/Users/valencia/Desktop/Flow/download/544124/AOD/002";

Float_t min_y = -4;
Float_t max_y = -2.5;
// const int numPtRanges = 1;
// Float_t min_pt[numPtRanges] = {0.0};
// Float_t max_pt[numPtRanges] = {12.0};
const int numPtRanges = 1;
Float_t min_pt[numPtRanges] = {0.0};
Float_t max_pt[numPtRanges] = {100.0};
Float_t min_cent = 0;
Float_t max_cent = 50;
Float_t min_mass = 2.;
Float_t max_mass = 4.;

Etails p_tails = kEMB;

const int numberFitRanges = 1;
Double_t min_fit[numberFitRanges] = {2.2};
Double_t max_fit[numberFitRanges] = {4.5};

Efunction configurations[1][4] = {
    // {kPol2Pol3, kCBExtended, kv2BackgroundPol3, kv2Pol3vsMassCB2Pol2Pol3},
    // {kVWG, kCBExtended, kv2BackgroundPol3, kv2Pol3vsMassCB2VWG},
    // {kVWGQuadratic, kCBExtended, kv2BackgroundPol3, kv2Pol3vsMassCB2QVWG},

    // {kPol2Pol3, kNA60, kv2BackgroundPol3, kv2Pol3vsMassNA60Pol2Pol3},
    // {kVWG, kNA60, kv2BackgroundPol3, kv2Pol3vsMassNA60VWG},
    // {kVWGQuadratic, kNA60, kv2BackgroundPol3, kv2Pol3vsMassNA60QVWG},

    // {kPol2Pol3, kCBExtended, kv2Background, kv2Pol2vsMassCB2Pol2Pol3},
    {kVWG, kCBExtended, kv2Background, kv2Pol2vsMassCB2VWG},
    // {kVWGQuadratic, kCBExtended, kv2Background, kv2Pol2vsMassCB2QVWG},

    // {kPol2Pol3, kNA60, kv2Background, kv2Pol2vsMassNA60Pol2Pol3},
    // {kVWG, kNA60, kv2Background, kv2Pol2vsMassNA60VWG},
    // {kVWGQuadratic, kNA60, kv2Background, kv2Pol2vsMassNA60QVWG},

    // {kPol2Pol3, kCBExtended, kv2BackgroundPolExp, kv2PolExpvsMassCB2Pol2Pol3},
    // {kVWG, kCBExtended, kv2BackgroundPolExp, kv2PolExpvsMassCB2VWG}
    // {kVWGQuadratic, kCBExtended, kv2BackgroundPolExp, kv2PolExpvsMassCB2QVWG},

    // {kPol2Pol3, kNA60, kv2BackgroundPolExp, kv2PolExpvsMassNA60Pol2Pol3},
    // {kVWG, kNA60, kv2BackgroundPolExp, kv2PolExpvsMassNA60VWG},
    // {kVWGQuadratic, kNA60, kv2BackgroundPolExp, kv2PolExpvsMassNA60QVWG}

};

void SysTotSignalExtractionPtNSAMPLES(const char *file = fileName2, Float_t minCent = min_cent, Float_t maxCent = max_cent, Efunction fBackgroundv2 = kv2BackgroundPol3, Efunction fSignalv2 = kv2Pol3vsMassCB2Pol2Pol3, Efunction fBackground = kPol2Pol3, Efunction fSignal = kCBExtended, Etails pTails = p_tails)
{

    std::ofstream outFile("output_extraction_txt/output_v2jpsi_sys.txt");

    for (int ptIndex = 0; ptIndex < numPtRanges; ptIndex++)
    {
        Float_t currentMinPt = min_pt[ptIndex];
        Float_t currentMaxPt = max_pt[ptIndex];
        TString ptRangeDirectory = Form("%s/%g-%gGeV", fileLocation2, currentMinPt, currentMaxPt);
        gSystem->mkdir(ptRangeDirectory, true); // Create a directory for each pT range

        for (int j = 0; j < numberFitRanges; j++)
        {
            for (int i = 0; i < 1; i++)
            {
                std::vector<Double_t> resultsTest;
                resultsTest.clear();
                resultsTest = SignalExtraction(Form("%s/%s", fileLocation2, file), min_y, max_y, currentMinPt, currentMaxPt, minCent, maxCent, min_mass, max_mass, configurations[i][0], configurations[i][1], configurations[i][2], configurations[i][3], pTails, min_fit[j], max_fit[j], kTRUE, kTRUE);

                if (resultsTest.size() == 0)
                    std::cout << "No signal extraction for configuration " << i + 1 << std::endl;
            }
        }
    }
}
