// Author: Victor Valencia
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include "TVector.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include <cctype>
#include <map>
#include <unordered_map>

#include "FitFunctions.C"

double par_v2Pol2vsMassCB2VWG[15];
double par_v2Pol2vsMassCB2QVWG[16];
double par_v2Pol2vsMassCB2Pol2Pol3[17];

double par_v2Pol2vsMassNA60QVWG[20];
double par_v2Pol2vsMassNA60VWG[19];
double par_v2Pol2vsMassNA60Pol2Pol3[21];

double par_v2Pol3vsMassCB2VWG[16];
double par_v2Pol3vsMassCB2QVWG[17];
double par_v2Pol3vsMassCB2Pol2Pol3[18];

double par_v2Pol3vsMassNA60QVWG[21];

double par_v2Pol3vsMassNA60VWG[20];

double par_v2Pol3vsMassNA60Pol2Pol3[22];

double par_v2PolExpvsMassCB2VWG[17];

double par_v2PolExpvsMassCB2Pol2Pol3[19];

double par_v2PolExpvsMassNA60QVWG[22];

double par_v2PolExpvsMassNA60VWG[21];

double par_v2PolExpvsMassNA60Pol2Pol3[23];

double par_v2PolExpvsMassCB2QVWG[18];

enum class V2MassModel
{
    PolExpvsMassCB2VWG = 0,
    PolExpvsMassCB2Pol2Pol3 = 1,
    PolExpvsMassNA60QVWG = 2,
    PolExpvsMassNA60VWG = 3,
    PolExpvsMassNA60Pol2Pol3 = 4,
    PolExpvsMassCB2QVWG = 5,

    Pol2vsMassCB2VWG = 6,
    Pol2vsMassCB2Pol2Pol3 = 7,
    Pol2vsMassNA60QVWG = 8,
    Pol2vsMassNA60VWG = 9,
    Pol2vsMassNA60Pol2Pol3 = 10,
    Pol2vsMassCB2QVWG = 11,

    Pol3vsMassCB2VWG = 12,
    Pol3vsMassCB2Pol2Pol3 = 13,
    Pol3vsMassNA60QVWG = 14,
    Pol3vsMassNA60VWG = 15,
    Pol3vsMassNA60Pol2Pol3 = 16,
    Pol3vsMassCB2QVWG = 17
};

map<string, V2MassModel> V2MassModelDictionary = {
    {"PolExpvsMassCB2VWG", V2MassModel::PolExpvsMassCB2VWG},
    {"PolExpvsMassCB2Pol2Pol3", V2MassModel::PolExpvsMassCB2Pol2Pol3},
    {"PolExpvsMassNA60QVWG", V2MassModel::PolExpvsMassNA60QVWG},
    {"PolExpvsMassNA60VWG", V2MassModel::PolExpvsMassNA60VWG},
    {"PolExpvsMassNA60Pol2Pol3", V2MassModel::PolExpvsMassNA60Pol2Pol3},
    {"PolExpvsMassCB2QVWG", V2MassModel::PolExpvsMassCB2QVWG},

    {"Pol2vsMassCB2VWG", V2MassModel::Pol2vsMassCB2VWG},
    {"Pol2vsMassCB2Pol2Pol3", V2MassModel::Pol2vsMassCB2Pol2Pol3},
    {"Pol2vsMassNA60QVWG", V2MassModel::Pol2vsMassNA60QVWG},
    {"Pol2vsMassNA60VWG", V2MassModel::Pol2vsMassNA60VWG},
    {"Pol2vsMassNA60Pol2Pol3", V2MassModel::Pol2vsMassNA60Pol2Pol3},
    {"Pol2vsMassCB2QVWG", V2MassModel::Pol2vsMassCB2QVWG},

    {"Pol3vsMassCB2VWG", V2MassModel::Pol3vsMassCB2VWG},
    {"Pol3vsMassCB2Pol2Pol3", V2MassModel::Pol3vsMassCB2Pol2Pol3},
    {"Pol3vsMassNA60QVWG", V2MassModel::Pol3vsMassNA60QVWG},
    {"Pol3vsMassNA60VWG", V2MassModel::Pol3vsMassNA60VWG},
    {"Pol3vsMassNA60Pol2Pol3", V2MassModel::Pol3vsMassNA60Pol2Pol3},
    {"Pol3vsMassCB2QVWG", V2MassModel::Pol3vsMassCB2QVWG}

};

void addV2MassModel(V2MassModel model, vector<double> params)
{

    switch (model)
    {
    case (V2MassModel::PolExpvsMassCB2VWG):
        if (params.size() != 17)
        {
            cout << "Not right nb of parameters ! 17vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 17; i++)
        {
            par_v2PolExpvsMassCB2VWG[i] = params[i];
            cout << "Setting PolExpvsMassCB2VWG param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::PolExpvsMassCB2Pol2Pol3):
        if (params.size() != 19)
        {
            cout << "Not right nb of parameters ! 19vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 19; i++)
        {
            par_v2PolExpvsMassCB2Pol2Pol3[i] = params[i];
            cout << "Setting PolExpvsMassCB2Pol2Pol3 param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::PolExpvsMassNA60QVWG):
        if (params.size() != 22)
        {
            cout << "Not right nb of parameters ! 22vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 22; i++)
        {
            par_v2PolExpvsMassNA60QVWG[i] = params[i];
            cout << "Setting PolExpvsMassNA60QVWG param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::PolExpvsMassNA60VWG):
        if (params.size() != 21)
        {
            cout << "Not right nb of parameters ! 21vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 21; i++)
        {
            par_v2PolExpvsMassNA60VWG[i] = params[i];
            cout << "Setting PolExpvsMassNA60VWG param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::PolExpvsMassNA60Pol2Pol3):
        if (params.size() != 23)
        {
            cout << "Not right nb of parameters ! 23vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 23; i++)
        {
            par_v2PolExpvsMassNA60Pol2Pol3[i] = params[i];
            cout << "Setting PolExpvsMassNA60Pol2Pol3 param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::PolExpvsMassCB2QVWG):
        if (params.size() != 18)
        {
            cout << "Not right nb of parameters ! 18vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 18; i++)
        {
            par_v2PolExpvsMassCB2QVWG[i] = params[i];
            cout << "Setting PolExpvsMassCB2QVWG param n°" << i << "=" << params[i] << endl;
        }
        break;

    case (V2MassModel::Pol3vsMassCB2VWG):
        if (params.size() != 16)
        {
            cout << "Not right nb of parameters ! 17vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 16; i++)
        {
            par_v2Pol3vsMassCB2VWG[i] = params[i];
            cout << "Setting Pol3vsMassCB2VWG param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol3vsMassCB2Pol2Pol3):
        if (params.size() != 18)
        {
            cout << "Not right nb of parameters ! 19vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 18; i++)
        {
            par_v2Pol3vsMassCB2Pol2Pol3[i] = params[i];
            cout << "Setting Pol3vsMassCB2Pol2Pol3 param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol3vsMassNA60QVWG):
        if (params.size() != 21)
        {
            cout << "Not right nb of parameters ! 22vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 21; i++)
        {
            par_v2Pol3vsMassNA60QVWG[i] = params[i];
            cout << "Setting Pol3vsMassNA60QVWG param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol3vsMassNA60VWG):
        if (params.size() != 20)
        {
            cout << "Not right nb of parameters ! 21vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 20; i++)
        {
            par_v2Pol3vsMassNA60VWG[i] = params[i];
            cout << "Setting Pol3vsMassNA60VWG param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol3vsMassNA60Pol2Pol3):
        if (params.size() != 22)
        {
            cout << "Not right nb of parameters ! 23vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 22; i++)
        {
            par_v2Pol3vsMassNA60Pol2Pol3[i] = params[i];
            cout << "Setting Pol3vsMassNA60Pol2Pol3 param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol3vsMassCB2QVWG):
        if (params.size() != 17)
        {
            cout << "Not right nb of parameters ! 18vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 17; i++)
        {
            par_v2Pol3vsMassCB2QVWG[i] = params[i];
            cout << "Setting Pol3vsMassCB2QVWG param n°" << i << "=" << params[i] << endl;
        }
        break;

    case (V2MassModel::Pol2vsMassCB2VWG):
        if (params.size() != 15)
        {
            cout << "Not right nb of parameters ! 17vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 15; i++)
        {
            par_v2Pol2vsMassCB2VWG[i] = params[i];
            cout << "Setting Pol2sMassCB2VWG param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol2vsMassCB2Pol2Pol3):
        if (params.size() != 17)
        {
            cout << "Not right nb of parameters ! 19vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 17; i++)
        {
            par_v2Pol2vsMassCB2Pol2Pol3[i] = params[i];
            cout << "Setting Pol2sMassCB2Pol2Pol3 param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol2vsMassNA60QVWG):
        if (params.size() != 20)
        {
            cout << "Not right nb of parameters ! 22vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 20; i++)
        {
            par_v2Pol2vsMassNA60QVWG[i] = params[i];
            cout << "Setting Pol2sMassNA60QVWG param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol2vsMassNA60VWG):
        if (params.size() != 19)
        {
            cout << "Not right nb of parameters ! 21vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 19; i++)
        {
            par_v2Pol2vsMassNA60VWG[i] = params[i];
            cout << "Setting Pol2sMassNA60VWG param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol2vsMassNA60Pol2Pol3):
        if (params.size() != 21)
        {
            cout << "Not right nb of parameters ! 23vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 21; i++)
        {
            par_v2Pol2vsMassNA60Pol2Pol3[i] = params[i];
            cout << "Setting Pol2sMassNA60Pol2Pol3 param n°" << i << "=" << params[i] << endl;
        }
        break;
    case (V2MassModel::Pol2vsMassCB2QVWG):
        if (params.size() != 16)
        {
            cout << "Not right nb of parameters ! 18vs" << params.size() << endl;
            return;
        }
        for (int i = 0; i < 16; i++)
        {
            par_v2Pol2vsMassCB2QVWG[i] = params[i];
            cout << "Setting Pol2sMassCB2QVWG param n°" << i << "=" << params[i] << endl;
        }
        break;

    default:
        cout << "Not found" << endl;
        break;
    }
}

Double_t par_VWG[4] = {1., 1.235, 1.600, 2.928};
TString name_VWG[4] = {"Norm_{VWG}", "#mu_{VWG}", "#alpha_{VWG}", "#beta_{VWG}"};

Double_t par_VWGQuadratic[5] = {1., 2.3, 0.5, 0.3, -0.79}; // low pT : 0-10, 10-30, 30-50, 0-90

TString name_VWGQuadratic[5] = {"Norm_{VWG}", "#mu_{VWG}", "#alpha_{VWG}", "#beta_{VWG}", "#gamma_{VWG}"};

Double_t par_Pol2Pol3[6] = {1., -0.35, 0.035, 6.5, -5.1, 1.04}; // low pT : 0-10, 10-30, 30-50, 0-90

TString name_Pol2Pol3[6] = {"Norm", "a1", "a2", "b1", "b2", "b3"};
TString name_v2Background[3] = {"a1", "a2", "a3"};
TString name_v2BackgroundPol3[4] = {"a1", "a2", "a3", "a4"};
TString name_v2BackgroundPolExp[5] = {"a1", "a2", "a3", "a4", "a5"};
Double_t par_v2bkg[3] = {1., 1., 1.};
Double_t par_v2bkgPol3[4] = {1., 1., 1., 1.};
Double_t par_v2bkgPolExp[5] = {1., 1., 1., 1., 1.};

Double_t par_DoubleExp[4] = {1265, -1.4, -0.04, -0.04};
TString name_DoubleExp[4] = {"Norm1_{DoubleExp}", "#alpha1__{DoubleExp}", "Norm2_{DoubleExp}", "#alpha2__{DoubleExp}"};

Double_t par_Exp[2] = {1, 1};
TString name_Exp[2] = {"Norm1_{DoubleExp}", "#alpha1__{DoubleExp}"};

// Initial parameters to fit signal
Double_t par_jpsi[3] = {1., 3.1, 0.065};
Double_t par_psi[3] = {1.7, par_jpsi[1] + 3.686109 - 3.096916, par_jpsi[2] * 1.05};

TString name_CBjpsi[7] = {"Norm_{J/#psi}", "#mu_{J/#psi}", "#sigma_{J/#psi}", "#alpha^{L}_{J/#psi}", "n^{L}_{J/#psi}", "#alpha^{R}_{J/#psi}", "n^{R}_{J/#psi}"};
TString name_CBpsi[3] = {"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}"};
TString name_CBpsi_full[7] = {"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}", "#alpha^{L}_{#psi'}", "n^{L}_{#psi'}", "#alpha^{R}_{#psi'}", "n^{R}_{#psi'}"};

TString name_NA60jpsi[11] = {"Norm_{J/#psi}", "#mu_{J/#psi}", "#sigma_{J/#psi}", "#alpha^{L}_{J/#psi}", "p1^{L}_{J/#psi}", "p2^{L}_{J/#psi}", "p3^{L}_{J/#psi}", "#alpha^{R}_{J/#psi}", "p1^{R}_{J/#psi}", "p2^{R}_{J/#psi}", "p3^{R}_{J/#psi}"};
TString name_NA60psi[3] = {"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}"};
TString name_NA60psi_full[11] = {"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}", "#alpha^{L}_{#psi'}", "p1^{L}_{#psi'}", "p2^{L}_{#psi'}", "p3^{L}_{#psi'}", "#alpha^{R}_{#psi'}", "p1^{R}_{#psi'}", "p2^{R}_{#psi'}", "p3^{R}_{#psi'}"};

TString name_v2Pol2vsMassCB2VWG[15] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"};
TString name_v2Pol2vsMassCB2Pol2Pol3[17] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"};
TString name_v2Pol2vsMassCB2QVWG[16] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
TString name_v2Pol2vsMassNA60VWG[19] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"};
TString name_v2Pol2vsMassNA60QVWG[20] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"};
TString name_v2Pol2vsMassNA60Pol2Pol3[21] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"};

TString name_v2Pol3vsMassCB2VWG[16] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
TString name_v2Pol3vsMassCB2Pol2Pol3[18] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"};
TString name_v2Pol3vsMassCB2QVWG[17] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"};
TString name_v2Pol3vsMassNA60VWG[20] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"};
TString name_v2Pol3vsMassNA60QVWG[21] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"};
TString name_v2Pol3vsMassNA60Pol2Pol3[22] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"};

TString name_v2PolExpvsMassCB2VWG[17] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"};
TString name_v2PolExpvsMassCB2Pol2Pol3[19] = {
    "0",
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
};
TString name_v2PolExpvsMassCB2QVWG[18] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17"};
TString name_v2PolExpvsMassNA60VWG[21] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"};
TString name_v2PolExpvsMassNA60QVWG[22] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"};
TString name_v2PolExpvsMassNA60Pol2Pol3[23] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};

// READ   output.txt file
std::vector<double> GetTailsParameter(TString nameTails)
{
    std::vector<double> param_tails;

    // Construct the file path.
    std::string filePath = "input_extraction_txt/params_input.txt";

    // Open the input file.
    std::ifstream inFile(filePath);
    if (!inFile)
    {
        std::cerr << "Unable to open file: " << filePath << std::endl;
        return param_tails; // Return an empty vector if file cannot be opened.
    }

    std::string line;
    bool collectData = false; // Flag to start collecting data

    // Read lines from the file.
    while (std::getline(inFile, line))
    {
        // Check if the line matches the section we are interested in.
        if (line == nameTails)
        {
            collectData = true;
            continue; // Skip the header line
        }

        if (collectData && !line.empty() && std::isalpha(line[0]))
        {
            // Stop collecting data if another section starts
            collectData = false;
            continue;
        }

        // Convert line to a double and add to the vector if we're in the right section.
        if (collectData)
        {
            try
            {
                double value = std::stod(line);
                param_tails.push_back(value);
            }
            catch (const std::invalid_argument &e)
            {
                std::cerr << "Invalid number found in file: " << line << std::endl;
                continue;
            }
        }
    }

    inFile.close(); // Close the input file.

    // Output the parameters for debugging or verification.
    for (double param : param_tails)
    {
        std::cout << "Parameter: " << param << std::endl;
    }

    return param_tails;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TF1 *BackGroundFunction(Efunction fName, Double_t xmin, Double_t xmax)
{
    Double_t *par = nullptr;
    TString *name_par = nullptr;
    TF1 *BGFunction = nullptr;
    Int_t nbPar = GetNPar(fName);

    switch (fName)
    {
    case kVWG:
        par = par_VWG;
        name_par = name_VWG;
        BGFunction = new TF1("fitBG", VWG, xmin, xmax, GetNPar(kVWG));
        break;

    case kVWGQuadratic:
        par = par_VWGQuadratic;
        name_par = name_VWGQuadratic;
        BGFunction = new TF1("fitBG", VWGQuadratic, xmin, xmax, GetNPar(kVWGQuadratic));
        break;

    case kPol2Pol3:
        par = par_Pol2Pol3;
        name_par = name_Pol2Pol3;
        BGFunction = new TF1("fitBG", Pol2Pol3, xmin, xmax, GetNPar(kPol2Pol3));
        break;

    case kDoubleExp:
        par = par_DoubleExp;
        name_par = name_DoubleExp;
        BGFunction = new TF1("fitBG", DoubleExp, xmin, xmax, GetNPar(kDoubleExp));
        break;

    case kExp:
        par = par_Exp;
        name_par = name_Exp;
        BGFunction = new TF1("fitBG", Exp, xmin, xmax, GetNPar(kExp));
        break;

    case kv2Background:
        par = par_v2bkg;
        name_par = name_v2Background;
        BGFunction = new TF1("fitBG", v2Background, xmin, xmax, GetNPar(kv2Background));

        break;

    case kv2BackgroundPol3:
        par = par_v2bkgPol3;
        name_par = name_v2BackgroundPol3;
        BGFunction = new TF1("fitBG", v2BackgroundPol3, xmin, xmax, GetNPar(kv2BackgroundPol3));
        break;

    case kv2BackgroundPolExp:
        par = par_v2bkgPolExp;
        name_par = name_v2BackgroundPolExp;
        BGFunction = new TF1("fitBG", v2BackgroundPolExp, xmin, xmax, GetNPar(kv2BackgroundPolExp));
        break;
    }

    for (int i = 0; i < nbPar; i++)
    {
        BGFunction->SetParameter(i, par[i]);
        BGFunction->SetParName(i, name_par[i]);
    }
    return BGFunction;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TF1 *SignalFunction(Efunction fSig, TString nTails, Epart fPart, Double_t xmin, Double_t xmax)
{
    Double_t *par_signal = nullptr;
    TString *name_par_signal = nullptr;
    std::vector<Double_t> par_tails = GetTailsParameter(nTails);

    Int_t nb_par_signal = GetNPar(fSig);

    TF1 *fSignal = nullptr;

    switch (fSig)
    {
    case kCBExtended:
        if (fPart == kJPsi)
        {
            par_signal = par_jpsi;
            name_par_signal = name_CBjpsi;
        }
        else if (fPart == kPsi2S)
        {
            par_signal = par_psi;
            name_par_signal = name_CBpsi_full;
        }

        fSignal = new TF1("fitSignal", CrystalBallExtended, xmin, xmax, GetNPar(kCBExtended));
        break;

    case kNA60:
        if (fPart == kJPsi)
        {
            par_signal = par_jpsi;
            name_par_signal = name_NA60jpsi;
        }
        else if (fPart == kPsi2S)
        {
            par_signal = par_psi;
            name_par_signal = name_NA60psi_full;
        }

        fSignal = new TF1("fitSignal", NA60, xmin, xmax, GetNPar(kNA60));

        break;
    }

    for (int i = 0; i < nb_par_signal; i++)
    {
        // Norm, mean and sigma JPsi
        if (i < 3)
            fSignal->SetParameter(i, par_signal[i]);
        // Tails Jpsi
        else
            fSignal->FixParameter(i, par_tails[i - 3]);

        fSignal->SetParName(i, name_par_signal[i]);
    }

    return fSignal;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TF1 *DistributionFunction(Efunction fBG, Efunction fSig, TString nTails, Double_t xmin, Double_t xmax)
{
    Double_t *par_bg = nullptr;
    Double_t *par_signal1 = nullptr;
    Double_t *par_signal2 = nullptr;
    std::vector<Double_t> par_tails = GetTailsParameter(nTails);

    TString *name_par_bg = nullptr;
    TString *name_par_signal1 = nullptr;
    TString *name_par_signal2 = nullptr;

    Int_t nb_par_bg = GetNPar(fBG);
    Int_t nb_par_signal = GetNPar(fSig);
    Int_t nb_par = nb_par_bg + nb_par_signal + 1;

    Int_t fchoice = fBG + fSig;

    TF1 *fDistrib = nullptr;

    switch (fBG)
    {
    case kVWGQuadratic:
        par_bg = par_VWGQuadratic;
        name_par_bg = name_VWGQuadratic;
        break;

    case kPol2Pol3:
        par_bg = par_Pol2Pol3;
        name_par_bg = name_Pol2Pol3;
        break;

    case kDoubleExp:
        par_bg = par_DoubleExp;
        name_par_bg = name_DoubleExp;
        break;

    case kExp:
        par_bg = par_Exp;
        name_par_bg = name_Exp;
        break;
    case kVWG:
        par_bg = par_VWG;
        name_par_bg = name_VWG;
        break;
    }

    switch (fSig)
    {
    case kCBExtended:
        par_signal1 = par_jpsi;
        par_signal2 = par_psi;
        name_par_signal1 = name_CBjpsi;
        name_par_signal2 = name_CBpsi;
        break;

    case kNA60:
        par_signal1 = par_jpsi;
        par_signal2 = par_psi;
        name_par_signal1 = name_NA60jpsi;
        name_par_signal2 = name_NA60psi;
        break;
    }

    switch (fchoice)
    {
    case 220:
        fDistrib = new TF1("fitDistrib", VWGquad_DoubleCBext, xmin, xmax, GetNPar(kVWGQuadratic) + GetNPar(kCBExtended) + 1);
        break;

    case 230:
        fDistrib = new TF1("fitDistrib", VWGquad_DoubleNA60, xmin, xmax, GetNPar(kVWGQuadratic) + GetNPar(kNA60) + 1);
        break;

    case 320:
        fDistrib = new TF1("fitDistrib", Pol2Pol3_DoubleCBext, xmin, xmax, GetNPar(kPol2Pol3) + GetNPar(kCBExtended) + 1);
        break;

    case 330:
        fDistrib = new TF1("fitDistrib", Pol2Pol3_NA60, xmin, xmax, GetNPar(kPol2Pol3) + GetNPar(kNA60) + 1);
        break;

    case 420:
        fDistrib = new TF1("fitDistrib", DoubleExp_DoubleCBext, xmin, xmax, GetNPar(kDoubleExp) + GetNPar(kCBExtended) + 1);
        break;

    case 430:
        fDistrib = new TF1("fitDistrib", DoubleExp_DoubleNA60, xmin, xmax, GetNPar(kDoubleExp) + GetNPar(kNA60) + 1);
        break;

    case 520:
        fDistrib = new TF1("fitDistrib", Exp_DoubleCBext, xmin, xmax, GetNPar(kExp) + GetNPar(kCBExtended) + 1);
        break;

    case 530:
        fDistrib = new TF1("fitDistrib", Exp_DoubleNA60, xmin, xmax, GetNPar(kExp) + GetNPar(kNA60) + 1);
        break;

    case 120:
        fDistrib = new TF1("fitDistrib", VWG_CBext, xmin, xmax, GetNPar(kVWG) + GetNPar(kCBExtended) + 1);
        break;
    case 130:
        fDistrib = new TF1("fitDistrib", VWG_NA60, xmin, xmax, GetNPar(kVWG) + GetNPar(kNA60) + 1);
        break;
    }

    // Background
    for (int i = 0; i < nb_par_bg; i++)
    {
        fDistrib->SetParameter(i, par_bg[i]);
        fDistrib->SetParName(i, name_par_bg[i]);
    }
    // Signal
    for (int j = 0; j < nb_par_signal; j++)
    {
        // Norm, mean and sigma JPsi
        if (j < 3)
        {
            fDistrib->SetParameter(j + nb_par_bg, par_signal1[j]);
            if (j == 1)
            {
                fDistrib->SetParLimits(j + nb_par_bg, 2.9, 3.3);
            }
            if (j == 2)
            {
                fDistrib->SetParLimits(j + nb_par_bg, 0.04, 0.1);
            }
        }
        // Tails JPsi
        else
            fDistrib->FixParameter(j + nb_par_bg, par_tails[j - 3]);
        fDistrib->SetParName(j + nb_par_bg, name_par_signal1[j]);
    }
    // psi2s
    fDistrib->SetParameter(nb_par_bg + nb_par_signal, par_signal2[0]);
    fDistrib->SetParLimits(nb_par_bg + nb_par_signal, 0, 200000000);
    fDistrib->SetParName(nb_par_bg + nb_par_signal, name_par_signal2[0]);
    return fDistrib;
}

TF1 *DistributionFunction(Efunction fBG, Efunction fSig, Efunction fBGv2, Efunction fSigv2, Double_t xmin, Double_t xmax)
{
    // Calculate total number of parameters for the original and v2 function types
    Int_t nb_par_bg = GetNPar(fBG);
    Int_t nb_par_signal = GetNPar(fSig);
    Int_t nb_par_bgv2 = GetNPar(fBGv2);
    Int_t nb_par_signalv2 = GetNPar(fSigv2);

    // Allocate parameter and name arrays based on function type
    Double_t *par_bg = nullptr, *par_bgv2 = nullptr;
    Double_t *par_signal1 = nullptr, *par_signal1v2 = nullptr;
    TString *name_par_bg = nullptr, *name_par_bgv2 = nullptr;
    TString *name_par_signal1 = nullptr, *name_par_signal1v2 = nullptr;

    // Assign pointers based on function types
    switch (fBG)
    {
    case kVWGQuadratic:
        par_bg = par_VWGQuadratic;
        name_par_bg = name_VWGQuadratic;
        break;
    case kPol2Pol3:
        par_bg = par_Pol2Pol3;
        name_par_bg = name_Pol2Pol3;
        break;
    case kDoubleExp:
        par_bg = par_DoubleExp;
        name_par_bg = name_DoubleExp;
        break;
    case kExp:
        par_bg = par_Exp;
        name_par_bg = name_Exp;
        break;
    case kVWG:
        par_bg = par_VWG;
        name_par_bg = name_VWG;
        break;
    }
    switch (fBGv2)
    {
    case kVWG:
        par_bgv2 = par_VWG;
        name_par_bgv2 = name_VWG;
        break;
    case kVWGQuadratic:
        par_bg = par_VWGQuadratic;
        name_par_bg = name_VWGQuadratic;
        break;
    case kPol2Pol3:
        par_bg = par_Pol2Pol3;
        name_par_bg = name_Pol2Pol3;
        break;
    case kv2Background:
        par_bg = par_v2bkg;
        name_par_bg = name_v2Background;

        break;

    case kv2BackgroundPol3:
        par_bg = par_v2bkgPol3;
        name_par_bg = name_v2BackgroundPol3;

        break;
    case kv2BackgroundPolExp:
        par_bg = par_v2bkgPolExp;
        name_par_bg = name_v2BackgroundPolExp;

        break;
    }

    switch (fSig)
    {
    case kCBExtended:
    case kNA60:
        par_signal1 = par_jpsi;
        name_par_signal1 = name_CBjpsi;
        break;
    }

    // Create the distribution function with a function type and total number of parameters
    TF1 *fDistrib = new TF1();

    switch (fSigv2)
    {
    case kv2Pol2vsMassCB2VWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol2vsMassCB2VWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol2vsMassCB2VWG;
        name_par_signal1v2 = name_v2Pol2vsMassCB2VWG;
        break;
    case kv2Pol2vsMassNA60QVWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol2vsMassNA60QVWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol2vsMassNA60QVWG;
        name_par_signal1v2 = name_v2Pol2vsMassNA60QVWG;
        break;
    case kv2Pol2vsMassNA60Pol2Pol3:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol2vsMassNA60Pol2Pol3, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol2vsMassNA60Pol2Pol3;
        name_par_signal1v2 = name_v2Pol2vsMassNA60Pol2Pol3;
        break;

    case kv2Pol2vsMassCB2Pol2Pol3:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol2vsMassCB2Pol2Pol3, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol2vsMassCB2Pol2Pol3;
        name_par_signal1v2 = name_v2Pol2vsMassCB2Pol2Pol3;
        break;

    case kv2Pol2vsMassCB2QVWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol2vsMassCB2QVWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol2vsMassCB2QVWG;
        name_par_signal1v2 = name_v2Pol2vsMassCB2QVWG;
        break;

    case kv2Pol2vsMassNA60VWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol2vsMassNA60VWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol2vsMassNA60VWG;
        name_par_signal1v2 = name_v2Pol2vsMassNA60VWG;
        break;

    ////////////////pol3
    case kv2Pol3vsMassCB2VWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol3vsMassCB2VWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol3vsMassCB2VWG;
        name_par_signal1v2 = name_v2Pol3vsMassCB2VWG;
        break;

    case kv2Pol3vsMassNA60QVWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol3vsMassNA60QVWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol3vsMassNA60QVWG;
        name_par_signal1v2 = name_v2Pol3vsMassNA60QVWG;
        break;
    case kv2Pol3vsMassNA60Pol2Pol3:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol3vsMassNA60Pol2Pol3, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol3vsMassNA60Pol2Pol3;
        name_par_signal1v2 = name_v2Pol3vsMassNA60Pol2Pol3;
        break;

    case kv2Pol3vsMassCB2Pol2Pol3:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol3vsMassCB2Pol2Pol3, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol3vsMassCB2Pol2Pol3;
        name_par_signal1v2 = name_v2Pol3vsMassCB2Pol2Pol3;
        break;

    case kv2Pol3vsMassCB2QVWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol3vsMassCB2QVWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol3vsMassCB2QVWG;
        name_par_signal1v2 = name_v2Pol3vsMassCB2QVWG;
        break;

    case kv2Pol3vsMassNA60VWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2Pol3vsMassNA60VWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2Pol3vsMassNA60VWG;
        name_par_signal1v2 = name_v2Pol3vsMassNA60VWG;
        break;

        /////////////////////////////////POL EXP

    case kv2PolExpvsMassCB2VWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2PolExpvsMassCB2VWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2PolExpvsMassCB2VWG;
        name_par_signal1v2 = name_v2PolExpvsMassCB2VWG;
        break;

    case kv2PolExpvsMassNA60QVWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2PolExpvsMassNA60QVWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2PolExpvsMassNA60QVWG;
        name_par_signal1v2 = name_v2PolExpvsMassNA60QVWG;
        break;
    case kv2PolExpvsMassNA60Pol2Pol3:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2PolExpvsMassNA60Pol2Pol3, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2PolExpvsMassNA60Pol2Pol3;
        name_par_signal1v2 = name_v2PolExpvsMassNA60Pol2Pol3;
        break;

    case kv2PolExpvsMassCB2Pol2Pol3:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2PolExpvsMassCB2Pol2Pol3, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2PolExpvsMassCB2Pol2Pol3;
        name_par_signal1v2 = name_v2PolExpvsMassCB2Pol2Pol3;
        break;

    case kv2PolExpvsMassCB2QVWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2PolExpvsMassCB2QVWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2PolExpvsMassCB2QVWG;
        name_par_signal1v2 = name_v2PolExpvsMassCB2QVWG;
        break;

    case kv2PolExpvsMassNA60VWG:
        // Create the distribution function with a function type and total number of parameters
        fDistrib = new TF1("fitV2", v2PolExpvsMassNA60VWG, xmin, xmax, nb_par_signalv2);

        par_signal1v2 = par_v2PolExpvsMassNA60VWG;
        name_par_signal1v2 = name_v2PolExpvsMassNA60VWG;
        break;
    }
    /////////////////////////

    // Set parameters for both backgrounds and signals

    for (int i = 0; i < nb_par_signalv2; ++i)
    {

        // Check if the current index is among the last four
        if (i >= (nb_par_bg + nb_par_signal))
        {

            fDistrib->SetParameter(i, par_signal1v2[i]);
            fDistrib->SetParName(i, name_par_signal1v2[i]);
        }
        else
        {

            // fDistrib->SetParameter(i, par_signal1v2[i]);
            fDistrib->FixParameter(i, par_signal1v2[i]); //
            fDistrib->SetParName(i, name_par_signal1v2[i]);
        }
    }

    return fDistrib;
}