// Author: Victor Valencia
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TColor.h>

void plotmerged(std::vector<std::string> filenames, std::vector<int> colors, std::vector<std::string> identifiers, double minCent, double maxCent)
{
    // Create a shared canvas
    TCanvas *canvas = new TCanvas("sharedCanvas", "Combined RUN2 AND RUN3", 1200, 600);

    ///////////////////////////////////////////////////// RUN 3

    // Check if the number of filenames, colors, and identifiers are the same
    if (filenames.size() != colors.size() || filenames.size() != identifiers.size())
    {
        std::cerr << "Error: The number of filenames, colors, and identifiers must be equal." << std::endl;
        return;
    }

    // Function to read data from a file and fill vectors
    auto readData = [](const std::string &filename, std::vector<double> &centers, std::vector<double> &lows, std::vector<double> &highs, std::vector<double> &contents, std::vector<double> &staterrors, std::vector<double> &syserrors)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: File " << filename << " could not be opened." << std::endl;
            return false;
        }
        std::string line;
        double low_edge, high_edge, content, staterror, syserror;
        while (getline(file, line))
        {
            std::stringstream ss(line);
            ss >> low_edge >> high_edge >> content >> staterror >> syserror;
            centers.push_back((low_edge + high_edge) / 2.0);
            lows.push_back(low_edge);
            highs.push_back(high_edge);
            contents.push_back(content);
            staterrors.push_back(staterror);
            syserrors.push_back(syserror);
        }
        file.close();
        return true;
    };

    TLegend *leg = new TLegend(0.77, 0.66, 0.85, 0.88);

    TLegend *legendInfo = new TLegend(0.15, 0.84, 0.25, 0.88);
    legendInfo->SetBorderSize(0);
    legendInfo->SetFillColor(0);
    legendInfo->SetTextSize(0.03);

    int minCentInt = static_cast<int>(minCent);
    int maxCentInt = static_cast<int>(maxCent);

    // Add the GeneralInfo to the legend
    std::string centString = "ALICE, Pb-Pb #sqrt{s} = 5.36 TeV, Inclusive J/#psi, 2.5 < y < 4, " + std::to_string(minCentInt) + "-" + std::to_string(maxCentInt) + "%";
    legendInfo->AddEntry((TObject *)0, centString.c_str(), "");

    // Process each file
    for (size_t idx = 0; idx < filenames.size(); ++idx)
    {
        std::vector<double> bin_center, bin_low, bin_high, bin_content, bin_staterror, bin_syserror;
        if (!readData(filenames[idx], bin_center, bin_low, bin_high, bin_content, bin_staterror, bin_syserror))
            continue;

        int n = bin_center.size();
        TGraphAsymmErrors *gr_stat = new TGraphAsymmErrors(n);
        TGraphAsymmErrors *gr_sys = new TGraphAsymmErrors(n);

        for (int i = 0; i < n; ++i)
        {
            gr_stat->SetPoint(i, bin_center[i], bin_content[i]);
            gr_stat->SetPointError(i, bin_center[i] - bin_low[i], bin_high[i] - bin_center[i], bin_staterror[i], bin_staterror[i]);
            gr_sys->SetPoint(i, bin_center[i], bin_content[i]);
            gr_sys->SetPointError(i, bin_center[i] - bin_low[i], bin_high[i] - bin_center[i], bin_syserror[i], bin_syserror[i]);
        }

        gr_stat->SetMarkerStyle(20);
        gr_stat->SetMarkerColor(colors[idx]);
        gr_stat->SetLineColor(colors[idx]);
        gr_sys->SetLineColor(colors[idx]);
        gr_sys->SetLineWidth(2);
        gr_sys->SetFillStyle(0);

        gr_stat->SetTitle(" ;#it{p}_{T} GeV/c ; v_{2}^{J/#psi}");
        gr_stat->GetYaxis()->SetTitleOffset(0.9);

        gr_stat->Draw(idx == 0 ? "AP" : "P same");
        gr_sys->Draw("2 same");

        leg->AddEntry(gr_stat, ("Stat - " + identifiers[idx]).c_str(), "pe");
        leg->AddEntry(gr_sys, ("Sys - " + identifiers[idx]).c_str(), "f");

        leg->SetTextSize(0.03); // Sets the text size. Adjust the value to fit your needs.
        leg->SetTextFont(42);   // Sets the font style to a standard.
        leg->SetMargin(0.2);    // Adjust the spacing between entries if necessary.
        leg->SetBorderSize(0);
    }

    ///////////////////////////////////////////////////// RUN 2

    std::ifstream filerun2("input_FinalPlot_txt/Run2Values/RUN2_v2.txt"); // Update the path as necessary
    if (!filerun2.is_open())
    {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    std::string header;
    std::getline(filerun2, header);

    const int nPointsrun2 = 4;
    double xrun2[nPointsrun2], yrun2[nPointsrun2], statErrorrun2[nPointsrun2], systErrorrun2[nPointsrun2];
    for (int i = 0; i < nPointsrun2 && !filerun2.eof(); i++)
    {
        filerun2 >> xrun2[i] >> yrun2[i] >> statErrorrun2[i] >> systErrorrun2[i];
    }
    filerun2.close();

    TGraphErrors *graphStatrun2 = new TGraphErrors(nPointsrun2, xrun2, yrun2, nullptr, statErrorrun2);
    graphStatrun2->SetTitle("Data with Errors;X-axis;Y-axis");
    graphStatrun2->SetMarkerStyle(20);
    graphStatrun2->SetMarkerColor(kGreen);
    graphStatrun2->SetLineColor(kGreen);
    graphStatrun2->Draw("P");
    TBox *boxrun2 = new TBox();
    boxrun2->SetFillStyle(0);
    boxrun2->SetFillColorAlpha(kGreen, 0.35);
    boxrun2->SetLineColor(kGreen);
    boxrun2->SetLineWidth(2);

    for (int i = 0; i < nPointsrun2; i++)
    {
        TBox *boxrun2 = new TBox(xrun2[i] - 0.1, yrun2[i] - systErrorrun2[i], xrun2[i] + 0.1, yrun2[i] + systErrorrun2[i]);
        boxrun2->SetFillStyle(0);
        boxrun2->SetFillColorAlpha(kGreen, 0.35);
        boxrun2->SetLineColor(kGreen);
        boxrun2->SetLineWidth(2);
        boxrun2->Draw("l SAME");
    }

    // Create a generic box just for the legend
    TBox *boxForLegend = new TBox();
    boxForLegend->SetFillStyle(0);
    boxForLegend->SetFillColorAlpha(kGreen, 0.35);
    boxForLegend->SetLineColor(kGreen);
    boxForLegend->SetLineWidth(2);
    leg->AddEntry(graphStatrun2, "Stat - v_{2}^{J/#psi} Run 2", "pe");
    leg->AddEntry(boxForLegend, "Sys - v_{2}^{J/#psi} Run 2", "f");

    legendInfo->Draw();
    leg->Draw();
    canvas->Modified();
    canvas->Update();

    canvas->SaveAs("output_FinalPlot_pdf/run3vsrun2comparison.pdf");
}

int plotv2jpsi_comparison()
{

    double minCent = 0.0;
    double maxCent = 50.0;

    // Setup parameters for RUN3 V2 JPSI
    std::vector<std::string> files = {"input_FinalPlot_txt/finalv2jpsiSPTOT.txt"};
    std::vector<int> colors = {kRed}; // ROOT color constants
    std::vector<std::string> identifiers = {"v_{2}^{J/#psi}{SP} Run 3"};

    plotmerged(files, colors, identifiers, minCent, maxCent);

    return 0;
}
