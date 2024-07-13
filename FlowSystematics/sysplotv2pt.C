// Author: Victor Valencia
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <cmath>
#include <regex>

struct Data
{
    std::string model;
    double v2jpsi;
    double errorV2jpsi;
    double minPt; // Added to store the minimum Pt value
    double maxPt; // Added to store the maximum Pt value
};

std::vector<Data> readData(const std::string &filename)
{
    std::vector<Data> data;
    std::ifstream file(filename);
    std::regex ptRegex(R"(\[Pt([\d.]+)-([\d.]+)\])"); // Regex to find the Pt range
    std::smatch matches;

    if (!file.is_open())
    {
        std::cerr << "Failed to open the file." << std::endl;
        return data;
    }

    std::string modelLine, dataLine;
    Data d;

    while (getline(file, modelLine))
    {
        if (std::regex_search(modelLine, matches, ptRegex))
        {
            d.minPt = std::stod(matches[1].str());
            d.maxPt = std::stod(matches[2].str());
        }

        if (!getline(file, dataLine))
        {
            std::cerr << "Expected data line after model line, none found." << std::endl;
            break;
        }

        std::istringstream iss(dataLine);
        std::string temp;
        char colon, comma;

        if (!(iss >> temp >> colon >> d.v2jpsi >> comma >> temp >> colon >> d.errorV2jpsi))
        {
            std::cerr << "Failed to parse data line: " << dataLine << std::endl;
            continue;
        }

        d.model = modelLine;
        data.push_back(d);
    }

    file.close();
    std::cout << "Total data points loaded: " << data.size() << std::endl;
    return data;
}

void plotData(const std::vector<Data> &data, const std::string &canvasName)
{
    double sumV2jpsi = 0;
    double sumErrorSquared = 0;
    int Ntot = data.size();

    for (const auto &entry : data)
    {
        sumV2jpsi += entry.v2jpsi;
        sumErrorSquared += pow(entry.errorV2jpsi, 2);
    }

    double meanV2jpsi = sumV2jpsi / Ntot;
    double statError = sqrt(sumErrorSquared);

    std::cout << "Mean v2jpsi: " << meanV2jpsi << std::endl;
    std::cout << "Statistical Error: " << statError << std::endl;

    TCanvas *c1 = new TCanvas(canvasName.c_str(), "v2jpsi Plot", 1600, 400);
    c1->SetRightMargin(0.25);
    c1->SetLeftMargin(0.12);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.2);

    TH1F *frame = new TH1F("frame", ";;v_{2}^{J/#psi}", Ntot, 0, Ntot);
    frame->SetStats(false);
    frame->GetYaxis()->SetTitleOffset(0.6);
    TGraphErrors *gr = new TGraphErrors(Ntot);

    double maxV2jpsi = std::numeric_limits<double>::lowest();
    double minV2jpsi = std::numeric_limits<double>::max();

    for (int i = 0; i < Ntot; ++i)
    {
        double upperLimit = data[i].v2jpsi + data[i].errorV2jpsi;
        double lowerLimit = data[i].v2jpsi - data[i].errorV2jpsi;

        gr->SetPoint(i, i + 0.5, data[i].v2jpsi);
        gr->SetPointError(i, 0, data[i].errorV2jpsi);
        frame->GetXaxis()->SetBinLabel(i + 1, data[i].model.c_str());

        if (upperLimit > maxV2jpsi)
            maxV2jpsi = upperLimit;
        if (lowerLimit < minV2jpsi)
            minV2jpsi = lowerLimit;
    }

    frame->SetMaximum(1.1 * maxV2jpsi);
    frame->SetMinimum(1.1 * minV2jpsi);

    frame->Draw();
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1);
    gr->SetLineColor(kRed);
    gr->SetLineWidth(1);
    gr->SetMarkerColor(kRed);
    gr->Draw("P SAME");

    c1->Modified();
    c1->Update();

    double sysError = maxV2jpsi - minV2jpsi;

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.DrawLatexNDC(0.56, 0.9, Form("v_{2}^{J/#psi} = %.4f #pm %.3f (stat)  #pm %.3f (sys)", meanV2jpsi, statError, sysError));
    latex.DrawLatexNDC(0.56, 0.86, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", data.front().minPt, data.front().maxPt));

    std::string filename = Form("output_Pt_%.1f_%.1f.pdf", data.front().minPt, data.front().maxPt);

    ////// Save pt v2jpsi  val stat sys inside new directory + .pdf files

    c1->SaveAs((std::string("output_sys_pdf/") + filename).c_str());

    // c1->SaveAs(filename.c_str());

    ////// Save pt v2jpsi  val stat sys in .txt file

    //  std::ofstream outFile("finalv2jpsipt.txt");
    std::ofstream outFile("output_sys_txt/finalv2jpsiSPTOT.txt", std::ios_base::app);

    // Check if file is open
    if (outFile.is_open())
    {
        // Write the formatted title to the file
        outFile << data.front().minPt << " " << data.front().maxPt << " " << meanV2jpsi << " " << statError << " " << sysError << std::endl;

        // Close the file
        outFile.close();
        std::cout << "Data written successfully!" << std::endl;
    }
    else
    {
        std::cout << "Unable to open file!" << std::endl;
    }
}

int sysplotv2pt()
{
    std::ofstream outFile("output_sys_txt/finalv2jpsiSPTOT.txt");
    outFile << "bin_range bin_content bin_stat_error bin_syserror" << std::endl;

    const std::vector<std::string> filenames = {"input_sys_txt/output_v2jpsi_sys_pt1.txt", "input_sys_txt/output_v2jpsi_sys_pt2.txt", "input_sys_txt/output_v2jpsi_sys_pt3.txt"};
    int idx = 0;
    for (const auto &file : filenames)
    {
        auto data = readData(file);
        std::string canvasName = "canvas_" + std::to_string(idx++);
        plotData(data, canvasName);
    }
    return 0;
}
