#include "extra_tools.cc"

parseOptions                                    options;
TH1F*                                           parseContribString(std::string tmpContribString, TH1* negHist = nullptr);
void                                            removeHistNegativeContent(TH1* procHist, TH1* negHist = nullptr);

void plotYields(std::string optFile) {

    clearHeap();

    gROOT->SetBatch();
    gStyle->SetOptStat(0);
    gStyle->SetLineScalePS(1);
    gStyle->SetCanvasPreferGL(1);
    gStyle->SetOptTitle(0);

    options.parseIt(optFile, "===");

    mkdir(options.get("writeDir"));

    Bool_t scaleByWidth                   = options.getInt("scaleByWidth");

    TLegend legend(options.getDouble("legx1"), options.getDouble("legy1"), options.getDouble("legx2"), options.getDouble("legy2"));
    legend.SetTextSize(options.getDouble("legTextSize"));
    legend.SetNColumns(options.getInt("legNcols"));
    legend.SetFillStyle(options.getInt("legFillStyle"));
    legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
    legend.SetLineColor(options.getTColFromHex("legFillColor"));
    legend.SetBorderSize(0);
    legend.SetMargin(options.getDouble("legMargin"));
    legend.SetColumnSeparation(options.getDouble("legColSep"));
    legend.SetEntrySeparation(options.getDouble("legEntrySep"));
    legend.SetCornerRadius(0.);

    std::vector<TH1F*> mainHists;

    TH1F* obs                             = parseContribString(options.get("obs"));
    Double_t obsIntegErr;
    Double_t obsInteg                     = obs->IntegralAndError(1, obs->GetNbinsX(), obsIntegErr);
    std::string obsTitle                  = "#splitline{" + std::string(obs->GetTitle()) + "}{" + std::to_string(Int_t(obsInteg)) + "}";
    // + "#pm" + to_string_with_precision(obsIntegErr, options.getInt("yieldPrec"))
    TGraphAsymmErrors* obsPoisGraph       = getPoissonErrorGraph(obs, scaleByWidth);
    obsPoisGraph->SetMarkerColor(obs->GetMarkerColor());
    obsPoisGraph->SetLineColor(obs->GetMarkerColor());
    obsPoisGraph->SetLineWidth(options.getInt("obsLineWidth"));
    obsPoisGraph->SetMarkerStyle(20);
    obsPoisGraph->SetMarkerSize(options.getFloat("obsMkrSize"));
    legend.AddEntry(obsPoisGraph, obsTitle.c_str(), "LPE");
    mainHists.push_back(obs);

    TH1F* pred                            = parseContribString(options.get("pred"));
    pred->SetLineWidth(options.getInt("totLineWidth"));
    Double_t predIntegErr;
    Double_t predInteg                    = pred->IntegralAndError(1, pred->GetNbinsX(), predIntegErr);
    predIntegErr = pred->GetBinContent(0);
    std::string predTitle                 = "#splitline{" + std::string(pred->GetTitle()) + "}{" + to_string_with_precision(predInteg, options.getInt("yieldPrec")) + "#pm" + to_string_with_precision(predIntegErr, options.getInt("yieldPrec")) + "}";
    TH1F* predFilledHist                  = (TH1F*) pred->Clone("predFilledHist");
    predFilledHist->SetFillStyle(options.getInt("totFillStyle"));
    legend.AddEntry(predFilledHist, predTitle.c_str(), "LF");
    pred->ResetAttFill();
    mainHists.push_back(pred);


    TH1F* signal                            = parseContribString(options.get("signal"));
    if(signal == nullptr) {
        signal = new TH1F();
    }
    signal->SetLineWidth(options.getInt("totLineWidth"));
    Double_t sigIntegErr;
    Double_t sigInteg                    = signal->IntegralAndError(1, signal->GetNbinsX(), sigIntegErr);
    sigIntegErr = signal->GetBinContent(0);
    Bool_t signalPresent = (sigInteg > 1E-3);
    signalPresent = 0;
    std::string sigTitle                 = "#splitline{" + std::string(signal->GetTitle()) + "}{" + to_string_with_precision(sigInteg, options.getInt("yieldPrec")) + "#pm" + to_string_with_precision(sigIntegErr, options.getInt("yieldPrec")) + "}";
    TH1F* signalFilledHist                  = (TH1F*) signal->Clone("signalFilledHist");
    signalFilledHist->SetFillStyle(options.getInt("totFillStyle"));
    if(signalPresent)legend.AddEntry(signal, sigTitle.c_str(), "L");
    signal->ResetAttFill();
    if(signalPresent) mainHists.push_back(signal);

    TGraphAsymmErrors* obsRatioHist = getPoissonErrorRatioGraph(obs, pred);
    obsRatioHist->SetMarkerColor(obs->GetMarkerColor());
    obsRatioHist->SetLineColor(obs->GetMarkerColor());
    obsRatioHist->SetLineWidth(options.getInt("obsLineWidth"));
    obsRatioHist->SetMarkerStyle(20);
    obsRatioHist->SetMarkerSize(options.getFloat("obsMkrSize"));

    if (scaleByWidth) {
        obs->Scale(1., "width");
        pred->Scale(1., "width");
        predFilledHist->Scale(1., "width");
        signal->Scale(1., "width");
    }

    TH1F* predRatio                 = (TH1F*) pred->Clone("predRatio");
    divideSelf(predRatio);
    TH1F* predRatioFilled           = (TH1F*) predFilledHist->Clone("predRatioFilled");
    divideSelf(predRatioFilled);

    TH1F* negHist                   = (TH1F*) obs->Clone("negHist");
    negHist->Reset();
    negHist->ResetAttFill();
    negHist->ResetAttLine();

    TH1F* sumHist                   = (TH1F*) negHist->Clone("sumHist");

    std::vector<TH1F*> contribHists;
    for (UInt_t i = 0; i < options.getList("breakDown", ";").size(); i++) {
        std::string iContrib = options.getList("breakDown", ";")[i];
        TH1F* iHist                        = parseContribString(iContrib, negHist);
        if (!iHist) {
            std::cout << "Warning! " << iContrib << " not found. Skipping..." << std::endl;
            delete iHist;
            continue;
        }
        if ((iHist->Integral() < 1E-2)) {
            std::cout << "Warning! " << iContrib << " has integral < 1E-2. Skipping..." << std::endl;
            delete iHist;
            continue;
        }
        contribHists.push_back(iHist);
        mainHists.push_back(iHist);
    }

    // for (Int_t i = (contribHists.size() - 1); i >= 0; i--) {
    for (UInt_t i = 0; i < contribHists.size(); i++) {
        TH1F* iHist                        = contribHists[i];
        Double_t iIntegErr;
        Double_t iInteg                    = iHist->IntegralAndError(1, iHist->GetNbinsX(), iIntegErr);
        iIntegErr = iHist->GetBinContent(0);
        std::string iTitle                 = "#splitline{" + std::string(iHist->GetTitle()) + "}{" + to_string_with_precision(iInteg, options.getInt("yieldPrec")) + "#pm" + to_string_with_precision(iIntegErr, options.getInt("yieldPrec")) + "}";
        if (scaleByWidth) iHist->Scale(1., "width");
        legend.AddEntry(iHist, iTitle.c_str(), "F");
        sumHist->Add(iHist);
    }

    if (scaleByWidth) negHist->Scale(1., "width");

    //// rescale plots to adjust for negative bin removal
    // for (Int_t iBin = 1; iBin <= pred->GetNbinsX(); iBin++) {
    //     Double_t iBinPred = pred->GetBinContent(iBin);
    //     Double_t iBinSum = sumHist->GetBinContent(iBin);
    //     Double_t iRescale = iBinPred / iBinSum;
    //     // std::cout<<"i="<<iBin<<"\trescale="<<iRescale;
    //     Double_t iBinPostRescaleSum = 0.;
    //     for (UInt_t iContrib = 0; iContrib < contribHists.size(); iContrib++) {
    //         contribHists[iContrib]->SetBinContent(iBin, iRescale * contribHists[iContrib]->GetBinContent(iBin));
    //         iBinPostRescaleSum += contribHists[iContrib]->GetBinContent(iBin);
    //     }
    //     // std::cout<<"\tpostRescaleSum="<<iBinPostRescaleSum<<"\tpred= "<<iBinPred<<std::endl;
    // }

    THStack hStack("hStack", "");
    // hStack.Add(negHist, "HIST");
    for (UInt_t i = 0; i < contribHists.size(); i++) {
        std::cout<<"Adding to stack "<<contribHists[i]->GetName()<<std::endl;
        hStack.Add(contribHists[i], "HIST");
    }

    TCanvas canvas("canvas", "", options.getDouble("canvasX"), options.getDouble("canvasY"));
    canvas.SetFillStyle(4000);

    TPad pad0("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
    pad0.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
    pad0.SetFillStyle(4000);
    pad0.SetFrameFillStyle(4000);
    pad0.SetGrid(options.getIntList("padGrid")[0], options.getIntList("padGrid")[1]);

    canvas.Draw();
    canvas.cd();

    pad0.Draw();
    pad0.cd();

    pred->Draw("hist");
    hStack.Draw("hist same");
    predFilledHist->Draw("e2 same");
    pred->Draw("hist same");
    if(signalPresent) {
        // signalFilledHist->Draw("e2 same");
        signal->Draw("hist same");
    }
    obsPoisGraph->Draw("pe0 same");

    pred->GetXaxis()->SetTitle("");
    pred->GetXaxis()->SetTitleSize(0);
    pred->GetXaxis()->SetLabelSize(0);
    pred->GetXaxis()->SetTitleOffset(0);
    pred->GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));
    pred->GetXaxis()->SetNoExponent(1);
    pred->GetYaxis()->SetTitle(scaleByWidth ? "Events/GeV" : "Events");
    pred->GetYaxis()->CenterTitle();
    pred->GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
    pred->GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
    pred->GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));
    pred->GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));

    cmsLegend(&pad0, options.get("lumitext"), options.get("CMStext"), options.getFloat("titleScale"));

    std::vector<Double_t> histMinMax = getHistMinMax(mainHists, 1, 1);
    Double_t yMin = histMinMax[0];
    Double_t yMax = histMinMax[1];
    Float_t linUnit             =   (yMax - yMin) / (1. - options.getDouble("yMinPadding") - options.getDouble("yMaxPadding"));
    if (yMin - options.getDouble("yMinPadding")*linUnit < 0) {
        yMin              =   options.getDouble("yMinRockBottom");
        linUnit                 =   yMax / (1. - options.getDouble("yMaxPadding"));
    } else {
        yMin              =   yMin - options.getDouble("yMinPadding") * linUnit;
    }
    yMax                  =   yMax + options.getDouble("yMaxPadding") * linUnit;
    pred->GetYaxis()->SetRangeUser(yMin, yMax);
    //// Ratio plot
    TPad pad1("pad1", "", options.getDouble("pad1x1"), options.getDouble("pad1y1"), options.getDouble("pad1x2"), options.getDouble("pad1y2"));
    pad1.SetMargin(options.getDouble("pad1marginL"), options.getDouble("pad1marginR"), options.getDouble("pad1marginB"), options.getDouble("pad1marginT"));
    pad1.SetFillStyle(4000);
    pad1.SetFillColor(0);
    pad1.SetFrameFillStyle(4000);
    pad1.SetGrid(1, 1);
    canvas.cd();
    pad1.Draw();
    pad1.cd();

    predRatio->Draw("hist");
    predRatioFilled->Draw("e2 same");
    predRatio->Draw("hist same");
    obsRatioHist->Draw("pe0 same");

    predRatio->GetYaxis()->SetRangeUser(options.getFloatList("ratioRange")[0], options.getFloatList("ratioRange")[1]);
    predRatio->GetXaxis()->SetTitle("V candidate jet soft-drop mass (GeV)");
    predRatio->GetXaxis()->CenterTitle();
    predRatio->GetXaxis()->SetTitleSize(options.getDouble("pad1axisTitleSize"));
    predRatio->GetXaxis()->SetLabelSize(options.getDouble("pad1axisLabelSize"));
    predRatio->GetXaxis()->SetTitleOffset(options.getDouble("pad1XtitleOffset"));
    predRatio->GetXaxis()->SetNdivisions(options.getInt("pad1xNdivs"));
    predRatio->GetXaxis()->SetNoExponent(1);
    predRatio->GetYaxis()->SetTitle("Obs/Pred");
    predRatio->GetYaxis()->CenterTitle();
    predRatio->GetYaxis()->SetTitleSize(options.getDouble("pad1axisTitleSize"));
    predRatio->GetYaxis()->SetLabelSize(options.getDouble("pad1axisLabelSize"));
    predRatio->GetYaxis()->SetTitleOffset(options.getDouble("pad1YtitleOffset"));
    predRatio->GetYaxis()->SetNdivisions(options.getInt("pad1yNdivs"));

    gPad->Update();
    gPad->Modified();
    gPad->Update();
    pad0.Update();
    pad0.Modified();
    pad0.Update();
    pad1.Modified();
    pad1.Update();
    gPad->RedrawAxis();
    gPad->RedrawAxis("G");
    pad0.RedrawAxis();
    pad0.RedrawAxis("G");
    pad1.RedrawAxis();
    pad1.RedrawAxis("G");

    pad0.cd();
    legend.Draw();

    canvas.SaveAs((options.get("writeDir") + "/" + options.get("postfix") + ".png").c_str());
    canvas.SaveAs((options.get("writeDir") + "/" + options.get("postfix") + ".pdf").c_str());

    yMin = histMinMax[0];
    yMax = histMinMax[1];
    Float_t logUnit = std::log10(yMax / yMin) / (1. - options.getDouble("yMinPaddingLog") - options.getDouble("yMaxPaddingLog"));
    yMin      = yMin / std::pow(10., options.getDouble("yMinPaddingLog") * logUnit);
    yMax      = yMax * std::pow(10., options.getDouble("yMaxPaddingLog") * logUnit);
    if (std::log10(yMax / yMin) < options.getFloat("setMoreLogLabels")) pred->GetYaxis()->SetMoreLogLabels();
    pred->SetMinimum(yMin);
    pred->SetMaximum(yMax);

    gPad->Update();
    gPad->Modified();
    gPad->Update();
    pad0.SetLogy();
    pad0.Update();
    pad0.Modified();
    pad0.Update();
    pad1.Modified();
    pad1.Update();
    gPad->RedrawAxis();
    gPad->RedrawAxis("G");
    pad0.RedrawAxis();
    pad0.RedrawAxis("G");
    pad1.RedrawAxis();
    pad1.RedrawAxis("G");

    pad0.cd();
    legend.Draw();

    canvas.SaveAs((options.get("writeDir") + "/logY_" + options.get("postfix") + ".png").c_str());
    canvas.SaveAs((options.get("writeDir") + "/logY_" + options.get("postfix") + ".pdf").c_str());

    delete obsRatioHist;
    delete obsPoisGraph;
};


TH1F*                                           parseContribString(std::string tmpContribString, TH1* negHist) {
    std::vector<std::string> contribOpts  = split_string(tmpContribString, ":");
    sanityCheck(contribOpts.size() > 3);
    TH1F* hist                            = nullptr;
    if (contribOpts[1].empty()) hist       = (TH1F*) getHistFromFile(options["tdir"] + "/" + contribOpts[0], options["sourceFile"]);
    else {
        std::vector<std::string> contribParts = split_string(contribOpts[1], ",");
        for (std::string iPart : contribParts) {
            if (!hist) {
                hist                                  = (TH1F*) getHistFromFile(options["tdir"] + "/" + iPart, options["sourceFile"]);
                hist->SetName(contribParts[0].c_str());
            } else {
                TH1F* iPartHist                       = (TH1F*) getHistFromFile(options["tdir"] + "/" + iPart, options["sourceFile"]);
                if (!iPartHist) {
                    std::cout << "Warning! " << iPart << " not found. Skipping..." << std::endl;
                    continue;
                }

                hist->Add(iPartHist);
                delete iPartHist;
            }
        }
    }
    if (!hist) return hist;
    removeHistNegativeContent(hist, negHist);
    hist->SetTitle(contribOpts[2].c_str());
    Int_t col                             = hex2rootColor(contribOpts[3]);
    hist->SetLineColor(col);
    hist->SetFillColor(col);
    hist->SetFillStyle(1001);
    hist->SetMarkerColor(col);
    return hist;
};


void           removeHistNegativeContent(TH1* procHist, TH1* negHist) {
    sanityCheck(procHist);
    Bool_t tmpNegFound = 0;
    std::string tmpStr = "Negative content found in hist" + std::string(procHist->GetName()) + "\n";
    for (Int_t iBin = 1; iBin <= procHist->GetNbinsX(); iBin++) {
        Double_t iBinErr                = procHist->GetBinError(iBin);
        Double_t iBinContent            = procHist->GetBinContent(iBin);
        if (iBinContent < 0.) {
            procHist->SetBinContent(iBin, 1E-20);
            procHist->SetBinError(iBin, 0.);
            tmpStr += std::string("\tBin ") + std::to_string(iBin) + ", content = " + to_string_with_precision(iBinContent, 3) + ", error = " + to_string_with_precision(iBinErr, 3) + "\n";
            if (negHist) {
                negHist->AddBinContent(iBin, iBinContent);
                // negHist->SetBinError(iBin, std::sqrt(std::pow(negHist->GetBinError(iBin),2)+iBinErr*iBinErr));
            }
            tmpNegFound = 1;
        }
    }
    if (tmpNegFound)             std::cout << "\t(L" << __LINE__ << ")\t" << tmpStr << std::endl;
}