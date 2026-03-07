void AnalysisPtReal() {
    // 1. Open the file
    TFile *fIn = TFile::Open("AnalysisResults.root");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Error: Cannot open AnalysisResults.root!" << std::endl;
        return;
    }
    
    // 2. Access the directory where histograms are stored by task
    TDirectoryFile *dir = (TDirectoryFile*)fIn->Get("my-upc-mass-01");
    if (!dir) {
        std::cerr << "Error: Directory 'my-upc-01' not found!" << std::endl;
        return;
    }
    
    // 3. Retrieve the 2D Histograms
    TH2D *h2Unlike = (TH2D*)dir->Get("MassVsPtUnlike");
    TH2D *h2Like   = (TH2D*)dir->Get("MassVsPtLike");
    if (!h2Unlike || !h2Like) {
        std::cerr << "Error: 2D histograms not found!" << std::endl;
        return;
    }

    //---------------------------------------------------------
    // Set the Mass Window Selection (e.g., J/psi region)
    //---------------------------------------------------------
    double massMin = 2.8; 
    double massMax = 3.4;

    // Convert mass values to bin numbers of the X-axis
    int binMin = h2Unlike->GetXaxis()->FindBin(massMin);
    int binMax = h2Unlike->GetXaxis()->FindBin(massMax);

    // 4. Project onto the Y-axis (Pt) for the specified X (Mass) bins
    // This creates 1D Pt histograms
    TH1D *hPtUnlikeRegion = h2Unlike->ProjectionY("hPtUnlikeRegion", binMin, binMax);
    TH1D *hPtLikeRegion   = h2Like->ProjectionY("hPtLikeRegion", binMin, binMax);

    // Styling
    hPtUnlikeRegion->SetLineColor(kBlack);
    hPtUnlikeRegion->SetMarkerColor(kBlack);
    hPtUnlikeRegion->SetMarkerStyle(20);
    hPtUnlikeRegion->SetMarkerSize(0.8);
    hPtUnlikeRegion->SetTitle(Form("Dimuon p_{T} (%.1f < M_{#mu#mu} < %.1f GeV/c^{2})", massMin, massMax));
    hPtUnlikeRegion->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPtUnlikeRegion->GetYaxis()->SetTitle(Form("Counts / %.3f GeV/c", hPtUnlikeRegion->GetBinWidth(1)));

    hPtLikeRegion->SetLineColor(kRed);
    hPtLikeRegion->SetFillColorAlpha(kRed, 0.3); // Semi-transparent red for background
    hPtLikeRegion->SetFillStyle(1001);
    
    // 5. Calculate the Signal (Unlike - Like)
    TH1D *hPtSignal = (TH1D*)hPtUnlikeRegion->Clone("hPtSignal");
    hPtSignal->Add(hPtLikeRegion, -1.0); // signal = Unlike - Like
    hPtSignal->SetLineColor(kBlue);
    hPtSignal->SetMarkerColor(kBlue);
    hPtSignal->SetMarkerStyle(21);
    
    //---------------------------------------------------------
    // Draw raw histograms (Unlike and Like together)
    //---------------------------------------------------------
    TCanvas *cRaw = new TCanvas("cRawPt", "Raw Dimuon pT", 800, 600);
    gPad->SetLogy(); // pT spectra typically fall exponentially, log scale is good
    
    hPtUnlikeRegion->Draw("E");
    hPtLikeRegion->Draw("HIST SAME");
    
    TLegend *leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg1->AddEntry(hPtUnlikeRegion, "Unlike Sign (+-)", "lep");
    leg1->AddEntry(hPtLikeRegion, "Like Sign (++, --)", "f");
    leg1->Draw();
    
    cRaw->SaveAs("RawPt_MassRegion.png");

    //---------------------------------------------------------
    // Draw Signal histogram (Unlike - Like)
    //---------------------------------------------------------
    TCanvas *cSig = new TCanvas("cSignalPt", "Signal Dimuon pT", 800, 600);
    gPad->SetLogy();
    
    hPtSignal->SetTitle(Form("Signal Dimuon p_{T} (%.1f < M_{#mu#mu} < %.1f GeV/c^{2})", massMin, massMax));
    // To avoid problems drawing Logy with negative bins after subtraction
    hPtSignal->SetMinimum(0.1); 
    hPtSignal->Draw("E");
    
    TLegend *leg2 = new TLegend(0.6, 0.8, 0.9, 0.9);
    leg2->AddEntry(hPtSignal, "Signal (Unlike - Like)", "lep");
    leg2->Draw();

    cSig->SaveAs("SignalPt_MassRegion.png");

    // Print some quick stats
    std::cout << "--- Yield in Region [" << massMin << ", " << massMax << "] GeV/c^2 ---" << std::endl;
    std::cout << "Unlike Sign : " << hPtUnlikeRegion->Integral() << std::endl;
    std::cout << "Like Sign   : " << hPtLikeRegion->Integral() << std::endl;
    std::cout << "Signal      : " << hPtSignal->Integral() << std::endl;
}
