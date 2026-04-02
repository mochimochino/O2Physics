// ============================================================
//  Step2_CrystalBallFit.C
//  ▼ Crystal Ball function (one-sided tail):
//   f(x) = N * exp(-(x-mu)^2/(2*sigma^2))        for (x-mu)/sigma > -alpha
//          N * (n/|alpha|)^n * exp(-|alpha|^2/2)
//              * (n/|alpha| - |alpha| - (x-mu)/sigma)^(-n)    otherwise
//   Parameters: [0]=N(normalization), [1]=mu(mean), [2]=sigma, [3]=alpha, [4]=n
// ============================================================
#include <iostream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMath.h"

// ------------------------------------------------------------------
// Crystal Ball function
//   par[0] = N      (normalization)
//   par[1] = mu     (mean / peak position)
//   par[2] = sigma  (core width → resolution)
//   par[3] = alpha  (tail start point, > 0)
//   par[4] = n      (tail power index, > 1)
// ------------------------------------------------------------------
double CrystalBall(double* x, double* par)
{
    double N     = par[0];
    double mu    = par[1];
    double sigma = par[2];
    double alpha = par[3];
    double n     = par[4];

    double t = (x[0] - mu) / sigma;
    if (t >= -TMath::Abs(alpha)) {
        return N * TMath::Exp(-0.5 * t * t);
    } else {
        double A = TMath::Power(n / TMath::Abs(alpha), n)
                 * TMath::Exp(-0.5 * alpha * alpha);
        double B = n / TMath::Abs(alpha) - TMath::Abs(alpha);
        return N * A * TMath::Power(B - t, -n);
    }
}

void Step2_CrystalBallFit()
{
    // ----------------------------------------------------------
    // Settings
    // ----------------------------------------------------------
    const TString inFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_merged.root";
    const TString histName = "h2Rebinned";

    const int minEntries = 30;  // minimum entries to fit
    const double fitRange = 0.12; // ±fitRange

    TFile* fIn = TFile::Open(inFile, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "[Error] Cannot open: " << inFile << std::endl;
        return;
    }
    TH2D* h2 = (TH2D*)fIn->Get(histName);
    if (!h2) {
        std::cerr << "[Error] h2Rebinned not found in " << inFile << std::endl;
        return;
    }
    h2->SetDirectory(0);
    fIn->Close();

    int nBinsX = h2->GetNbinsX();

    // ----------------------------------------------------------
    // Style
    // ----------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    const int fontCode = 42;
    for (const char* ax : {"X", "Y", "Z"}) {
        gStyle->SetLabelFont(fontCode, ax);
        gStyle->SetTitleFont(fontCode, ax);
        gStyle->SetTitleSize(0.048, ax);
        gStyle->SetLabelSize(0.042, ax);
    }

    // ----------------------------------------------------------
    // Graphs to store fit results
    // ----------------------------------------------------------
    // gr_mu  : μ (mean → bias) vs pT_true
    // gr_sig : σ (resolution)  vs pT_true
    TGraphErrors* gr_mu  = new TGraphErrors();
    TGraphErrors* gr_sig = new TGraphErrors();
    gr_mu ->SetName("gr_mu");
    gr_sig->SetName("gr_sig");
    int nPts = 0;

    // ----------------------------------------------------------
    // Fit each bin
    // ----------------------------------------------------------
    for (int i = 0; i < nBinsX; ++i) {
        TH1D* hProj = h2->ProjectionY(Form("hProj_step2_%d", i), i+1, i+1);

        double ptLo  = h2->GetXaxis()->GetBinLowEdge(i+1);
        double ptHi  = h2->GetXaxis()->GetBinUpEdge(i+1);
        double ptCen = h2->GetXaxis()->GetBinCenter(i+1);
        double ptErr = h2->GetXaxis()->GetBinWidth(i+1) / 2.0;

        if (hProj->GetEntries() < minEntries) {
            std::cout << Form("[Skip] pT bin %d (%.2f-%.2f GeV/c): entries=%.0f < %d",
                i, ptLo, ptHi, hProj->GetEntries(), minEntries) << std::endl;
            delete hProj;
            continue;
        }

        // Initial parameter estimation
        int    maxBin = hProj->GetMaximumBin();
        double peakX  = hProj->GetBinCenter(maxBin);
        double peakY  = hProj->GetMaximum();

        // 1st pass: Gauss for rough peak/sigma estimation
        TF1* fGaus = new TF1(Form("fGaus_%d", i), "gaus",
                              peakX - fitRange, peakX + fitRange);
        fGaus->SetParameters(peakY, peakX, 0.02);
        hProj->Fit(fGaus, "Q0R");
        double mu0    = fGaus->GetParameter(1);
        double sigma0 = std::abs(fGaus->GetParameter(2));
        delete fGaus;

        if (sigma0 < 1e-4 || sigma0 > 0.5) { mu0 = peakX; sigma0 = 0.03; }

        // 2nd pass: Crystal Ball fit
        // Fit range: ±3σ (estimated)
        double fitLo = mu0 - 3.0 * sigma0;
        double fitHi = mu0 + 3.0 * sigma0;

        TF1* fCB = new TF1(Form("fCB_%d", i), CrystalBall, fitLo, fitHi, 5);
        fCB->SetParNames("N", "#mu", "#sigma", "#alpha", "n");
        fCB->SetParameters(peakY, mu0, sigma0, 1.2, 5.0);
        fCB->SetParLimits(0, 0, 1e9);        // N > 0
        fCB->SetParLimits(2, 1e-5, 10.0);     // sigma > 0
        fCB->SetParLimits(3, 0.1, 10.0);      // alpha > 0
        fCB->SetParLimits(4, 1.1, 50.0);     // n > 1

       // int fitStatus = hProj->Fit(fCB, "Q0RSL");
        TFitResultPtr r = hProj->Fit(fCB, "Q0RSL");

        double mu_fit    = fCB->GetParameter(1);
        double sigma_fit = std::abs(fCB->GetParameter(2));
        double mu_err    = fCB->GetParError(1);
        double sigma_err = fCB->GetParError(2);
        // double chi2ndf   = (fCB->GetNDF() > 0) ? fCB->GetChisquare()/fCB->GetNDF() : -1;
        double chi2ndf = -1;
        if (r->IsValid() && r->Ndf() > 0) {
            chi2ndf = r->Chi2() / r->Ndf();
        }

        TCanvas* c2 = new TCanvas(Form("c2_bin%d", i),
                                  Form("Step2: Crystal Ball Fit, pT bin %d", i),
                                  700, 550);
        c2->SetLeftMargin(0.13);
        c2->SetBottomMargin(0.13);
        c2->SetTopMargin(0.05);
        c2->SetRightMargin(0.05);

        hProj->SetLineColor(kAzure+1);
        hProj->SetLineWidth(2);
        hProj->SetFillColorAlpha(kAzure+1, 0.25);
        hProj->SetFillStyle(1001);
        hProj->GetXaxis()->SetTitle("p_{T} Resolution: (p_{T}^{reco} - p_{T}^{true}) / p_{T}^{true}");
        hProj->GetYaxis()->SetTitle("Counts");
        hProj->GetXaxis()->SetRangeUser(-0.5, 0.5);
        hProj->GetYaxis()->SetRangeUser(0, hProj->GetMaximum()*1.5);
        hProj->Draw("HIST");

        fCB->SetLineColor(kRed+1);
        fCB->SetLineWidth(2);
        fCB->Draw("SAME");

        TLegend* leg = new TLegend(0.68, 0.70, 0.90, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(fontCode);
        leg->SetTextSize(0.038);
        leg->AddEntry(hProj, "Data", "f");
        leg->AddEntry(fCB,   "Crystal Ball fit", "l");
        leg->Draw();

        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(fontCode);
        tex.SetTextSize(0.036);
        tex.DrawLatex(0.16, 0.90, "#bf{p_{T} Resolution After Coherent J/#psi Cut}");
        tex.DrawLatex(0.16, 0.85, "#bf{LHC26b8} #font[52]{(MC dataset)}");
        tex.DrawLatex(0.16, 0.80, Form("#bf{p_{T}^{true} #in [%.2f, %.2f] GeV/c}", ptLo, ptHi));
        tex.DrawLatex(0.16, 0.75, Form("#bf{Entries = %.0f}", hProj->GetEntries()));
        tex.DrawLatex(0.16, 0.70, Form("#mu = %.4f #pm %.4f", mu_fit, mu_err));
        tex.DrawLatex(0.16, 0.65, Form("#sigma = %.4f #pm %.4f", sigma_fit, sigma_err));
        if (chi2ndf >= 0)
            tex.DrawLatex(0.16, 0.60, Form("#chi^{2}/ndf = %.2f", chi2ndf));
        

        TString outPng = Form("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step2_CrystalBallFit_bin%02d.png", i);
        c2->SaveAs(outPng);
        std::cout << "[Info] Saved: " << outPng
                  << Form("  mu=%.4f  sigma=%.4f  chi2/ndf=%.2f", mu_fit, sigma_fit, chi2ndf)
                  << std::endl;


        gr_mu ->SetPoint(nPts, ptCen, mu_fit);
        gr_mu ->SetPointError(nPts, ptErr, mu_err);
        gr_sig->SetPoint(nPts, ptCen, sigma_fit);
        gr_sig->SetPointError(nPts, ptErr, sigma_err);
        ++nPts;

        delete leg;
        delete fCB;
        delete c2;
        delete hProj;
    }

    // ----------------------------------------------------------
    // Save fit results
    // ----------------------------------------------------------
    TFile* fOut = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step2_fitresults.root", "RECREATE");
    gr_mu ->Write("gr_mu");
    gr_sig->Write("gr_sigma");
    fOut->Close();
    std::cout << "[Info] Saved fit results: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step2_fitresults.root" << std::endl;

    delete gr_mu;
    delete gr_sig;
    delete h2;

    std::cout << "[Done] Step 2 finished. " << nPts << " bins fitted." << std::endl;
}
