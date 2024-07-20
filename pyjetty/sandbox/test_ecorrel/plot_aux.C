
void format_legend(TLegend* leg, const char* header) {
	leg->SetBorderSize(0);
	leg->SetTextSize(0.031);
	leg->SetFillStyle(0);
	leg->SetMargin(0.3);
	leg->SetHeader(header);
}

void plot_aux() {
	TFile* fin_noue = new TFile("background_20.root");
	TFile* fin_ue = new TFile("background_ue_20.root");
	
	gStyle->SetOptStat(0);

	
	// number background tracks
	TCanvas* c1 = new TCanvas();
	TH1F* hbgtrk = (TH1F*)(fin_ue->Get("hbgtrks")->Clone());
	hbgtrk->SetLineWidth(2);
	hbgtrk->SetLineColor(kBlue);
	hbgtrk->Draw("same");    
	c1->Print("output/bgtrks.pdf");

   
	// jet pT spectrum
	TCanvas* c2 = new TCanvas();
	double leg_y_hi = 0.88;
	double leg_y_lo = leg_y_hi-0.1;
	double leg_x_lo = 0.6;
	double leg_x_hi = leg_x_lo+0.37;
	TLegend* leg = new TLegend(leg_x_lo,leg_y_lo,leg_x_hi,leg_y_hi);
	format_legend(leg, "jet pT spectrum");

	TH1F* hjpt = (TH1F*)(fin_ue->Get("hjpt")->Clone());
	TH1F* hjpt_sub = (TH1F*)(fin_ue->Get("hjpt_sub")->Clone());
	TH1F* hjpt_noue = (TH1F*)(fin_noue->Get("hjpt")->Clone());
	hjpt_noue->Scale(1/hjpt_noue->Integral());
	hjpt_noue->SetLineWidth(2);
	hjpt_noue->SetLineColor(kBlack);
	hjpt_noue->Draw("same");
	hjpt_sub->Scale(1/hjpt_sub->Integral());
	hjpt_sub->SetLineWidth(2);
	hjpt_sub->SetLineColor(kRed);
	hjpt_sub->Draw("same");
	hjpt->Scale(1/hjpt->Integral());
	hjpt->SetLineWidth(2);
	hjpt->SetLineColor(kBlue);
	hjpt->Draw("same");

	leg->AddEntry(hjpt,"with UE","pl");
	leg->AddEntry(hjpt_sub,"with UE, sub","pl");
	leg->AddEntry(hjpt_noue,"no UE","pl");
	leg->Draw("same");
	c2->Print("output/ptspec.pdf");

	// event multiplicity
	TCanvas* c3 = new TCanvas();
	leg_x_lo = 0.65; leg_x_hi = leg_x_lo + 0.37;
	TLegend* leg2 = new TLegend(leg_x_lo,leg_y_lo,leg_x_hi,leg_y_hi);
	format_legend(leg2, "# charged tracks");

	TH1F* hNt = (TH1F*)(fin_ue->Get("hNtracks")->Clone());
	TH1F* hNt_noue = (TH1F*)(fin_noue->Get("hNtracks")->Clone());
	hNt_noue->SetLineWidth(2);
	hNt_noue->SetLineColor(kRed);
	hNt_noue->Draw("same");
	hNt->SetLineWidth(2);
	hNt->SetLineColor(kBlue);
	hNt->Draw("same");

	leg2->AddEntry(hNt,"with UE","pl");
	leg2->AddEntry(hNt_noue,"no UE","pl");
	leg2->Draw("same");
	c3->Print("output/Ntracks.pdf");

	// jet pt vs perp cone multiplicity
	TCanvas* c4 = new TCanvas();
	TH2F* hjpt_nperp = (TH2F*)(fin_ue->Get("hjpt_nperp")->Clone());
	hjpt_nperp->Draw("colz");
	c4->Print("output/jpt_nperp.pdf");

	// jet multiplicity vs perp cone multiplicity
	TCanvas* c5 = new TCanvas();
	TH2F* hnjet_nperp = (TH2F*)(fin_ue->Get("hnjet_nperp_2040")->Clone());
	hnjet_nperp->Draw("colz");
	c5->Print("output/njet_nperp_2040.pdf");
	
}