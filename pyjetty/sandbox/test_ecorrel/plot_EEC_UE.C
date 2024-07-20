
const double jetR = 0.4;
const double xmin = 0.01;

void format_pad(TPad* p, double top, double bottom) {
	p->Draw();
	p->SetFillStyle(0);
	p->SetLeftMargin(0.1);
	p->SetRightMargin(0.05);
	p->SetTopMargin(top);
	p->SetBottomMargin(bottom);
	p->cd();
	gPad->SetLogx();
}

void format_legend(TLegend* leg, const char* header) {
	leg->SetBorderSize(0);
	leg->SetTextSize(0.031);
	leg->SetFillStyle(0);
	leg->SetMargin(0.3);
	leg->SetHeader(header);
}

void format_hist(TH1D* h, const char* title, bool no_x = false) {
	h->SetTitle(title);
	h->SetMarkerSize(0.7);
	h->SetMarkerStyle(20);
	h->SetLineWidth(2);
	int binmax = h->GetXaxis()->FindBin(jetR);
	float xmax = h->GetXaxis()->GetBinLowEdge(binmax);
	h->GetXaxis()->SetRangeUser(xmin,xmax);
	if (no_x) {
		h->GetXaxis()->SetLabelOffset(999);
		h->GetXaxis()->SetLabelSize(0);
	}
}

TH1D* norm_hist(TH2F* EEC, TH1F* pTspec, const char* label, int jetptlo, int jetpthi, bool jetnorm=1) {
	int eec_binlo = EEC->GetXaxis()->FindBin(jetptlo);
	int pt_binlo = pTspec->FindBin(jetptlo);
	int eec_binhi = EEC->GetXaxis()->FindBin(jetpthi);
	int pt_binhi = pTspec->FindBin(jetpthi);
	TH1D* EEC_slice = EEC->ProjectionY(Form("%s_%d%d",label,jetptlo,jetptlo), eec_binlo, eec_binhi);
	if (jetnorm) {
		cout << "norm factor: " << pTspec->Integral(pt_binlo, pt_binhi) << endl;
		EEC_slice->Scale(1/pTspec->Integral(pt_binlo, pt_binhi));
	} else {
		EEC_slice->Scale(1/EEC_slice->Integral());
	}
	return EEC_slice;
}

void do_ratio(TH1D* h_ratio, TH1D* h_ratio_ref, const char* ytitle, Color_t c){
	h_ratio->Divide(h_ratio_ref);
	h_ratio->GetXaxis()->SetRangeUser(xmin,jetR);
	h_ratio->GetYaxis()->CenterTitle(true);
	h_ratio->GetYaxis()->SetTitle(ytitle);
	format_hist(h_ratio, "");
	h_ratio->SetMarkerColor(c);
	h_ratio->SetLineColor(c);
	h_ratio->GetXaxis()->SetTitle("R_{L}");
	h_ratio->Draw("psame");

	int binmax = h_ratio->GetXaxis()->FindBin(jetR);
	float xmax = h_ratio->GetXaxis()->GetBinLowEdge(binmax);
 
	TLine *line = new TLine(xmin, 1, xmax, 1);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw();
}

void plot_EEC() {

	bool jetnorm = 1;
	int jetptlo = 20;
	int jetpthi = jetptlo + 20;

	TFile* fin_noue = new TFile(Form("background_%d.root",jetptlo));
	TFile* fin_ue = new TFile(Form("background_ue_%d.root",jetptlo));

	// get histograms
	
	TH1F* hpt = (TH1F*)(fin_ue->Get("hjpt")->Clone());
	TH1F* hpt_sub = (TH1F*)(fin_ue->Get("hjpt_sub")->Clone());
	
	TH1F* hpt_noue = (TH1F*)(fin_noue->Get("hjpt")->Clone());
	TH1F* hpt_bg = (TH1F*)(fin_noue->Get("hjpt_bg")->Clone());
	TH1F* hpt_m = (TH1F*)(fin_noue->Get("hjpt_m")->Clone()); // matched
	TH1F* hpt_bg_sub = (TH1F*)(fin_noue->Get("hjpt_bg_sub")->Clone());
	
	// thermal plots
	
	// const char* label_hEEC_bg = "hec1_bg_2";
	// TH2F* hEEC_bg_2d = (TH2F*)(fin_noue->Get(label_hEEC_bg)->Clone());    
	// TH1D* hEEC_bg = norm_hist(hEEC_bg_2d, hpt_bg, label_hEEC_bg, jetptlo, jetpthi, jetnorm);
	
	// const char* label_hEEC_bg_matched = "hecm_bg_2";
	// TH2F* hEEC_bg_matched_2d = (TH2F*)(fin_noue->Get(label_hEEC_bg_matched)->Clone());
	// TH1D* hEEC_bg_matched = norm_hist(hEEC_bg_matched_2d, hpt_m, label_hEEC_bg_matched, jetptlo, jetpthi, jetnorm);

	// const char* label_hEEC_bg_sub = "hec1_bg_sub_2";
	// TH2F* hEEC_bg_sub_2d = (TH2F*)(fin_noue->Get(label_hEEC_bg_sub)->Clone());
	// TH1D* hEEC_bg_sub = norm_hist(hEEC_bg_sub_2d, hpt_bg_sub, label_hEEC_bg_sub, jetptlo, jetpthi, jetnorm);

	// const char* label_hEEC_bg_perp = "hec1_bg_perp_2";
	// TH2F* hEEC_bg_perp_2d = (TH2F*)(fin_noue->Get(label_hEEC_bg_perp)->Clone());
	// TH1D* hEEC_bg_perp = norm_hist(hEEC_bg_perp_2d, hpt_bg, label_hEEC_bg_perp, jetptlo, jetpthi, jetnorm);
	// hEEC_bg_perp->Scale(0.5);
	
	// UE plots
	
	const char* label_hEEC = "hec1_2";
	TH2F* hEEC_2d = (TH2F*)(fin_ue->Get(label_hEEC)->Clone());
	TH1D* hEEC = norm_hist(hEEC_2d, hpt, label_hEEC, jetptlo, jetpthi, jetnorm);

	const char* label_hEEC_sub = "hec1_sub_2";
	TH2F* hEEC_sub_2d = (TH2F*)(fin_ue->Get(label_hEEC_sub)->Clone());
	TH1D* hEEC_sub = norm_hist(hEEC_sub_2d, hpt_sub, label_hEEC_sub, jetptlo, jetpthi, jetnorm);

	const char* label_hEEC_perp = "hec1_perp_2";
	TH2F* hEEC_perp_2d = (TH2F*)(fin_ue->Get(label_hEEC_perp)->Clone());
	TH1D* hEEC_perp = norm_hist(hEEC_perp_2d, hpt_sub, label_hEEC_perp, jetptlo, jetpthi, jetnorm);
	hEEC_perp->Scale(0.5);
	
	TH2F* hEEC_noue_2d = (TH2F*)(fin_noue->Get("hec1_2")->Clone());
	const char* label_hEEC_noue = "hec1_2_noue";
	TH1D* hEEC_noue = norm_hist(hEEC_noue_2d, hpt_noue, label_hEEC_noue, jetptlo, jetpthi, jetnorm);

	cout << "HISTOGRAMS LOADED" << endl;
	
	// set ratio hists
	TH1D* h_ratio_1 = (TH1D*)hEEC->Clone();
	TH1D* h_ratio_2 = (TH1D*)hEEC_sub->Clone();
	TH1D* h_ratio_3 = (TH1D*)hEEC_perp->Clone();
	TH1D* h_ratio_ref = (TH1D*)hEEC_noue->Clone();


	TCanvas* c1 = new TCanvas("c","c",800,1000);
	gStyle->SetOptStat(0);

	// EEC comparison plot 

	TPad* p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.4);
	
	double leg_y_hi = 0.93;
	double leg_y_lo = leg_y_hi-0.1;
	TLegend* leg = new TLegend(0.1,leg_y_lo,0.47,leg_y_hi);
	format_legend(leg, "");
	
	format_hist(hEEC_sub, Form("EECs and background, %d-%d GeV jets", jetptlo, jetpthi), true);
	hEEC_sub->SetMarkerColor(kRed);
	hEEC_sub->SetLineColor(kRed);
	hEEC_sub->Draw("psame");

	format_hist(hEEC_perp, "");
	hEEC_perp->SetMarkerColor(kGreen+2);
	hEEC_perp->SetLineColor(kGreen+2);
	hEEC_perp->Draw("psame");
	
	format_hist(hEEC, "");
	hEEC->SetMarkerColor(kBlue);
	hEEC->SetLineColor(kBlue);
	hEEC->Draw("psame");

	format_hist(hEEC_noue, "");
	hEEC_noue->SetMarkerColor(kBlack);
	hEEC_noue->SetLineColor(kBlack);
	hEEC_noue->Draw("psame");
	
	leg->AddEntry(hEEC,"EEC with bkg","pl");
	leg->AddEntry(hEEC_sub,"EEC with bkg, rho sub","pl");
	leg->AddEntry(hEEC_perp,"EEC with bkg, perp sub","pl");
	leg->AddEntry(hEEC_noue,"EEC","pl");
	leg->Draw("same");

	// RATIO PANEL
	TPad* p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.62, 0.1);
	
	do_ratio(h_ratio_2, h_ratio_ref, "UE/none", kRed);
	do_ratio(h_ratio_1, h_ratio_ref, "", kBlue);
	do_ratio(h_ratio_3, h_ratio_ref, "", kGreen+2);
	
	c1->SaveAs(Form("EEC_UE_sub_%d%d.pdf",jetptlo,jetpthi));
	
}