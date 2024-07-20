
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

void format_legend(TLegend* leg, const char* header, float size) {
	leg->SetBorderSize(0);
	leg->SetTextSize(size);
	leg->SetFillStyle(0);
	leg->SetMargin(0.3);
	leg->SetHeader(header);
}

void format_hist(TH1D* h, const char* title, int style, bool no_x = false) {
	h->SetTitle(title);
	h->SetMarkerSize(0.8);
	h->SetMarkerStyle(style);
	h->SetLineWidth(2);
	int binmax = h->GetXaxis()->FindBin(2*jetR);
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
		EEC_slice->Scale(1/pTspec->Integral(pt_binlo, pt_binhi), "width");
	} else {
		EEC_slice->Scale(1/EEC_slice->Integral(), "width");
	}
	return EEC_slice;
}

void do_ratio(TH1D* h_ratio, TH1D* h_ratio_ref, const char* ytitle, Color_t c, bool xlabel=true, bool line=true){
	h_ratio->Divide(h_ratio_ref);
	h_ratio->GetXaxis()->SetRangeUser(xmin,2*jetR);
	h_ratio->GetYaxis()->CenterTitle(true);
	h_ratio->GetYaxis()->SetTitle(ytitle);
	h_ratio->GetYaxis()->SetTitleSize(0.025);
	h_ratio->GetYaxis()->SetLabelSize(0.025);
	format_hist(h_ratio, "", 20);
	h_ratio->SetMarkerColor(c);
	h_ratio->SetLineColor(c);
	if (xlabel) {
		h_ratio->GetXaxis()->SetTitle("R_{L}");
		h_ratio->GetXaxis()->SetLabelSize(0.025);
	} else {
		h_ratio->GetXaxis()->SetLabelSize(0);
	}
	h_ratio->Draw("psame");

	int binmax = h_ratio->GetXaxis()->FindBin(2*jetR);
	float xmax = h_ratio->GetXaxis()->GetBinLowEdge(binmax);
 
	if (line) {
		TLine *line = new TLine(xmin, 1, xmax, 1);
		line->SetLineWidth(1);
		line->SetLineColor(kBlack);
		line->Draw();
	}
}

void plot_EEC_randbkg() {

	bool jetnorm = 1;
	int jetptlo = 20;
	int jetpthi = jetptlo + 20;

	TFile* fin = new TFile(Form("randtrk_%d.root",jetptlo));
	// fin = new TFile(Form("randtrk_%d_noptmod.root",jetptlo));

	// get histograms
	
	TH1F* hpt = (TH1F*)(fin->Get("hjpt")->Clone());

	const char* label_hEEC = "hec1_2";
	TH2F* hEEC_2d = (TH2F*)(fin->Get(label_hEEC)->Clone());
	TH1D* hEEC = norm_hist(hEEC_2d, hpt, label_hEEC, jetptlo, jetpthi, jetnorm);

	// 1 track added
	const char* label_hEEC_trk = "hec1_trk_2";
	TH2F* hEEC_trk_2d = (TH2F*)(fin->Get(label_hEEC_trk)->Clone());
	TH1D* hEEC_trk = norm_hist(hEEC_trk_2d, hpt, label_hEEC_trk, jetptlo, jetpthi, jetnorm);

	// 1 track added to every other jet
	const char* label_hEEC_htrk = "hec1_htrk_2";
	TH2F* hEEC_htrk_2d = (TH2F*)(fin->Get(label_hEEC_htrk)->Clone());
	TH1D* hEEC_htrk = norm_hist(hEEC_htrk_2d, hpt, label_hEEC_htrk, jetptlo, jetpthi, jetnorm);

	// 2 tracks added to every jet
	const char* label_hEEC_trk2 = "hec1_trk2_2";
	TH2F* hEEC_trk2_2d = (TH2F*)(fin->Get(label_hEEC_trk2)->Clone());
	TH1D* hEEC_trk2 = norm_hist(hEEC_trk2_2d, hpt, label_hEEC_trk2, jetptlo, jetpthi, jetnorm);
	
	// pairs with added track
	const char* label_hEEC_trkonly = "h_trk_2";
	TH2F* hEEC_trkonly_2d = (TH2F*)(fin->Get(label_hEEC_trkonly)->Clone());
	TH1D* hEEC_trkonly = norm_hist(hEEC_trkonly_2d, hpt, label_hEEC_trkonly, jetptlo, jetpthi, jetnorm);

	// pairs with 50% added track
	const char* label_hEEC_htrkonly = "h_htrk_2";
	TH2F* hEEC_htrkonly_2d = (TH2F*)(fin->Get(label_hEEC_htrkonly)->Clone());
	TH1D* hEEC_htrkonly = norm_hist(hEEC_htrkonly_2d, hpt, label_hEEC_htrkonly, jetptlo, jetpthi, jetnorm);

	// pairs with 2 added tracks
	const char* label_hEEC_2trkonly = "h_trk2_2";
	TH2F* hEEC_2trkonly_2d = (TH2F*)(fin->Get(label_hEEC_2trkonly)->Clone());
	TH1D* hEEC_2trkonly = norm_hist(hEEC_2trkonly_2d, hpt, label_hEEC_2trkonly, jetptlo, jetpthi, jetnorm);

	cout << "HISTOGRAMS LOADED" << endl;
	
	// set ratio hists
	TH1D* h_ratio_1 = (TH1D*)hEEC_trk->Clone();
	TH1D* h_ratio_2 = (TH1D*)hEEC_htrk->Clone();
	TH1D* h_ratio_3 = (TH1D*)hEEC_trk2->Clone();
	TH1D* h_ratio_ref = (TH1D*)hEEC->Clone();


	TCanvas* c1 = new TCanvas("c","c",800,1000);
	gStyle->SetOptStat(0);

	// EEC comparison plot 

	TPad* p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.4);
	
	double leg_y_hi = 0.95;
	double leg_y_lo = leg_y_hi-0.2;
	TLegend* leg = new TLegend(0.1,leg_y_lo,0.3,leg_y_hi);
	format_legend(leg, "", 0.025);
	
	format_hist(hEEC_trk2, Form("EECs with random 1 GeV tracks, %d-%d GeV jets", jetptlo, jetpthi), 20, true);
	hEEC_trk2->SetMarkerColor(kGreen+2);
	hEEC_trk2->SetLineColor(kGreen+2);
	hEEC_trk2->Draw("psame");

	format_hist(hEEC_trk, "", 20);
	hEEC_trk->SetMarkerColor(kAzure-3);
	hEEC_trk->SetLineColor(kAzure-3);
	hEEC_trk->Draw("psame");
	
	format_hist(hEEC_trkonly, "", 24);
	hEEC_trkonly->SetMarkerColor(kAzure-3);
	hEEC_trkonly->SetLineColor(kAzure-3);
	hEEC_trkonly->Draw("psame");

	format_hist(hEEC_htrk, "", 20);
	hEEC_htrk->SetMarkerColor(kRed-3);
	hEEC_htrk->SetLineColor(kRed-3);
	hEEC_htrk->Draw("psame");
	
	format_hist(hEEC_htrkonly, "", 24);
	hEEC_htrkonly->SetMarkerColor(kRed-3);
	hEEC_htrkonly->SetLineColor(kRed-3);
	hEEC_htrkonly->Draw("psame");

	format_hist(hEEC, "", 20);
	hEEC->SetMarkerColor(kBlack);
	hEEC->SetLineColor(kBlack);
	hEEC->Draw("psame");
	
	format_hist(hEEC_2trkonly, "", 24);
	hEEC_2trkonly->SetMarkerColor(kGreen+2);
	hEEC_2trkonly->SetLineColor(kGreen+2);
	hEEC_2trkonly->Draw("psame");
	
	leg->AddEntry(hEEC,"EEC","p");
	leg->AddEntry(hEEC_trk,"EEC with track","p");
	leg->AddEntry(hEEC_trkonly,"EEC from track","p");
	leg->AddEntry(hEEC_htrk,"EEC with 50% track","p");
	leg->AddEntry(hEEC_htrkonly,"EEC from 50% track","p");
	leg->AddEntry(hEEC_trk2,"EEC with 2 tracks","p");
	leg->AddEntry(hEEC_2trkonly,"EEC from 2 tracks","p");
	leg->Draw("same");

	// RATIO PANEL
	TPad* p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.61, 0.21);

	TH2F* htemp = new TH2F("htemp","",10,xmin,0.6812,10,0.8,1.5);
	htemp->GetYaxis()->SetLabelSize(0.025);
	htemp->GetXaxis()->SetLabelSize(0);
	htemp->GetYaxis()->SetTitleSize(0.025);
	htemp->GetYaxis()->SetTitle("trk / none");
	htemp->Draw();
	
	do_ratio(h_ratio_1, h_ratio_ref, "", kAzure-3, false);
	do_ratio(h_ratio_2, h_ratio_ref, "", kRed-3, false);
	do_ratio(h_ratio_3, h_ratio_ref, "", kGreen+2, false);
	
	// c1->SaveAs(Form("EEC_randtrks_%d%d.pdf",jetptlo,jetpthi));

	
	// RATIO PANEL 2
	// set ratio hists
	h_ratio_1 = (TH1D*)hEEC_htrkonly->Clone();
	h_ratio_2 = (TH1D*)hEEC_2trkonly->Clone();
	h_ratio_ref = (TH1D*)hEEC_trkonly->Clone();

	TPad* p3 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p3, 0.8, 0.1);

	htemp = new TH2F("htemp","",10,xmin,0.6812,10,0,2.5);
	htemp->GetYaxis()->SetLabelSize(0.025);
	htemp->GetYaxis()->SetTitleSize(0.025);
	htemp->GetYaxis()->SetTitle("trkonly");
	htemp->GetXaxis()->SetTitle("R_{L}");
	htemp->Draw();
	
	do_ratio(h_ratio_1, h_ratio_ref, "trkonly   ", kRed-3, true, false);
	do_ratio(h_ratio_2, h_ratio_ref, "", kGreen+2, true, false);
	do_ratio(h_ratio_ref, h_ratio_ref, "", kAzure-3, true, true);
	
	TLine *line = new TLine(xmin, 0.5, 0.6812, 0.5);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw();
	line = new TLine(xmin, 2, 0.6812, 2);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw();


	// c1->SaveAs(Form("EEC_randtrks_noptmod_%d%d.pdf",jetptlo,jetpthi));
	c1->SaveAs(Form("EEC_randtrks_%d%d.pdf",jetptlo,jetpthi));
	
}