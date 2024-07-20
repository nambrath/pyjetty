
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

void do_ratio(TH1D* h_ratio, TH1D* h_ratio_ref, const char* ytitle, Color_t c){
	h_ratio->Divide(h_ratio_ref);
	h_ratio->GetXaxis()->SetRangeUser(xmin,2*jetR);
	h_ratio->GetYaxis()->CenterTitle(true);
	h_ratio->GetYaxis()->SetTitle(ytitle);
	format_hist(h_ratio, "");
	h_ratio->SetMarkerColor(c);
	h_ratio->SetLineColor(c);
	h_ratio->GetXaxis()->SetTitle("R_{L}");
	h_ratio->Draw("psame");

	int binmax = h_ratio->GetXaxis()->FindBin(2*jetR);
	float xmax = h_ratio->GetXaxis()->GetBinLowEdge(binmax);
 
	TLine *line = new TLine(xmin, 1, xmax, 1);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw();
}

void plot_EEC_ptshift() {

	bool jetnorm = 1;
	int jetptlo = 20;
	int jetpthi = jetptlo + 20;

	TFile* fin = new TFile(Form("randtrk_%d.root",jetptlo));

	// get histograms
	
	TH1F* hpt = (TH1F*)(fin->Get("hjpt")->Clone());

	const char* label_hEEC = "hec1_2";
	TH2F* hEEC_2d = (TH2F*)(fin->Get(label_hEEC)->Clone());
	TH1D* hEEC = norm_hist(hEEC_2d, hpt, label_hEEC, jetptlo, jetpthi, jetnorm);
	
	TH1D* hEEC_ptdown = norm_hist(hEEC_2d, hpt, label_hEEC, jetptlo-2, jetpthi-2, jetnorm);
	TH1D* hEEC_ptup = norm_hist(hEEC_2d, hpt, label_hEEC, jetptlo+2, jetpthi+2, jetnorm);

	cout << "HISTOGRAMS LOADED" << endl;
	
	// set ratio hists
	TH1D* h_ratio_1 = (TH1D*)hEEC_ptdown->Clone();
	TH1D* h_ratio_2 = (TH1D*)hEEC_ptup->Clone();
	TH1D* h_ratio_ref = (TH1D*)hEEC->Clone();


	TCanvas* c1 = new TCanvas("c","c",800,1000);
	gStyle->SetOptStat(0);

	// EEC comparison plot 

	TPad* p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.4);
	
	double leg_y_hi = 0.93;
	double leg_y_lo = leg_y_hi-0.1;
	TLegend* leg = new TLegend(0.1,leg_y_lo,0.47,leg_y_hi);
	format_legend(leg, "");
	
	format_hist(hEEC, Form("EECs with shifted #it{p}_{T} bins, %d-%d GeV jets", jetptlo, jetpthi), true);
	hEEC->SetMarkerColor(kBlack);
	hEEC->SetLineColor(kBlack);
	hEEC->Draw("psame");

	format_hist(hEEC_ptdown, "");
	hEEC_ptdown->SetMarkerColor(kRed);
	hEEC_ptdown->SetLineColor(kRed);
	hEEC_ptdown->Draw("psame");
	
	format_hist(hEEC_ptup, "");
	hEEC_ptup->SetMarkerColor(kBlue);
	hEEC_ptup->SetLineColor(kBlue);
	hEEC_ptup->Draw("psame");
	
	leg->AddEntry(hEEC,"20-40 GeV","pl");
	leg->AddEntry(hEEC_ptdown,"18-38 GeV","pl");
	leg->AddEntry(hEEC_ptup,"22-42 GeV","pl");
	leg->Draw("same");

	// RATIO PANEL
	TPad* p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.62, 0.1);
	
	do_ratio(h_ratio_2, h_ratio_ref, "shift / none", kBlue);
	do_ratio(h_ratio_1, h_ratio_ref, "shift / none", kRed);
	
	c1->SaveAs(Form("EEC_ptshift_%d%d.pdf",jetptlo,jetpthi));
	
}