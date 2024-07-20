
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

void plot_EEC_z() {

	bool jetnorm = 1;
	int jetptlo = 20;
	int jetpthi = jetptlo + 20;

	TFile* fin = new TFile(Form("z_%d.root",jetptlo));

	// get histograms
	
	TH1F* hpt = (TH1F*)(fin->Get("hjpt")->Clone());

	const char* label_hEEC = "hec1_2";
	TH2F* hEEC_2d = (TH2F*)(fin->Get(label_hEEC)->Clone());
	TH1D* hEEC = norm_hist(hEEC_2d, hpt, label_hEEC, jetptlo, jetpthi, jetnorm);

	// z weight
	const char* label_hEEC_z = "hec1_z_2";
	TH2F* hEEC_z_2d = (TH2F*)(fin->Get(label_hEEC_z)->Clone());
	TH1D* hEEC_z = norm_hist(hEEC_z_2d, hpt, label_hEEC_z, jetptlo, jetpthi, jetnorm);

	cout << "HISTOGRAMS LOADED" << endl;
	
	// set ratio hists
	TH1D* h_ratio_1 = (TH1D*)hEEC_z->Clone();
	TH1D* h_ratio_ref = (TH1D*)hEEC->Clone();


	TCanvas* c1 = new TCanvas("c","c",800,1000);
	gStyle->SetOptStat(0);

	// EEC comparison plot 

	TPad* p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.4);
	
	double leg_y_hi = 0.95;
	double leg_y_lo = leg_y_hi-0.1;
	TLegend* leg = new TLegend(0.1,leg_y_lo,0.3,leg_y_hi);
	format_legend(leg, "", 0.03);
	
	format_hist(hEEC, Form("EECs with #it{z}_{ch} reweighting, %d-%d GeV jets", jetptlo, jetpthi), 20, true);
	hEEC->SetMarkerColor(kBlack);
	hEEC->SetLineColor(kBlack);
	hEEC->Draw("psame");

	format_hist(hEEC_z, "", 20);
	hEEC_z->SetMarkerColor(kAzure-3);
	hEEC_z->SetLineColor(kAzure-3);
	hEEC_z->Draw("psame");
	
	leg->AddEntry(hEEC,"EEC","p");
	leg->AddEntry(hEEC_z,"#it{z} reweighted","p");
	leg->Draw("same");

	// RATIO PANEL
	TPad* p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.62, 0.1);

	TH2F* htemp = new TH2F("htemp","",10,xmin,0.6812,10,0.8,1.25);
	htemp->GetYaxis()->SetLabelSize(0.025);
	htemp->GetYaxis()->SetTitleSize(0.025);
	htemp->GetYaxis()->SetTitle("z / none");
	htemp->GetXaxis()->SetTitle("R_{L}");
	htemp->Draw();
	
	do_ratio(h_ratio_1, h_ratio_ref, "", kAzure-3);
	
	c1->SaveAs(Form("EEC_z_%d%d.pdf",jetptlo,jetpthi));
	
}