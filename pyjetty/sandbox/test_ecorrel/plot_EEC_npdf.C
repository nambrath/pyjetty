
const double jetR = 0.4;
const double xmin = 0.01;

void format_pad(TPad* p, double top, double bottom, bool logx=true) {
	p->Draw();
	p->SetFillStyle(0);
	p->SetLeftMargin(0.1);
	p->SetRightMargin(0.05);
	p->SetTopMargin(top);
	p->SetBottomMargin(bottom);
	p->cd();
	if (logx) gPad->SetLogx();
}

void format_legend(TLegend* leg, const char* header, double size=0.031) {
	leg->SetBorderSize(0);
	leg->SetTextSize(size);
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

void draw_Z(TH1D* hZ, Color_t c, double xmin=0, double xmax=1) {
	hZ->SetMarkerSize(1.);
	hZ->SetMarkerStyle(20);
	hZ->SetLineWidth(2);
	hZ->GetXaxis()->SetRangeUser(xmin,xmax);
	hZ->SetMarkerColor(c);
	hZ->SetLineColor(c);
	hZ->GetXaxis()->SetLabelOffset(999);
	hZ->GetXaxis()->SetLabelSize(0);
	hZ->Draw("pesame");
}

TH1D* norm_Z(TH1F* hpt, TH1D* Z, int jetptlo, int jetpthi) {
	TH1D* hZ = (TH1D*)Z->Clone();
	int pt_binlo = hpt->FindBin(jetptlo);
	int pt_binhi = hpt->FindBin(jetpthi);
	hZ->Scale(1/hpt->Integral(pt_binlo, pt_binhi));
	return hZ;
}

void do_Z_ratio(TH1D* hZ_rpa, TH1D* hZ, const char* ytitle, Color_t c, double ymin, double ymax, double xmin=0, double xmax=1) {
	TH1D* h_ratio = (TH1D*)hZ_rpa->Clone();
	h_ratio->Divide((TH1D*)hZ->Clone());
	h_ratio->GetXaxis()->SetRangeUser(xmin,xmax);
	h_ratio->GetYaxis()->SetRangeUser(ymin,ymax);
	h_ratio->GetYaxis()->CenterTitle(true);
	h_ratio->GetYaxis()->SetTitle(ytitle);
	h_ratio->SetTitle("");
	h_ratio->GetYaxis()->SetNdivisions(509);
	h_ratio->SetMarkerSize(1.);
	h_ratio->SetMarkerStyle(20);
	h_ratio->SetLineWidth(2);
	h_ratio->SetMarkerColor(c);
	h_ratio->SetLineColor(c);
	h_ratio->Draw("psame");
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

void do_ratio(TH1D* h_ratio, TH1D* h_ratio_ref, const char* ytitle, Color_t c, double ymin, double ymax){
	h_ratio->Divide(h_ratio_ref);
	h_ratio->GetXaxis()->SetRangeUser(xmin,2*jetR);
	h_ratio->GetYaxis()->SetRangeUser(ymin,ymax);
	h_ratio->GetYaxis()->SetTitleSize(0.02);
	h_ratio->GetYaxis()->SetLabelSize(0.025);
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

void plot_EEC_npdf() {

	bool jetnorm = 1;
	int jetptlo = 20;
	int jetpthi = jetptlo + 20;

	TFile* fin_npdf = new TFile("npdf_big_20.root"); // 500K events
	TFile* fin_pp = new TFile("z_big_20.root"); // 500K events

	double ymin[3] = {0.85, 0.8, 0.7};
	double ymax[3] = {1.2, 1.2, 1.2};

	// get histograms
	
	TH1F* hpt_npdf = (TH1F*)(fin_npdf->Get("hjpt")->Clone());
	TH1F* hpt_pp = (TH1F*)(fin_pp->Get("hjpt")->Clone());
	
	// UE plots
	
	const char* label_hEEC = "hec1_2";
	TH2F* hEEC_pp_2d = (TH2F*)(fin_pp->Get(label_hEEC)->Clone());
	TH1D* hEEC_pp = norm_hist(hEEC_pp_2d, hpt_pp, label_hEEC, jetptlo, jetpthi, jetnorm);

	const char* label_hEEC_npdf = "hec1_npdf_2";
	TH2F* hEEC_npdf_2d = (TH2F*)(fin_npdf->Get(label_hEEC)->Clone());
	TH1D* hEEC_npdf = norm_hist(hEEC_npdf_2d, hpt_npdf, label_hEEC_npdf, jetptlo, jetpthi, jetnorm);

	cout << "HISTOGRAMS LOADED" << endl;
	
	// set ratio hists
	TH1D* h_ratio_1 = (TH1D*)hEEC_npdf->Clone();
	TH1D* h_ratio_ref = (TH1D*)hEEC_pp->Clone();


	TCanvas* c1 = new TCanvas("c","c",800,1000);
	gStyle->SetOptStat(0);

	// EEC comparison plot 

	TPad* p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.4);
	
	double leg_y_hi = 0.96;
	double leg_y_lo = leg_y_hi-0.1;
	double leg_x_lo = 0.1;
	if (jetptlo > 20) leg_x_lo = 0.6;
	TLegend* leg = new TLegend(leg_x_lo,leg_y_lo,leg_x_lo+0.27,leg_y_hi);
	format_legend(leg, "");
	
	format_hist(hEEC_npdf, Form("EEC with nPDF, %d-%d GeV jets", jetptlo, jetpthi), true);
	hEEC_npdf->SetMarkerColor(kAzure-3);
	hEEC_npdf->SetLineColor(kAzure-3);
	hEEC_npdf->Draw("psame");

	format_hist(hEEC_pp, "");
	hEEC_pp->SetMarkerColor(kBlack);
	hEEC_pp->SetLineColor(kBlack);
	hEEC_pp->Draw("psame");
	
	leg->AddEntry(hEEC_pp,"EEC","p");
	leg->AddEntry(hEEC_npdf,"EEC with EPPS16","p");
	leg->Draw("same");

	// RATIO PANEL
	TPad* p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.62, 0.1);
	
	do_ratio(h_ratio_1, h_ratio_ref, "nPDF/PDF", kAzure-3, ymin[jetptlo/20-1], ymax[jetptlo/20-1]);
	
	c1->SaveAs(Form("EEC_npdf_%d%d.pdf",jetptlo,jetpthi));

	//=============== Z plots =================//
	TFile* frpa = new TFile("rpa_20.root"); // 50K events
	TFile* fsmall = new TFile("npdf_20.root"); // 50K events
	TH1D* hZ = (TH1D*)(frpa->Get("z_ch")->Clone());
	TH1D* hZ_npdf = (TH1D*)(fsmall->Get("z_ch")->Clone());

	hZ = norm_Z(hpt_pp, hZ, jetptlo, jetpthi);
	hZ_npdf = norm_Z(hpt_npdf, hZ_npdf, jetptlo, jetpthi);

	TH1D* hZ_npdf_ratio = (TH1D*)hZ_npdf->Clone();

	TCanvas* c2 = new TCanvas("c2","c2",800,1000);
	gStyle->SetOptStat(0);
	p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.45, false);
	gPad->SetLogy();
	leg_y_hi = 0.96;
	leg_y_lo = leg_y_hi-0.2;
	TLegend* legz = new TLegend(0.63,leg_y_lo,0.93,leg_y_hi);
	format_legend(legz, "", 0.03);
	
	hZ->SetTitle(Form("#it{z}_{ch} with nPDF, %d-%d GeV jets", jetptlo, jetpthi));
	draw_Z(hZ, kBlack);
	draw_Z(hZ_npdf, kAzure-3);

	legz->AddEntry(hZ,"#it{z}_{ch}","p");
	legz->AddEntry(hZ_npdf,"#it{z}_{ch} with EPPS16","p");
	legz->Draw("same");

	// RATIO PANEL
	p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.57, 0.1, false);

	do_Z_ratio(hZ_npdf_ratio, hZ, "nPDF / none", kAzure-3, 0.5, 1.3);

	TLine *line = new TLine(0, 1, 1, 1);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw();

	c2->SaveAs(Form("zch_npdf_%d%d.pdf",jetptlo,jetpthi));
	
}