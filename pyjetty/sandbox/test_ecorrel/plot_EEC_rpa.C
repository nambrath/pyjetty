
const double jetR = 0.4;
const double xmin = 0.01;

void format_pad(TPad* p, double top, double bottom, bool log=true) {
	p->Draw();
	p->SetFillStyle(0);
	p->SetLeftMargin(0.1);
	p->SetRightMargin(0.05);
	p->SetTopMargin(top);
	p->SetBottomMargin(bottom);
	p->cd();
	if (log) gPad->SetLogx();
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

TH1D* norm_Z(TH1F* hpt, TH1D* Z, int jetptlo, int jetpthi) {
	TH1D* hZ = (TH1D*)Z->Clone();
	int pt_binlo = hpt->FindBin(jetptlo);
	int pt_binhi = hpt->FindBin(jetpthi);
	hZ->Scale(1/hpt->Integral(pt_binlo, pt_binhi));
	return hZ;
}

void do_ratio(TH1D* h_ratio, TH1D* h_ratio_ref, const char* ytitle, Color_t c, double ymin, double ymax){
	h_ratio->Divide(h_ratio_ref);
	h_ratio->GetXaxis()->SetRangeUser(xmin,2*jetR);
	h_ratio->GetYaxis()->SetRangeUser(ymin,ymax);
	h_ratio->GetYaxis()->CenterTitle(true);
	h_ratio->GetYaxis()->SetTitle(ytitle);
	h_ratio->GetYaxis()->SetNdivisions(505);
	format_hist(h_ratio, "", 20);
	h_ratio->SetMarkerColor(c);
	h_ratio->SetLineColor(c);
	h_ratio->Draw("psame");

	int binmax = h_ratio->GetXaxis()->FindBin(2*jetR);
	float xmax = h_ratio->GetXaxis()->GetBinLowEdge(binmax);
 
	TLine *line = new TLine(xmin, 1, xmax, 1);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw();
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

void plot_EEC_rpa() {

	bool jetnorm = 1;
	int jetptlo = 20;
	int jetpthi = jetptlo + 20;

	TFile* fin = new TFile(Form("rpa_%d.root",jetptlo));
	// TFile* fin = new TFile(Form("rpa_rand_%d.root",jetptlo));

	// get histograms
	
	TH1F* hpt = (TH1F*)(fin->Get("hjpt")->Clone());
	TH1F* hpt_rpa_br = (TH1F*)(fin->Get("hjpt_rpa_br")->Clone());
	TH1F* hpt_rpa_bc = (TH1F*)(fin->Get("hjpt_rpa_bc")->Clone());
	TH1F* hpt_rpa_ar = (TH1F*)(fin->Get("hjpt_rpa_ar")->Clone());
	TH1F* hpt_rpa_ac = (TH1F*)(fin->Get("hjpt_rpa_ac")->Clone());

	const char* label_hEEC = "hec1_2";
	TH2F* hEEC_2d = (TH2F*)(fin->Get(label_hEEC)->Clone());
	TH1D* hEEC = norm_hist(hEEC_2d, hpt, label_hEEC, jetptlo, jetpthi, jetnorm);

	const char* label_hEEC_rpa_br = "hec1_rpa_br_2";
	TH2F* hEEC_rpa_br_2d = (TH2F*)(fin->Get(label_hEEC_rpa_br)->Clone());
	TH1D* hEEC_rpa_br = norm_hist(hEEC_rpa_br_2d, hpt_rpa_br, label_hEEC_rpa_br, jetptlo, jetpthi, jetnorm);

	const char* label_hEEC_rpa_bc = "hec1_rpa_bc_2";
	TH2F* hEEC_rpa_bc_2d = (TH2F*)(fin->Get(label_hEEC_rpa_bc)->Clone());
	TH1D* hEEC_rpa_bc = norm_hist(hEEC_rpa_bc_2d, hpt_rpa_bc, label_hEEC_rpa_bc, jetptlo, jetpthi, jetnorm);

	const char* label_hEEC_rpa_ar = "hec1_rpa_ar_2";
	TH2F* hEEC_rpa_ar_2d = (TH2F*)(fin->Get(label_hEEC_rpa_ar)->Clone());
	TH1D* hEEC_rpa_ar = norm_hist(hEEC_rpa_ar_2d, hpt_rpa_ar, label_hEEC_rpa_ar, jetptlo, jetpthi, jetnorm);

	const char* label_hEEC_rpa_ac = "hec1_rpa_ac_2";
	TH2F* hEEC_rpa_ac_2d = (TH2F*)(fin->Get(label_hEEC_rpa_ac)->Clone());
	TH1D* hEEC_rpa_ac = norm_hist(hEEC_rpa_ac_2d, hpt_rpa_ac, label_hEEC_rpa_ac, jetptlo, jetpthi, jetnorm);

	cout << "HISTOGRAMS LOADED" << endl;
	
	// set ratio hists
	TH1D* h_ratio_1 = (TH1D*)hEEC_rpa_br->Clone();
	TH1D* h_ratio_2 = (TH1D*)hEEC_rpa_bc->Clone();
	TH1D* h_ratio_3 = (TH1D*)hEEC_rpa_ar->Clone();
	TH1D* h_ratio_4 = (TH1D*)hEEC_rpa_ac->Clone();
	TH1D* h_ratio_ref = (TH1D*)hEEC->Clone();


	// EEC comparison plot 
	TCanvas* c1 = new TCanvas("c","c",800,1000);
	gStyle->SetOptStat(0);
	TPad* p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.4);
	double leg_y_hi = 0.97;
	double leg_y_lo = leg_y_hi-0.2;
	TLegend* leg = new TLegend(0.1,leg_y_lo,0.3,leg_y_hi);
	format_legend(leg, "", 0.03);
	
	format_hist(hEEC_rpa_bc, Form("EECs with R_{pA} efficiency, %d-%d GeV jets", jetptlo, jetpthi), 20, true);
	hEEC_rpa_bc->SetMarkerColor(kAzure-3);
	hEEC_rpa_bc->SetLineColor(kAzure-3);
	hEEC_rpa_bc->Draw("psame");
	
	format_hist(hEEC, "", 20);
	hEEC->SetMarkerColor(kBlack);
	hEEC->SetLineColor(kBlack);
	hEEC->Draw("psame");

	format_hist(hEEC_rpa_br, "", 20);
	hEEC_rpa_br->SetMarkerColor(kRed-3);
	hEEC_rpa_br->SetLineColor(kRed-3);
	hEEC_rpa_br->Draw("psame");

	format_hist(hEEC_rpa_ac, "", 20);
	hEEC_rpa_ac->SetMarkerColor(kGreen+2);
	hEEC_rpa_ac->SetLineColor(kGreen+2);
	hEEC_rpa_ac->Draw("psame");

	format_hist(hEEC_rpa_ar, "", 20);
	hEEC_rpa_ar->SetMarkerColor(kViolet-5);
	hEEC_rpa_ar->SetLineColor(kViolet-5);
	hEEC_rpa_ar->Draw("psame");

	leg->AddEntry(hEEC,"EEC","p");
	leg->AddEntry(hEEC_rpa_br,"EEC with R_{pA} br","p");
	leg->AddEntry(hEEC_rpa_bc,"EEC with R_{pA} bc","p");
	leg->AddEntry(hEEC_rpa_ar,"EEC with R_{pA} ar","p");
	leg->AddEntry(hEEC_rpa_ac,"EEC with R_{pA} ac","p");
	leg->Draw("same");

	// RATIO PANEL
	TPad* p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.62, 0.1);
	do_ratio(h_ratio_1, h_ratio_ref, "pA / pp", kRed-3, 0.8, 1.2);
	do_ratio(h_ratio_2, h_ratio_ref, "", kAzure-3, 0.8, 1.2);
	do_ratio(h_ratio_3, h_ratio_ref, "", kViolet-5, 0.8, 1.2);
	do_ratio(h_ratio_4, h_ratio_ref, "", kGreen+2, 0.8, 1.2);
	c1->SaveAs(Form("EEC_rpa_%d%d.pdf",jetptlo,jetpthi));
	

	//=============== Z plots =================//
	TH1D* hZ = (TH1D*)(fin->Get("z_ch")->Clone());
	TH1D* hZ_rpa_br = (TH1D*)(fin->Get("z_ch_rpa_br")->Clone());
	TH1D* hZ_rpa_bc = (TH1D*)(fin->Get("z_ch_rpa_bc")->Clone());
	TH1D* hZ_rpa_ar = (TH1D*)(fin->Get("z_ch_rpa_ar")->Clone());
	TH1D* hZ_rpa_ac = (TH1D*)(fin->Get("z_ch_rpa_ac")->Clone());

	hZ = norm_Z(hpt, hZ, jetptlo, jetpthi);
	hZ_rpa_br = norm_Z(hpt_rpa_br, hZ_rpa_br, jetptlo, jetpthi);
	hZ_rpa_bc = norm_Z(hpt_rpa_bc, hZ_rpa_bc, jetptlo, jetpthi);
	hZ_rpa_ar = norm_Z(hpt_rpa_ar, hZ_rpa_ar, jetptlo, jetpthi);
	hZ_rpa_ac = norm_Z(hpt_rpa_ac, hZ_rpa_ac, jetptlo, jetpthi);

	TH1D* hZ_rpa_br_ratio = (TH1D*)hZ_rpa_br->Clone();
	TH1D* hZ_rpa_bc_ratio = (TH1D*)hZ_rpa_bc->Clone();
	TH1D* hZ_rpa_ar_ratio = (TH1D*)hZ_rpa_ar->Clone();
	TH1D* hZ_rpa_ac_ratio = (TH1D*)hZ_rpa_ac->Clone();


	TCanvas* c2 = new TCanvas("c2","c2",800,1000);
	gStyle->SetOptStat(0);
	p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.45, false);
	gPad->SetLogy();
	leg_y_hi = 0.96;
	leg_y_lo = leg_y_hi-0.2;
	TLegend* legz = new TLegend(0.63,leg_y_lo,0.93,leg_y_hi);
	format_legend(legz, "", 0.03);
	
	hZ->SetTitle(Form("#it{z}_{ch} with R_{pA} efficiency, %d-%d GeV jets", jetptlo, jetpthi));
	draw_Z(hZ, kBlack);
	draw_Z(hZ_rpa_br, kRed-3);
	draw_Z(hZ_rpa_bc, kAzure-3);
	draw_Z(hZ_rpa_ar, kViolet-5);
	draw_Z(hZ_rpa_ac, kGreen+2);

	legz->AddEntry(hZ,"#it{z}_{ch}","p");
	legz->AddEntry(hZ_rpa_br,"#it{z}_{ch} with R_{pA} br","p");
	legz->AddEntry(hZ_rpa_bc,"#it{z}_{ch} with R_{pA} bc","p");
	legz->AddEntry(hZ_rpa_ar,"#it{z}_{ch} with R_{pA} ar","p");
	legz->AddEntry(hZ_rpa_ac,"#it{z}_{ch} with R_{pA} ac","p");
	legz->Draw("same");

	// RATIO PANEL
	p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.57, 0.1, false);

	do_Z_ratio(hZ_rpa_br_ratio, hZ, "pA / pp", kRed-3, 0.5, 1.3);
	do_Z_ratio(hZ_rpa_bc_ratio, hZ, "", kAzure-3, 0.5, 1.3);
	do_Z_ratio(hZ_rpa_ar_ratio, hZ, "", kViolet-5, 0.5, 1.3);
	do_Z_ratio(hZ_rpa_ac_ratio, hZ, "", kGreen+2, 0.5, 1.3);

	TLine *line = new TLine(0, 1, 1, 1);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw();

	c2->SaveAs(Form("zch_rpa_%d%d.pdf",jetptlo,jetpthi));

	//=============== Ntracks plots =================//
	TH1D* hNtracks = (TH1D*)(fin->Get("hNtracks")->Clone());
	TH1D* hNtracks_rpa_br = (TH1D*)(fin->Get("hNtracks_rpa_br")->Clone());
	TH1D* hNtracks_rpa_br_ratio = (TH1D*)hNtracks_rpa_br->Clone();
	TH1D* hNtracks_rpa_bc = (TH1D*)(fin->Get("hNtracks_rpa_bc")->Clone());
	TH1D* hNtracks_rpa_bc_ratio = (TH1D*)hNtracks_rpa_bc->Clone();

	TCanvas* c3 = new TCanvas("c3","c3",800,1000);
	gStyle->SetOptStat(0);
	p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.45, false);
	leg_y_hi = 0.96;
	leg_y_lo = leg_y_hi-0.2;
	TLegend* legN = new TLegend(0.63,leg_y_lo,0.93,leg_y_hi);
	format_legend(legN, "", 0.03);
	
	hNtracks->SetTitle(Form("# of tracks per event with R_{pA} efficiency, %d-%d GeV jets", jetptlo, jetpthi));
	draw_Z(hNtracks_rpa_br, kRed-3, 0, 100);
	draw_Z(hNtracks_rpa_bc, kAzure-3, 0, 100);
	draw_Z(hNtracks, kBlack, 0, 100);

	legN->AddEntry(hNtracks,"PYTHIA","p");
	legN->AddEntry(hNtracks_rpa_br,"R_{pA} br","p");
	legN->AddEntry(hNtracks_rpa_bc,"R_{pA} bc","p");
	legN->Draw("same");

	// RATIO PANEL
	p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.57, 0.1, false);

	do_Z_ratio(hNtracks_rpa_br_ratio, hNtracks, "pA / pp", kRed-3, 0, 5, 0, 100);
	do_Z_ratio(hNtracks_rpa_bc_ratio, hNtracks, "", kAzure-3, 0, 5, 0, 100);

	line = new TLine(0, 1, 1, 1);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw();

	c3->SaveAs(Form("ntracks_rpa_%d%d.pdf",jetptlo,jetpthi));

	//=============== RpA plots =================//
	TH1D* hyield = (TH1D*)(fin->Get("hyield")->Clone());
	TH1D* hyield_rpa_br = (TH1D*)(fin->Get("hyield_rpa_br")->Clone());
	TH1D* hyield_rpa_br_ratio = (TH1D*)hyield_rpa_br->Clone();
	TH1D* hyield_rpa_bc = (TH1D*)(fin->Get("hyield_rpa_bc")->Clone());
	TH1D* hyield_rpa_bc_ratio = (TH1D*)hyield_rpa_bc->Clone();

	TCanvas* c4 = new TCanvas("c3","c3",800,1000);
	gStyle->SetOptStat(0);
	p1 = new TPad("pad1", "", 0., 0., 1., 1.);
	format_pad(p1, 0.06, 0.45);
	gPad->SetLogy();
	leg_y_hi = 0.96;
	leg_y_lo = leg_y_hi-0.15;
	TLegend* legR = new TLegend(0.6,leg_y_lo,0.9,leg_y_hi);
	format_legend(legR, "", 0.03);
	
	hNtracks->SetTitle(Form("R_{pA}, %d-%d GeV jets", jetptlo, jetpthi));
	draw_Z(hyield_rpa_br, kRed-3, 0.5, 200);
	draw_Z(hyield_rpa_bc, kAzure-3, 0.5, 200);
	draw_Z(hyield, kBlack, 0.5, 200);

	legR->AddEntry(hNtracks,"charged hadron yield","p");
	legR->AddEntry(hNtracks_rpa_br,"ch h yield, R_{pA} br","p");
	legR->AddEntry(hNtracks_rpa_bc,"ch h yield, R_{pA} bc","p");
	legR->Draw("same");

	// RATIO PANEL
	p2 = new TPad("pad2", "", 0., 0., 1., 1.);
	format_pad(p2, 0.57, 0.1);

	do_Z_ratio(hyield_rpa_br_ratio, hyield, "RpA", kRed-3, 0.3, 1.5, 0.5, 200);
	do_Z_ratio(hyield_rpa_bc_ratio, hyield, "", kAzure-3, 0.3, 1.5, 0.5, 200);

	line = new TLine(0, 1, 200, 1);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw();

	c4->SaveAs(Form("rpa_rpa_%d%d.pdf",jetptlo,jetpthi));

}