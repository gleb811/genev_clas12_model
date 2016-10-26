void gist_diff_sig_contr_q2dep () {
gStyle->SetOptStat(0);
TCanvas *c = new TCanvas("c","c",1500,1500);
c->Divide(5,4);
ostringstream qqq;
TH1F *h, *h1, *h2;
//gStyle->Reset();
//gStyle->SetPalette(1);
//gStyle->SetOptLogz(1);
Float_t W;
gStyle->SetTitleSize(0.11,"t"); 

TFile *MyFile1 = new TFile("out_26Aug_w_q2_dep_diff_bins_diff_sigmas_from_thresh2.root","READ");



for (Int_t i=0;i<=18;i++){
MyFile1->cd();
qqq << "h_odn_q2_dep_tot_" <<10000*(1.3375+0.025*i) ;
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h);
qqq.str("");


W = 1.3375+0.025*(i-1);


c->cd(i+1);
c->cd(i+1)->SetBottomMargin(0.2);
//c->cd(i+1)->SetLeftMargin(0.2);
//c->cd()->SetRightMargin(0.2);


Float_t factor=1.;
factor = factor/100000000.*100.*21.;

//factor = factor/3.1415;
//factor = factor/2./3.1415;
//factor = factor*100./2./3.1415;

//factor = factor/((W-0.938)-(0.13957+0.13957)+0.04);

//factor = factor/((W-0.13957)-(0.938+0.13957)+0.04);

h->Scale(factor);
//h->SetAxisRange((h->GetMaximum()/2.), h->GetMaximum()+(h->GetMaximum())/5., "Y");
//h->SetAxisRange(0.8, h->GetXaxis()->GetXmax(), "X");
//h->GetXaxis()->SetLimits(h->GetXaxis()->GetXmin()+0.1,h->GetXaxis()->GetXmax()-0.01);
//h1->SetMinimum(0.001);
h->SetLineColor(kBlack);
h->SetLineWidth(2);
//h->SetMaximum( h->GetMaximum()+(h->GetMaximum())/2.);
h->Draw();

h->SetTitle(" ");
h->GetYaxis()->SetLabelSize(0.12);
h->GetXaxis()->SetLabelSize(0.12);
//h->GetXaxis()->SetTitle("M23");

h->GetXaxis()->SetTitle("Q^{2}, GeV^{2}");
//h->GetYaxis()->SetTitle("M23, GeV");
h->GetYaxis()->SetTitle("#sigma, #mub");
h->GetYaxis()->SetTitleSize(0.1);
//c->cd()->SetRightMargin(0.2);
h->GetXaxis()->SetTitleSize(0.1);
//h->SetAxisRange(1.,2.,"X");
//h->SetAxisRange(0,1.4,"Y");
//h->GetYaxis()->SetTitleSize(0.08);
//h_incl->GetXaxis()->SetNdivisions(6,0);
h->GetYaxis()->SetNdivisions(6);
h->GetXaxis()->SetNdivisions(5);
qqq <<  "  W = "<<1.3375+0.025*i <<" GeV" ;
h->SetTitle(qqq.str().c_str());
qqq.str("");


qqq << "h_odn_q2_dep_t_" <<10000*(1.3375+0.025*i);
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h1);
qqq.str("");

h1->SetLineColor(kBlue);
h1->SetLineWidth(2);
h1->Scale(factor);
h1->Draw("same");

qqq << "h_odn_q2_dep_l_" <<10000*(1.3375+0.025*i);
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h2);
qqq.str("");

h2->SetLineColor(kGreen);
h2->SetLineWidth(2);
h2->Scale(factor);
h2->Draw("same");



};
};
