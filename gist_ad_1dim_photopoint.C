void gist_ad_1dim_photopoint () {
gStyle->SetOptStat(0);
TCanvas *c = new TCanvas("c","c",1500,1500);
c->Divide(4,2);
ostringstream qqq;
TH1F *h, *h1;
//gStyle->Reset();
//gStyle->SetPalette(1);
//gStyle->SetOptLogz(1);
Float_t W;
gStyle->SetTitleSize(0.11,"t"); 

//TFile *MyFile1 = new TFile("qqq_26Apr_1dim_golovach_comp.root","READ");
//TFile *MyFile1 = new TFile("qqq_21May_test.root","READ");
//TFile *MyFile1 = new TFile("qqq_21May_photo_1dim_good.root","READ");
TFile *MyFile1 = new TFile("qqq_26Apr2.root","READ");
//TFile *MyFile1 = new TFile("qqq_31May.root","READ");
TFile *MyFile2 = new TFile("mok_model_1dim_plots_0_golovach.root","READ");
//TFile *MyFile2 = new TFile("mok_model_1dim_plots_0_golovach_2.root","READ");


for (Int_t i=0;i<8;i++){
MyFile1->cd();
//i=7;
qqq << "h_odn_theta_2_" << 10000*(1.6125+0.025*i);
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h);
qqq.str("");

MyFile2->cd();
qqq << "h_theta_" << i+1;
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h1);
qqq.str("");

W = 1.6125+0.025*i;

for (Int_t j=1; j<=100; j++) {
h1->SetBinError(j,0);
};


c->cd(i+1);
c->cd(i+1)->SetBottomMargin(0.2);
c->cd(i+1)->SetLeftMargin(0.2);
//c->cd()->SetRightMargin(0.2);
//h->SetMinimum(0.);
cout <<h->Integral(0,100) << " qqq \n";
//h->Scale(1/h->Integral(0,100));

Float_t factor=1.;

//factor = 100./2./3.1415;

factor = factor/50000000.*15.*100.;

//factor = factor/3.1415;
//factor = factor/2./3.1415;
//factor = factor*100./2./3.1415;
factor = factor/100.;
//factor = factor/((W-0.938)-(0.13957+0.13957)+0.04);
//factor = factor/((W-0.13957)-(0.938+0.13957)+0.04);


h->Scale(factor);
//h->Scale(100./2./3.1415);
//h->Scale(1./50000000.*15./100);

//h->Scale(1/h->GetMaximum());
//h->SetAxisRange(0., h->GetMaximum()+(h->GetMaximum())/10., "Y");
//h->SetAxisRange(0.8, h->GetXaxis()->GetXmax(), "X");
//h->GetXaxis()->SetLimits(h->GetXaxis()->GetXmin()+0.1,h->GetXaxis()->GetXmax()-0.01);
//h1->SetMinimum(0.001);
h->Draw();
h1->SetMarkerStyle(20);
//cout <<h1->Integral(0,100) << " qqq1 \n";
//h1->Scale(0.175/h1->Integral(0,100));
//h1->Scale(1/h1->GetMaximum());
//h1->SetAxisRange(0., h1->GetMaximum()+(h1->GetMaximum())/10., "Y");
//h1->SetAxisRange(0.8, h1->GetXaxis()->GetXmax(), "X");
//h1->SetMinimum(0.05);
//h1->GetXaxis()->SetRange(0.8,1.6);
//h1->GetXaxis()->SetLimits(h1->GetXaxis()->GetXmin()+0.1,h->GetXaxis()->GetXmax()-0.01);
//h1->Draw("EPy0 same");


h1->Scale(2.*3.1415);
h1->Draw("P same");
//h->SetTitle("Q^{2} vs W");
h->SetTitle(" ");
h->GetYaxis()->SetLabelSize(0.05);
h->GetXaxis()->SetLabelSize(0.07);
//h->GetXaxis()->SetTitle("M23");
h->GetYaxis()->SetTitle("d#sigma/d(-cos#theta), #mub/rad");
h->GetYaxis()->SetTitleSize(0.1);
h->GetXaxis()->SetTitle("theta, rad");
//h->GetYaxis()->SetTitle("M23, GeV");

//c->cd()->SetRightMargin(0.2);
h->GetXaxis()->SetTitleSize(0.08);
//h->SetAxisRange(1.,2.,"X");
//h->SetAxisRange(0,1.4,"Y");
h->GetYaxis()->SetTitleSize(0.07);
//h_incl->GetXaxis()->SetNdivisions(6,0);
h->GetYaxis()->SetNdivisions(6);
h->GetXaxis()->SetNdivisions(5);
qqq << "W = " << 1.6125+0.025*i <<" GeV";
h->SetTitle(qqq.str().c_str());
qqq.str("");

};
};
