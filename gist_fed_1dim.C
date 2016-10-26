void gist_fed_1dim () {
gStyle->SetOptStat(0);
TCanvas *c = new TCanvas("c","c",1500,1500);
c->Divide(5,4);
ostringstream qqq;
TH1F *h, *h1;
TH1F *t = new TH1F ("qq","qq",100,0.1,0.11);;
Float_t W;
gStyle->SetTitleSize(0.12,"t"); 


//TFile *MyFile1 = new TFile("out_new_11GeV_131_181_0525.root","READ");
TFile *MyFile1 = new TFile("out_27Aug_1dim_0325.root","READ");

for (Int_t i=0;i<=18;i++){
MyFile1->cd();
qqq << "h_odn_theta_2_" << 10000*(1.3375+0.025*i);
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h1);
qqq.str("");


W = 1.3375+0.025*i;



c->cd(i+1);
c->cd(i+1)->SetBottomMargin(0.2);

c->cd(i+1)->SetLeftMargin(0.2);


Float_t factor=1.;
//factor = factor/20000000.*19.*100.;
//factor = factor/20000000.*100;
//factor = factor/3.1415;
//factor = factor/2./3.1415;
//factor = factor*100./2./3.1415;

factor = factor/((W-0.938)-(0.13957+0.13957)+0.04);

//factor = factor/((W-0.13957)-(0.938+0.13957)+0.04);

//h1->Scale(factor);


h1->SetMarkerStyle(20);

h1->SetAxisRange(0., h1->GetMaximum()+(h1->GetMaximum())/10., "Y");
//h1->SetAxisRange(0.8, h1->GetXaxis()->GetXmax(), "X");
//h1->SetMinimum(0.05);
//h1->GetXaxis()->SetRange(0.8,1.6);
//h1->GetXaxis()->SetLimits(h1->GetXaxis()->GetXmin()+0.1,h->GetXaxis()->GetXmax()-0.01);
//h1->Draw("EPy0 same");

//h1->Scale(2.*2.*3.1415);
h1->SetLineWidth(2);
h1->Draw();
//h->SetTitle("Q^{2} vs W");
h1->SetTitle(" ");
h1->GetYaxis()->SetLabelSize(0.09);
h1->GetXaxis()->SetLabelSize(0.11);
//h->GetXaxis()->SetTitle("M23");

h1->GetXaxis()->SetTitle("theta, rad");
//h->GetYaxis()->SetTitle("M23, GeV");
//h1->GetYaxis()->SetTitle("d#sigma/dM, #mub/GeV");
h1->GetYaxis()->SetTitleSize(0.15);
//c->cd()->SetRightMargin(0.2);
h1->GetXaxis()->SetTitleSize(0.1);
//h->SetAxisRange(1.,2.,"X");
//h->SetAxisRange(0,1.4,"Y");
//h->GetYaxis()->SetTitleSize(0.08);
//h_incl->GetXaxis()->SetNdivisions(6,0);
h1->GetYaxis()->SetNdivisions(6);
h1->GetXaxis()->SetNdivisions(5);
qqq << "W = " << 1.3375+0.025*i <<" GeV";
h1->SetTitle(qqq.str().c_str());
qqq.str("");
t->Fill(0.);
};


//c->cd(20);
//c->cd(20)->SetBottomMargin(0.4);
//c->cd(20)->SetLeftMargin(0.4);
//c->cd(20)->SetRightMargin(0.4);
//c->cd(20)->SetTopMargin(0.4);
//t->SetTitle("444");
//t->Draw();
TLegend *leg = new TLegend(0.,0.7,0.98,0.4);
leg->AddEntry(t,"Q^{2} = 0.325 GeV^{2}","l");

leg->SetTextSize(0.18);
c->cd(20);
leg->Draw();

};
