void gist_ad_diff_q2_en_w_disrt () {
TCanvas *c = new TCanvas("c","c",1000,1000);
//c->Divide(4,4);
ostringstream qqq;
TH1F *h1,*h2,*h3,*h4,*h5,*h6,*h7, *h8, *h9, *h10, *h11, *h12, *h13;
gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.03,"t"); 

//TFile *MyFile1 = new TFile("out_new_1_5GeV_131_181_0325.root","READ");
//TFile *MyFile1 = new TFile("qqq_18Aug_Ebeam6_0_005_0_0051.root","READ");
TFile *MyFile1 = new TFile("out_26Aug_w_q2_dep_diff_bins_diff_sigmas_from_thresh2.root","READ");

MyFile1->cd();
qqq << "h_odn_w_dep_tot_15";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h1);
qqq.str("");

c->cd();
c->cd()->SetBottomMargin(0.1);
c->cd()->SetLeftMargin(0.09);
c->cd()->SetTopMargin(0.05);
c->cd()->SetRightMargin(0.2);

h1->GetYaxis()->SetLabelSize(0.05);
h1->GetXaxis()->SetLabelSize(0.05);

Float_t factor;
factor =1.;
//factor= factor/10000000.*100.; 
factor= factor/100000000.*100.*12.95; 

h1->Scale(factor);


h1->GetXaxis()->SetTitle("W, GeV");
h1->GetYaxis()->SetTitle("#sigma, #mub");

h1->SetLineWidth(2.);
h1->SetLineColor(kBlack);
h1->GetXaxis()->SetTitleSize(0.05);
h1->SetMaximum(73.);
//h1->SetAxisRange(0,1.4,"Y");
h1->GetYaxis()->SetTitleSize(0.05);
h1->GetYaxis()->SetNdivisions(6);
h1->GetXaxis()->SetNdivisions(6);
qqq << "E_{beam} = " << 6 <<" GeV";
h1->SetTitle(qqq.str().c_str());
qqq.str("");
h1->Draw();


//TFile *MyFile2 = new TFile("out_new_1_5GeV_131_181_0375.root","READ");
//TFile *MyFile2 = new TFile("qqq_18Aug_Ebeam6_005_01.root","READ");
//MyFile2->cd();
qqq << "h_odn_w_dep_tot_25";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h2);
qqq.str("");
h2->Scale(factor);
h2->SetLineWidth(2.);
h2->SetLineColor(kBlue);
h2->Draw("same");


//TFile *MyFile3 = new TFile("out_new_1_5GeV_131_181_0425.root","READ");
//TFile *MyFile3 = new TFile("qqq_18Aug_Ebeam6_01_015.root","READ");
//MyFile3->cd();
qqq << "h_odn_w_dep_tot_35";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h3);
qqq.str("");
h3->Scale(factor);
h3->SetLineWidth(2.);
h3->SetLineColor(kGreen);
h3->Draw("same");


//TFile *MyFile4 = new TFile("out_new_1_5GeV_131_181_0475.root","READ");
//TFile *MyFile4 = new TFile("qqq_18Aug_Ebeam6_015_02.root","READ");
//MyFile4->cd();
qqq << "h_odn_w_dep_tot_45";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h4);
qqq.str("");
h4->Scale(factor);
h4->SetLineWidth(2.);
h4->SetLineColor(kRed);
h4->Draw("same");


//TFile *MyFile4 = new TFile("out_new_1_5GeV_131_181_0525.root","READ");
//TFile *MyFile4 = new TFile("qqq_18Aug_Ebeam6_02_025.root","READ");
//MyFile4->cd();
qqq << "h_odn_w_dep_tot_55";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h5);
qqq.str("");
h5->Scale(factor);
h5->SetLineWidth(2.);
h5->SetLineColor(kMagenta);
h5->Draw("same");


//TFile *MyFile4 = new TFile("out_new_1_5GeV_141_181_0575.root","READ");
//TFile *MyFile4 = new TFile("qqq_18Aug_Ebeam6_025_035.root","READ");
//MyFile4->cd();
qqq << "h_odn_w_dep_tot_65";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h6);
qqq.str("");
h6->Scale(factor);
h6->SetLineWidth(2.);
h6->SetLineColor(kOrange);
h6->Draw("same");


//TFile *MyFile4 = new TFile("out_new_1_5GeV_141_181_0625.root","READ");
//TFile *MyFile4 = new TFile("qqq_18Aug_Ebeam6_035_045.root","READ");
//MyFile4->cd();
qqq << "h_odn_w_dep_tot_75";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h7);
qqq.str("");
h7->Scale(factor);
h7->SetLineWidth(2.);
h7->SetLineColor(kCyan);
h7->Draw("same");

//TFile *MyFile4 = new TFile("qqq_18Aug_Ebeam6_045_055.root","READ");
//MyFile4->cd();
qqq << "h_odn_w_dep_tot_85";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h8);
qqq.str("");
h8->Scale(factor);
h8->SetLineWidth(2.);
h8->SetLineColor(kYellow);
h8->Draw("same");

//TFile *MyFile4 = new TFile("qqq_18Aug_Ebeam6_055_065.root","READ");
//MyFile4->cd();
qqq << "h_odn_w_dep_tot_95";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h9);
qqq.str("");
h9->Scale(factor);
h9->SetLineWidth(2.);
h9->SetLineColor(kRed-9);
h9->Draw("same");


//TFile *MyFile4 = new TFile("qqq_18Aug_Ebeam6_08_11.root","READ");
//MyFile4->cd();
qqq << "h_odn_w_dep_tot_105";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h10);
qqq.str("");
h10->Scale(factor);
h10->SetLineWidth(2.);
h10->SetLineColor(kGreen-9);
h10->Draw("same");


//TFile *MyFile4 = new TFile("qqq_18Aug_Ebeam6_11_13.root","READ");
//MyFile4->cd();
qqq << "h_odn_w_dep_tot_115";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h11);
qqq.str("");
h11->Scale(factor);
h11->SetLineWidth(2.);
h11->SetLineColor(kMagenta+2);
h11->Draw("same");

qqq << "h_odn_w_dep_tot_125";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h12);
qqq.str("");
h12->Scale(factor);
h12->SetLineWidth(2.);
h12->SetLineColor(kCyan+3);
h12->Draw("same");

TFile *MyFile4 = new TFile("out_27Aug_w_q2_dep_0_05_from_thresh2.root","READ");
MyFile4->cd();
qqq << "W";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h13);
qqq.str("");
Float_t factor1;
factor1 =1.;

factor1= factor1/8000000.*100.; 
h13->Scale(factor1);
h13->SetLineWidth(2.);
h13->SetLineColor(kRed+2);
h13->Draw("same");




leg = new TLegend(0.1,0.3,0.4,0.4);
leg->AddEntry(h13,"Q^{2} = 0.05 GeV^{2}","l");
leg->AddEntry(h1,"Q^{2} = 0.15 GeV^{2}","l");
leg->AddEntry(h2,"Q^{2} = 0.25  GeV^{2}","l");
leg->AddEntry(h3,"Q^{2} = 0.35  GeV^{2}","l");
leg->AddEntry(h4,"Q^{2} = 0.45  GeV^{2}","l");
leg->AddEntry(h5,"Q^{2} = 0.55  GeV^{2}","l");
leg->AddEntry(h6,"Q^{2} = 0.65  GeV^{2}","l");
leg->AddEntry(h7,"Q^{2} = 0.75  GeV^{2}","l");
leg->AddEntry(h8,"Q^{2} = 0.85  GeV^{2}","l");
leg->AddEntry(h9,"Q^{2} = 0.95  GeV^{2}","l");
leg->AddEntry(h10,"Q^{2} = 1.05  GeV^{2}","l");
leg->AddEntry(h11,"Q^{2} = 1.15  GeV^{2}","l");
leg->AddEntry(h12,"Q^{2} = 1.25  GeV^{2}","l");


leg->Draw();
};
