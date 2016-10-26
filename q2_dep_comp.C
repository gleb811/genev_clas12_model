void q2_dep_comp () {
gStyle->SetOptStat(0);
TCanvas *c = new TCanvas("c","c",1000,1000);
//c->Divide(4,4);
ostringstream qqq;
TH1F *h,*h1;


TH1F *h_qqq065 = new TH1F ("h_qqq065","h_qqq065",150, 0, 1.35);

h_qqq065->Fill(0.001,61.07);
h_qqq065->Fill(0.65,20.92);
h_qqq065->Fill(0.95,17.43);
h_qqq065->Fill(1.3,8.05);


c->cd();
c->cd()->SetBottomMargin(0.2);
c->cd()->SetLeftMargin(0.15);
h_qqq065->SetMarkerStyle(20);
h_qqq065->SetMarkerSize(1.5);
//h_qqq065->Scale(0.15/h_qqq065->Integral(0,15));
h_qqq065->SetAxisRange(0., h_qqq065->GetMaximum()+(h_qqq065->GetMaximum())/10., "Y");
h_qqq065->SetTitle(" ");
h_qqq065->GetYaxis()->SetLabelSize(0.05);
h_qqq065->GetXaxis()->SetLabelSize(0.05);
h_qqq065->GetXaxis()->SetTitle("Q^{2}, GeV^{2}");
h_qqq065->GetYaxis()->SetTitle("#sigma, #mub");
h_qqq065->GetXaxis()->SetTitleSize(0.06);
h_qqq065->GetYaxis()->SetTitleSize(0.06);
//h_qqq065->SetTitle ("1.775 GeV < ");


h_qqq065->Draw("P");









TFile *MyFile1 = new TFile("qqq_mult.root","READ");
MyFile1->cd();
qqq << "Q2";
cout << qqq.str().c_str() <<"\n";
gDirectory->GetObject(qqq.str().c_str(),h);
qqq.str("");
h->Scale(1./13.1);
h->SetLineColor(kBlue);
h->SetLineWidth(2);
h->SetAxisRange(0., h->GetMaximum()+(h->GetMaximum())/10., "Y");
h->Draw("same ][");



leg = new TLegend(0.1,0.3,0.4,0.4);

leg->AddEntry(h,"1.775 GeV < W < 1.8 GeV","l");


leg->Draw();



};
