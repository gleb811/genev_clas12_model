#include "TROOT.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TObject.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <iostream>
#include <THnSparse.h>

#include "inp_file_read.h"
#include "read_xsect_files.h"
#include "read_fit_param_files.h"
#include "out_file_write.h"
#include "out_file_open.h"
#include "out_file_close.h"
#include "anti_rot.h"
#include "get_xsect_ripani.h"
#include "get_xsect_near_threshold.h"

#include "get_xsect_golovach.h"
#include "get_xsect_fedotov.h"
#include "get_xsect_rip_fed_join.h"
#include "get_xsect_14_18_lowq2_fit.h"
#include "get_xsect_scale_gol_18_25.h"
#include "rot.h"
#include "global.h"
#include <stdlib.h>
#include <time.h>
#include <TLorentzVector.h>
#include <sstream>
#include <TRandom3.h>

//for BOS-files creation
/*
extern "C" 
{

#include <signal.h>
#include <errno.h>
#include <ntypes.h>
#include <bostypes.h>
#include <ec.h>
#include <clas_cern.h>
#include <ctype.h>
#include <kinematics.h>
#include <map_manager.h>
#include <trk.h>
#include <clasmdl.h>
#include <utility.h>
#include <pid.h>
#include <makebanks.h>
#include <call.h>
#include <bosddl.h>
#include <tagtnorm.h>
#include <vertex.h>

  void initbos();
  int  getBOS(BOSbank *bcs, int lun, char *list);
  void cleanBanks(BOSbank *bcs);
  void dropAllBanks(BOSbank *bcs, char *list);
  void *getBank(BOSbank *,const char *);
  void open_fpack_unit(char *filename,char *dataname,int unitnum);
  void close_fpack_unit(char *dataname);
  BOSbank bcs_ ;
  BOSbank wcs_ ;

  float ranf_();
  int ranset_( float* );
  void ranlux_( float* , const int* ) ;
 


//  time_t   time( time_t );
//  int   clock();
}

*/
 using namespace std;      
     
     
     


//Byckling function declaration
     Float_t G_BYCKLING(Float_t x, Float_t y, Float_t z, Float_t u, Float_t v, Float_t w) {
     return x*x*y+x*y*y+z*z*u+z*u*u+v*v*w+v*w*w+x*z*w+x*u*v+y*z*v+y*u*w-x*y*(z+u+v+w)-z*u*(x+y+v+w)-v*w*(x+y+z+u);
     };



int main(int argc, char** argv) {




    vector<string> args;
   
    
    Float_t sigma_t_final = 0.;
    Float_t sigma_l_final = 0.; 
    Float_t sigma_c2f_final = 0.;
    Float_t sigma_s2f_final = 0.;
    Float_t sigma_cf_final = 0.;
   Float_t  sigma_sf_final = 0.;
   
   
   Float_t sigma_t_final_1, sigma_l_final_1 , sigma_c2f_final_1, sigma_s2f_final_1,sigma_cf_final_1 ,sigma_sf_final_1 ;
   
   Float_t sigma_t_final_2, sigma_l_final_2 , sigma_c2f_final_2, sigma_s2f_final_2,sigma_cf_final_2 ,sigma_sf_final_2; 
   
   
    
    Float_t  sigma_total, sigma_total_1, sigma_total_2;
    
    
     Float_t W,Q2,phi_e;
     const Float_t phi_e_min = 0;
     const Float_t phi_e_max = 2*M_PI;
     string* file = NULL;
     
    


//This needed for taking masses of the particles from pdg_table located in ROOT_DIR


    const char *HOME_ROOT;
    const char *HOME_ROOT1;
//HOME_ROOT = getenv("ROOT");
system("root_home=`root-config --etcdir`");
HOME_ROOT = getenv("root_home");
ostringstream ROOT_DIR;
ROOT_DIR << HOME_ROOT << "/pdg_table.txt";



TDatabasePDG *pdg = new TDatabasePDG();
pdg->ReadPDGTable(ROOT_DIR.str().c_str());
TParticlePDG *part1 = new TParticlePDG();
part1 = pdg->GetParticle("proton");
MP= part1->Mass();
part1 = pdg->GetParticle("pi+");
MPIP= part1->Mass();  
part1 = pdg->GetParticle("pi-");
MPIM= part1->Mass();  

   //Reading input parameters
   inp_file_read();
   
   //Reading diff cross section from the tables in .dat files (filling out GLOBAL arrays)
   read_xsect_files();

   read_fit_param_files();
   

     
   //Reasonably changing kinematical variables if needed
    if (Q2_max > 4.*E_beam*E_beam*sin(Theta_max*M_PI/180./2.)*sin(Theta_max*M_PI/180./2.)) {
    Q2_max = 4.*E_beam*E_beam*sin(Theta_max*M_PI/180./2.)*sin(Theta_max*M_PI/180./2.);
    cout << "maximum Q2 has been changed to " << Q2_max << "\n";
    };
    
     if (Q2_min < 4.*E_beam*E_eprime_min*sin(Theta_min*M_PI/180./2.)*sin(Theta_min*M_PI/180./2.)) {
    Q2_min = 4.*E_beam*E_eprime_min*sin(Theta_min*M_PI/180./2.)*sin(Theta_min*M_PI/180./2.);
    cout << "minimum Q2 has been changed to " << Q2_min << "\n";
    };   
    
    if (W_max*W_max > (MP*MP+2.*MP*(E_beam-E_eprime_min) - Q2_min)) {
    W_max = sqrt(MP*MP+2.*MP*(E_beam-E_eprime_min)- Q2_min);
    cout << "maximum W has been changed to " << W_max << "\n";
    };
    
 if (W_min < (MP + MPIP + MPIM + 0.01)) {
    W_min = MP + MPIP + MPIM + 0.01;
    cout << "minimum W has been changed to " << W_min << "\n";
    };
    
    
    
    
   
    
    TH1F *h_W = new TH1F("W","W",22,1.3,1.85);
    TH1F *h_Q2 = new TH1F("Q2","Q2",100,Q2_min,Q2_max);
    TH1F *h_phi_e = new TH1F("phi_e","phi_e",100,phi_e_min,phi_e_max);
    TH2F *h_Q2vsW = new TH2F("Q2vsW","Q2vsWW",100,W_min,W_max,100,Q2_min,Q2_max);
    
     TH2F *h_eps_l = new TH2F("eps_l","eps_l",100,W_min,W_max,100,Q2_min,Q2_max);
    
    TH1F *h_nu = new TH1F("nu","nu",100,-1.*E_beam,E_beam);
    TH1F *h_zel = new TH1F("Z_EL","Z_EL",100,Targ_off-Targ_len/2.-1.,Targ_off+Targ_len/2.+1.);
    
    TH1F *h_inv_m12 = new TH1F("h_inv_m12","h_inv_m12",100,(MPIP+MPIM)*(MPIP+MPIM)-0.02,(1.625-MP)*(1.625-MP)+0.02);
    TH1F *h_inv_m23 = new TH1F("h_inv_m23","h_inv_m23",100,(MPIP+MP)*(MPIP+MP)-0.02,(1.625-MPIM)*(1.625-MPIM)+0.02);
     TH1F *h_th_hadr = new TH1F("h_th_hadr","h_th_hadr",100,0,M_PI);
      TH1F *h_th_hadr_2 = new TH1F("h_th_hadr_2","h_th_hadr_2",100,0,M_PI);
     TH1F *h_ph_hadr = new TH1F("h_ph_hadr","h_ph_hadr",100,0,2.*M_PI);
     TH1F *h_alph_hadr = new TH1F("h_alph_hadr","h_alph_hadr",100,0,2.*M_PI);
    
    TH2F *h_dalitz = new TH2F("dalitz","dalitz",100,MPIM+MPIP,W_max-MP,100,MPIP+MP,W_max-MPIM);
    
    TH1F *h_odn_inv_m12[37];
     TH1F *h_odn_inv_m23[37];
      TH1F *h_odn_alpha[37];
      TH1F *h_odn_theta[37]; 
     TH1F *h_odn_theta_2[37];
     
     
     
      TH1F *h_odn_wwide_inv_m12[7];
     TH1F *h_odn_wwide_inv_m23[7];
      TH1F *h_odn_wwide_alpha[7];
      TH1F *h_odn_wwide_theta[7]; 
     TH1F *h_odn_wwide_theta_2[7];
     
     TH1F *h_odn_q2_dep_t[22];
     TH1F *h_odn_q2_dep_l[22];
     TH1F *h_odn_q2_dep_l2[22];
     TH1F *h_odn_q2_dep_tot[22];
     
     TH1F *h_odn_w_dep_t[12];
      TH1F *h_odn_w_dep_l[12];
      TH1F *h_odn_w_dep_l2[12];
       TH1F *h_odn_w_dep_tot[12];
       
       ostringstream qqq;
       
      for (Short_t ii=0; ii<=11; ii++) {
      
      qqq << "h_odn_w_dep_t_" << 100*(0.15+0.1*ii);
h_odn_w_dep_t[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");

  qqq << "h_odn_w_dep_l_" << 100*(0.15+0.1*ii);
h_odn_w_dep_l[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");


  qqq << "h_odn_w_dep_l2_" << 100*(0.15+0.1*ii);
h_odn_w_dep_l2[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");



 qqq << "h_odn_w_dep_tot_" << 1000*(0.275+0.05*ii);
h_odn_w_dep_tot[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");


      
      };
      
     for (Short_t ii=0; ii<=21; ii++) { 
     
     qqq << "h_odn_q2_dep_t_" << 10000*(1.2625+0.025*ii);
h_odn_q2_dep_t[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");

qqq << "h_odn_q2_dep_l_" << 10000*(1.2625+0.025*ii);
h_odn_q2_dep_l[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");

qqq << "h_odn_q2_dep_l2_" << 10000*(1.2625+0.025*ii);
h_odn_q2_dep_l2[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");


qqq << "h_odn_q2_dep_tot_" << 10000*(1.2625+0.025*ii);
h_odn_q2_dep_tot[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");

     
     };
      
      
      
       
    
    for (Short_t ii=0; ii<=36; ii++) {
qqq << "h_odn_inv_m12_" << 10000*(1.2375+0.025*ii);
h_odn_inv_m12[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MPIM-0.02 ,1.2375+0.025*ii-MP+0.02 );
qqq.str("");

qqq << "h_odn_inv_m23_" << 10000*(1.2375+0.025*ii);
h_odn_inv_m23[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MP-0.02 ,1.2375+0.025*ii-MPIM+0.02 );
qqq.str("");

qqq << "h_odn_alpha_" << 10000*(1.2375+0.025*ii);
h_odn_alpha[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,2*M_PI );
qqq.str("");

qqq << "h_odn_theta_" << 10000*(1.2375+0.025*ii);
h_odn_theta[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

qqq << "h_odn_theta_2_" << 10000*(1.2375+0.025*ii);
h_odn_theta_2[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

//qqq << "h_odn_q2_dep_" << 10000*(1.4375+0.025*ii);
///h_odn_q2_dep[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
//qqq.str("");


};


 for (Short_t ii=0; ii<7; ii++) {
qqq << "h_odn_wwide_inv_m12_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_inv_m12[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MPIM-0.02 ,2.1875+0.05*ii-MP+0.02 );
qqq.str("");

qqq << "h_odn_wwide_inv_m23_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_inv_m23[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MP-0.02 ,2.1875+0.05*ii-MPIM+0.02 );
qqq.str("");

qqq << "h_odn_wwide_alpha_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_alpha[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,2*M_PI );
qqq.str("");

qqq << "h_odn_wwide_theta_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_theta[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

qqq << "h_odn_wwide_theta_2_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_theta_2[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

//qqq << "h_odn_q2_dep_" << 10000*(1.4375+0.025*ii);
///h_odn_q2_dep[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
//qqq.str("");


};






    
    
    Float_t nu,E_E_prime,Theta_e_prime;
    Float_t M1,M2,M3;
    TLorentzVector P4_E_prime, P4_PIP,P4_Pfin,P4_PIM;
    Float_t alph_const = 1./137.;
   
    
    srand (time(NULL));
    Int_t k=0;
    Float_t z_EL;
   
   Float_t inv_m12,inv_m23,inv_m13,th_hadr,alph_hadr,ph_hadr,s12,s23; 
   Float_t W_tmp,dummy;
   Int_t dummy_int,ii;
   
     //FOR 1-PIM, 2-PIP, 3-P
     M1 = MPIM;
     M2 = MPIP;
     M3 = MP; 
     //FOR 1-P, 2-PIP, 3-PIM
//     M1 = MP;
//    M2 = MPIP;
//     M3 = MPIM;
     //FOR 1-PIP, 2-PIM, 3-P
 //    M1 = MPIP;
 //    M2 = MPIM;
 //    M3 = MP;
    
  //open input file
     out_file_open();
     
     TRandom3 ph_e_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 th_hadr_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 W_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 Q2_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 z_EL_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 alph_hadr_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.)); 
     TRandom3 s12_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 s23_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));     
     TRandom3 ph_hadr_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));          
    
     THnSparseD *h_empty1[22][22];
     THnSparseD *h_empty2[22][22];
     THnSparseD *h_empty3[22][22]; 
     
     Int_t *bins = new Int_t[5];
     Float_t W_bin[22];
     Float_t Q2_bin[22];
     Double_t Var_1[5],Var_2[5],Var_3[5];  
     
  //  Double_t xmin[5] = {1.028, 0.229,0.,0.,0.};
   Double_t xmin[5] = {(0.938272 + 0.13957), (0.13957 + 0.13957),0.,0.,0.};
  
    Double_t xmax[5];
     xmax[2] = 180.;
     xmax[3] = 360.;
     xmax[4] = 360.;
    
    Short_t m12m12_max;
    Short_t m23m23_max;
    Short_t thth_max;
    Short_t alphalph_max; 
     Short_t phph_max; 
     
 
  // Start to generate electrons    
    for (Int_t qq2=0; qq2<22; qq2++) {
    Q2_bin[qq2] = 0.275+0.05*qq2;
    for (Int_t ww=0; ww<22; ww++) { 
    
     W_bin[ww] = 1.3125+0.025*ww; 

 //Warning alpa and phi angles are mix up
     if ((ww>=0)&&(ww<=1)){
m12m12_max=8; //8
m23m23_max=8; //8
thth_max=6; //6
phph_max=5;
alphalph_max=5; // 5
};
if ((ww>=2)&&(ww<=3)){
m12m12_max=10;
m23m23_max=10;
thth_max=8;
phph_max=5;
alphalph_max=6;
};

if ((ww>=4)&&(ww<=6)){
m12m12_max=12;
m23m23_max=12;
thth_max=10;
phph_max=5;
alphalph_max=8;
};

if ((ww>=7)&&(ww<=14)){
m12m12_max=12;
m23m23_max=12;
thth_max=10;
phph_max=8;
alphalph_max=8;
};

if ((ww>=15)&&(ww<=21)){
m12m12_max=12;
m23m23_max=12;
thth_max=10;
phph_max=8;
alphalph_max=8;
};

     
    bins[0]=m12m12_max;
    bins[1]=m23m23_max;
    bins[2]=thth_max;
    bins[3]=phph_max;
    bins[4]=alphalph_max;
   
    
   //  xmax[0] =  W_bin[ww]+0.0125 - 0.13957 + 0.05;
   //  xmax[1] =  W_bin[ww]+0.0125 - 0.938272 + 0.05;
   
   xmax[0] =  (1.3125+0.025*ww - 0.13957)+((1.3125+0.025*ww- 0.13957)-(0.938272 + 0.13957))/(bins[0]-1);
   xmax[1] = (1.3125+0.025*ww - 0.938272)+((1.3125+0.025*ww - 0.938272)-(0.13957 + 0.13957))/(bins[1]-1);
     
    qqq.str("");
    qqq <<"h_5dim_empty1_" << Q2_bin[qq2]*1000 << "_w_" << 1000*W_bin[ww];
    h_empty1[qq2][ww] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),5,bins,xmin,xmax);
    qqq.str(""); 
    qqq <<"h_5dim_empty2_" << Q2_bin[qq2]*1000 << "_w_" << 1000*W_bin[ww];
    h_empty2[qq2][ww] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),5,bins,xmin,xmax);
    qqq.str(""); 
    qqq <<"h_5dim_empty3_" << Q2_bin[qq2]*1000 << "_w_" << 1000*W_bin[ww];
    h_empty3[qq2][ww] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),5,bins,xmin,xmax);
    qqq.str(""); 
    for (Int_t m12m12=0; m12m12<m12m12_max; m12m12++) {   
    Var_2[0] = xmin[0] + (xmax[0] - xmin[0])/m12m12_max/2.+ ( xmax[0] - xmin[0])/m12m12_max*m12m12;
    for (Int_t m23m23=0; m23m23<m23m23_max; m23m23++) { 
    Var_2[1] = xmin[1] + (xmax[1] - xmin[1])/m23m23_max/2.+ ( xmax[1] - xmin[1])/m23m23_max*m23m23;
    for (Int_t thth=0; thth<thth_max; thth++) {    
     Var_2[2] = xmin[2] + (xmax[2] - xmin[2])/thth_max/2.+ ( xmax[2] - xmin[2])/thth_max*thth; 
    for (Int_t phph=0; phph<phph_max; phph++) {     
     Var_2[3] = xmin[3] + (xmax[3] - xmin[3])/phph_max/2.+ ( xmax[3] - xmin[3])/phph_max*phph;    
    for (Int_t alphalph=0; alphalph<alphalph_max; alphalph++) {   
     Var_2[4] = xmin[4] + (xmax[4] - xmin[4])/alphalph_max/2.+ ( xmax[4] - xmin[4])/alphalph_max*alphalph; 
           

phi_e =ph_e_rndm.Uniform(phi_e_min,phi_e_max);

z_EL = z_EL_rndm.Uniform(Targ_off  - Targ_len/2.,Targ_off  + Targ_len/2.);
ph_hadr = Var_2[3]*M_PI/180.;

alph_hadr =  Var_2[4]*M_PI/180.; 

    k++;

    W = W_bin[ww];

    Q2 = Q2_bin[qq2];
    
    nu=(W*W+Q2-MP*MP)/2./MP;
    if (nu <= (E_beam-E_eprime_min))  {
    
   
    E_E_prime=E_beam-nu;
    Theta_e_prime = acos(1.-Q2/E_beam/E_E_prime/2.);
    P4_E_prime.SetXYZT(E_E_prime*cos(phi_e)*sin(Theta_e_prime),E_E_prime*sin(phi_e)*sin(Theta_e_prime),E_E_prime*cos(Theta_e_prime),E_E_prime);
  //  W_tmp = W;
   // W = 1.6125;
    dummy_int = 0;
   s12 = Var_2[1]*Var_2[1];
   s23 = Var_2[0]*Var_2[0];
   inv_m12 = sqrt(s12);
   inv_m23 = sqrt(s23); 
   dummy_int++;

    // cout << W << " " << Q2 << " " << inv_m12 << "  " << inv_m23 << "  " << sigma_total << "\n";
//   } while (G_BYCKLING(inv_m12*inv_m12,inv_m23*inv_m23,W*W,M2*M2,M1*M1,M3*M3) > -5e-07);

     //FOR 1-PIM, 2-PIP, 3-P
     M1 = MPIM;
     M2 = MPIP;
     M3 = MP; 
     //FOR 1-P, 2-PIP, 3-PIM
//     M1 = MP;
//    M2 = MPIP;
//     M3 = MPIM;
     //FOR 1-PIP, 2-PIM, 3-P
 //    M1 = MPIP;
 //    M2 = MPIM;
 //    M3 = MP;

   if (G_BYCKLING(inv_m12*inv_m12,inv_m23*inv_m23,W*W,M2*M2,M1*M1,M3*M3) <= -5e-07) {
   
//   h_inv_m12 ->Fill(s12,1.);
   

	th_hadr = Var_2[2]*M_PI/180.;
	
 
//Getting cross section in given generated (W, Q2, s12, s23, theta, alpha)-point

 if ((W>=1.4125)&&(W<=1.8125)&&(Q2>=0.65)&&(Q2<=1.3)) {
 get_xsect_ripani(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);

};



 if ((W>=1.6125)&&(W<=1.8125)&&(Q2>0.000001)&&(Q2<0.65)) {

 get_xsect_14_18_lowq2_fit(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);
};



if (((W>=1.3125)&&(W<=1.4375)&&(Q2>0.275)&&(Q2<0.575))||((W>=1.4125)&&(W<=1.4875)&&(Q2>0.275)&&(Q2<0.525))||((W>=1.4125)&&(W<=1.5125)&&(Q2>0.275)&&(Q2<0.425))||((W>=1.4125)&&(W<=1.5375)&&(Q2>0.225)&&(Q2<0.375))||((W>=1.4125)&&(W<=1.5625)&&(Q2>0.225)&&(Q2<0.275))){

 get_xsect_fedotov(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);


};

 if (((W>=1.4125)&&(W<=1.4375)&&(Q2>=0.575)&&(Q2<=0.65))||((W>=1.4375)&&(W<=1.4875)&&(Q2>=0.525)&&(Q2<=0.65))||((W>=1.4875)&&(W<=1.5125)&&(Q2>=0.425)&&(Q2<=0.65))||((W>=1.5125)&&(W<=1.5375)&&(Q2>=0.325)&&(Q2<=0.65))||((W>=1.5375)&&(W<=1.5625)&&(Q2>=0.275)&&(Q2<=0.65))||((W>=1.5625)&&(W<=1.5875)&&(Q2>=0.225)&&(Q2<=0.65))||((W>=1.3125)&&(W<1.4125)&&(Q2>=0.575)&&(Q2<=1.3))||((W>=1.3125)&&(W<1.5125)&&(Q2>=0.002)&&(Q2<=0.275))||((W>=1.5125)&&(W<=1.5875)&&(Q2>=0.002)&&(Q2<=0.225))){
 
 get_xsect_rip_fed_join(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final); 

 };
 
  

 if ((W>=1.2375)&&(W<1.3125)&&(Q2>0.002)&&(Q2<1.3)) {


get_xsect_near_threshold(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );


};


// if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.000001)&&(Q2<0.65)) {
//get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_1);
//get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_2); 
//sigma_total = 1./0.025;
//sigma_total = sigma_total*(sigma_total_1*fabs(1.5875-W)+sigma_total_2*fabs(1.6125-W));
//};

if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.000001)&&(Q2<0.65)) {
get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

 get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1); 

sigma_t_final = 1./0.025;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(1.5875-W)+sigma_t_final_1*fabs(1.6125-W));


sigma_l_final = 1./0.025;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(1.5875-W)+sigma_l_final_1*fabs(1.6125-W));



sigma_c2f_final = 1./0.025;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_1*fabs(1.5875-W)+sigma_c2f_final_2*fabs(1.6125-W));

sigma_s2f_final = 1./0.025;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_1*fabs(1.5875-W)+sigma_s2f_final_2*fabs(1.6125-W));

sigma_cf_final = 1./0.025;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_1*fabs(1.5875-W)+sigma_cf_final_2*fabs(1.6125-W));

sigma_sf_final = 1./0.025;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_1*fabs(1.5875-W)+sigma_sf_final_2*fabs(1.6125-W));

};

//sigma_t_final = 0.;
//sigma_l_final = 0.;

/*sigma_c2f_final = 0.;
sigma_s2f_final = 0.;
sigma_cf_final = 0.;
sigma_sf_final = 0.;*/


//calculating sigma_total from different sigmas, eps_l and eps_t
Float_t eps_l,eps_t,nu_g,theta_el;

nu_g = (W*W + Q2 - MP*MP)/2./MP;
theta_el = acos(1.- Q2/E_beam/(E_beam - nu_g)/2.);

eps_t = 1./(1.+ 2.*(1. + nu_g*nu_g/Q2)*tan(theta_el/2.)*tan(theta_el/2.));
eps_l = Q2*eps_t/nu_g/nu_g;
sigma_total =0.;

sigma_total = sigma_t_final;
sigma_total = sigma_total + eps_l*sigma_l_final;
sigma_total = sigma_total + eps_t*(sigma_c2f_final*cos(2.*ph_hadr) + sigma_s2f_final*sin(2.*ph_hadr));
sigma_total = sigma_total + sqrt(2.*eps_l*(eps_t+1))*(sigma_cf_final*cos(ph_hadr) + sigma_sf_final*sin(ph_hadr));


//multiply sigma_total by virtual photon flux
/*Float_t V_flux;

V_flux = alph_const/4./M_PI;
V_flux = V_flux/E_beam/E_beam/MP/MP;
V_flux = V_flux/(1.-eps_t)/Q2;
V_flux = V_flux*W*(W*W-MP*MP);
*/
//sigma_total = sigma_total*V_flux;
//--------

if ((W>=1.8125)&&(W<=2.5)&&(Q2>0.000001)&&(Q2<0.65)) {
get_xsect_scale_gol_18_25(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total);

};

// if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.000001)&&(Q2<0.65)) {
//get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_1);
//get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_2); 
//sigma_total = 1./0.025;
//sigma_total = sigma_total*(sigma_total_1*fabs(1.5875-W)+sigma_total_2*fabs(1.6125-W));
//};
 
 
 
    
    //FOR 1-PIM, 2-PIP, 3-P
//    anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr,  MPIM, MPIP, MP, P4_PIM, P4_PIP,  P4_Pfin);

   //FOR 1-P, 2-PIP, 3-PIM
     // anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MP,MPIP, MPIM, P4_Pfin,P4_PIP, P4_PIM);
   
     //FOR 1-PIP, 2-PIM, 3-P 
     //   anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MPIP,MPIM, MP, P4_PIP,P4_PIM, P4_Pfin);
	
//-------------------------------------------------------------------------------------  
    

    //FOR 1-PIM, 2-PIP, 3-P
//     rot(Q2, E_beam, P4_E_prime,P4_PIM, P4_PIP,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
   
     //FOR 1-P, 2-PIP, 3-PIM
  //   rot(Q2, E_beam, P4_E_prime,P4_Pfin, P4_PIP,  P4_PIM, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
  //FOR 1-PIP, 2-PIM, 3-P
     //  rot(Q2, E_beam, P4_E_prime,P4_PIP, P4_PIM,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
   


        //writing generated events into desired input file
//   out_file_write(i, P4_E_prime,P4_Pfin, P4_PIP,P4_PIM);
   

if (!(sigma_total > 0.)) {
if (!(sigma_total < 0.)) {
cout << sigma_total  << "  " << sigma_t_final << "  " << sigma_l_final  << " " << W << "   " << Q2 << "   " << alph_hadr*180./M_PI << "  " << ph_hadr*180./M_PI << "  "<<  th_hadr*180./M_PI   << endl;
};
};
    

//      h_empty2[qq2][ww]->Fill(Var_2,sigma_total*2.*Var_2[1]*2.*Var_2[0]*((cos((Var_2[2]-(xmax[2] - xmin[2])/thth_max/2.)*M_PI/180.))-(cos((Var_2[2]+(xmax[2] - xmin[2])/thth_max/2.)*M_PI/180.))));
  h_empty2[qq2][ww]->Fill(Var_2,sigma_total*2.*Var_2[1]*2.*Var_2[0]*(xmax[0] - xmin[0])*(xmax[1] - xmin[1])*(sin(Var_2[2]*M_PI/180.)));
  }; // end of Bycling function check  
   
  }; // enf of Q2 vs W tiangle check 
   
    }; //end of event PHI HADRON
    }; //end of event ALPHA HADRON
    };//end of event THETA HADRON
    };//end of event M23  
    };//end of event M12 
     h_empty2[qq2][ww]->Scale(1./m12m12_max/m23m23_max/thth_max/alphalph_max/phph_max); 
//       if (h_empty2[qq2][ww]->GetNbins()>0) {
//    if (qq2 == 3)h_W->Fill(W,((h_empty2[qq2][ww]->Projection(3))->Integral())/h_empty2[qq2][ww]->GetNbins());
//    if (qq2 == 7)h_W->Fill(W,((h_empty2[qq2][ww]->Projection(3))->Integral())); 
//    };



   // END with VAR SET 2




    for (Int_t m12m12=0; m12m12<m12m12_max; m12m12++) {   
    Var_2[0] = xmin[0] + (xmax[0] - xmin[0])/m12m12_max/2.+ ( xmax[0] - xmin[0])/m12m12_max*m12m12;
    for (Int_t m23m23=0; m23m23<m23m23_max; m23m23++) { 
    Var_2[1] = xmin[1] + (xmax[1] - xmin[1])/m23m23_max/2.+ ( xmax[1] - xmin[1])/m23m23_max*m23m23;
    for (Int_t thth=0; thth<thth_max; thth++) {    
     Var_2[2] = xmin[2] + (xmax[2] - xmin[2])/thth_max/2.+ ( xmax[2] - xmin[2])/thth_max*thth; 
    for (Int_t phph=0; phph<phph_max; phph++) {     
      Var_2[3] = xmin[3] + (xmax[3] - xmin[3])/phph_max/2.+ ( xmax[3] - xmin[3])/phph_max*phph;   
    for (Int_t alphalph=0; alphalph<alphalph_max; alphalph++) {   
     Var_2[4] = xmin[4] + (xmax[4] - xmin[4])/alphalph_max/2.+ ( xmax[4] - xmin[4])/alphalph_max*alphalph; 
             

phi_e =ph_e_rndm.Uniform(phi_e_min,phi_e_max);

z_EL = z_EL_rndm.Uniform(Targ_off  - Targ_len/2.,Targ_off  + Targ_len/2.);
ph_hadr = Var_2[3]*M_PI/180.;

alph_hadr =  Var_2[4]*M_PI/180.; 

    k++;

    W = W_bin[ww];

    Q2 = Q2_bin[qq2];
    
    nu=(W*W+Q2-MP*MP)/2./MP;
    if (nu <= (E_beam-E_eprime_min))  {
    
   
    E_E_prime=E_beam-nu;
    Theta_e_prime = acos(1.-Q2/E_beam/E_E_prime/2.);
    P4_E_prime.SetXYZT(E_E_prime*cos(phi_e)*sin(Theta_e_prime),E_E_prime*sin(phi_e)*sin(Theta_e_prime),E_E_prime*cos(Theta_e_prime),E_E_prime);
  //  W_tmp = W;
   // W = 1.6125;
    dummy_int = 0;
   s12 = Var_2[0]*Var_2[0];
   s23 = Var_2[1]*Var_2[1];
   inv_m12 = sqrt(s12);
   inv_m23 = sqrt(s23); 
   dummy_int++;

    // cout << W << " " << Q2 << " " << inv_m12 << "  " << inv_m23 << "  " << sigma_total << "\n";
//   } while (G_BYCKLING(inv_m12*inv_m12,inv_m23*inv_m23,W*W,M2*M2,M1*M1,M3*M3) > -5e-07);

   th_hadr = Var_2[2]*M_PI/180.;

   
        //FOR 1-PIM, 2-PIP, 3-P
//     M1 = MPIM;
//     M2 = MPIP;
//     M3 = MP; 
     //FOR 1-P, 2-PIP, 3-PIM
    M1 = MP;
    M2 = MPIP;
    M3 = MPIM;
     //FOR 1-PIP, 2-PIM, 3-P
 //    M1 = MPIP;
 //    M2 = MPIM;
 //    M3 = MP;

   if (G_BYCKLING(inv_m12*inv_m12,inv_m23*inv_m23,W*W,M2*M2,M1*M1,M3*M3) <= -5e-07) {

   
          //FOR 1-PIM, 2-PIP, 3-P
//    anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr,  MPIM, MPIP, MP, P4_PIM, P4_PIP,  P4_Pfin);

   //FOR 1-P, 2-PIP, 3-PIM
     anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MP,MPIP, MPIM, P4_Pfin,P4_PIP, P4_PIM);
     
//   cout  << "M12_2 = " << sqrt((P4_PIP+P4_Pfin)*(P4_PIP+P4_Pfin)) <<"\n";     
   
     //FOR 1-PIP, 2-PIM, 3-P 
     //   anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MPIP,MPIM, MP, P4_PIP,P4_PIM, P4_Pfin);
	
//------------------------------------------------------------------------------------- 

 
       
          //FOR 1-PIM, 2-PIP, 3-P
         M1 = MPIM;
         M2 = MPIP;
         M3 = MP; 
     //FOR 1-P, 2-PIP, 3-PIM
 //   M1 = MP;
 //   M2 = MPIP;
 //   M3 = MPIM;
     //FOR 1-PIP, 2-PIM, 3-P
 //    M1 = MPIP;
 //    M2 = MPIM;
 //    M3 = MP; 

    //FOR 1-PIM, 2-PIP, 3-P
     rot(Q2, E_beam, P4_E_prime,P4_PIM, P4_PIP,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
//      cout  << "M12_2 = " << inv_m23 <<"\n";   
          
   
   
     //FOR 1-P, 2-PIP, 3-PIM
  //   rot(Q2, E_beam, P4_E_prime,P4_Pfin, P4_PIP,  P4_PIM, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
  //FOR 1-PIP, 2-PIM, 3-P
     //  rot(Q2, E_beam, P4_E_prime,P4_PIP, P4_PIM,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
   
   	//cout  << "M12  = " << inv_m12 << " M23 = " << inv_m23 << " th_hadr = " << th_hadr <<"\n";   
   
//   h_inv_m12 ->Fill(s12,1.);
   
   s12 = inv_m12*inv_m12;
   s23 = inv_m23*inv_m23;
	
	
 
//Getting cross section in given generated (W, Q2, s12, s23, theta, alpha)-point

 if ((W>=1.4125)&&(W<=1.8125)&&(Q2>=0.65)&&(Q2<=1.3)) {
 get_xsect_ripani(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);

};







 if ((W>=1.6125)&&(W<=1.8125)&&(Q2>0.000001)&&(Q2<0.65)) {

 get_xsect_14_18_lowq2_fit(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);
};



if (((W>=1.3125)&&(W<=1.4375)&&(Q2>0.275)&&(Q2<0.575))||((W>=1.4125)&&(W<=1.4875)&&(Q2>0.275)&&(Q2<0.525))||((W>=1.4125)&&(W<=1.5125)&&(Q2>0.275)&&(Q2<0.425))||((W>=1.4125)&&(W<=1.5375)&&(Q2>0.225)&&(Q2<0.375))||((W>=1.4125)&&(W<=1.5625)&&(Q2>0.225)&&(Q2<0.275))){

 get_xsect_fedotov(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);


};

 if (((W>=1.4125)&&(W<=1.4375)&&(Q2>=0.575)&&(Q2<=0.65))||((W>=1.4375)&&(W<=1.4875)&&(Q2>=0.525)&&(Q2<=0.65))||((W>=1.4875)&&(W<=1.5125)&&(Q2>=0.425)&&(Q2<=0.65))||((W>=1.5125)&&(W<=1.5375)&&(Q2>=0.325)&&(Q2<=0.65))||((W>=1.5375)&&(W<=1.5625)&&(Q2>=0.275)&&(Q2<=0.65))||((W>=1.5625)&&(W<=1.5875)&&(Q2>=0.225)&&(Q2<=0.65))||((W>=1.3125)&&(W<1.4125)&&(Q2>=0.575)&&(Q2<=1.3))||((W>=1.3125)&&(W<1.5125)&&(Q2>=0.002)&&(Q2<=0.275))||((W>=1.5125)&&(W<=1.5875)&&(Q2>=0.002)&&(Q2<=0.225))){
 
 get_xsect_rip_fed_join(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final); 

 };
 
  

 if ((W>=1.2375)&&(W<1.3125)&&(Q2>0.002)&&(Q2<1.3)) {


get_xsect_near_threshold(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );


};


// if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.000001)&&(Q2<0.65)) {
//get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_1);
//get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_2); 
//sigma_total = 1./0.025;
//sigma_total = sigma_total*(sigma_total_1*fabs(1.5875-W)+sigma_total_2*fabs(1.6125-W));
//};

if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.000001)&&(Q2<0.65)) {
get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

 get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1); 

sigma_t_final = 1./0.025;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(1.5875-W)+sigma_t_final_1*fabs(1.6125-W));


sigma_l_final = 1./0.025;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(1.5875-W)+sigma_l_final_1*fabs(1.6125-W));



sigma_c2f_final = 1./0.025;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_1*fabs(1.5875-W)+sigma_c2f_final_2*fabs(1.6125-W));

sigma_s2f_final = 1./0.025;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_1*fabs(1.5875-W)+sigma_s2f_final_2*fabs(1.6125-W));

sigma_cf_final = 1./0.025;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_1*fabs(1.5875-W)+sigma_cf_final_2*fabs(1.6125-W));

sigma_sf_final = 1./0.025;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_1*fabs(1.5875-W)+sigma_sf_final_2*fabs(1.6125-W));

};

//sigma_t_final = 0.;
//sigma_l_final = 0.;

/*sigma_c2f_final = 0.;
sigma_s2f_final = 0.;
sigma_cf_final = 0.;
sigma_sf_final = 0.;*/


//calculating sigma_total from different sigmas, eps_l and eps_t
Float_t eps_l,eps_t,nu_g,theta_el;

nu_g = (W*W + Q2 - MP*MP)/2./MP;
theta_el = acos(1.- Q2/E_beam/(E_beam - nu_g)/2.);

eps_t = 1./(1.+ 2.*(1. + nu_g*nu_g/Q2)*tan(theta_el/2.)*tan(theta_el/2.));
eps_l = Q2*eps_t/nu_g/nu_g;
sigma_total =0.;

sigma_total = sigma_t_final;
sigma_total = sigma_total + eps_l*sigma_l_final;
sigma_total = sigma_total + eps_t*(sigma_c2f_final*cos(2.*ph_hadr) + sigma_s2f_final*sin(2.*ph_hadr));
sigma_total = sigma_total + sqrt(2.*eps_l*(eps_t+1))*(sigma_cf_final*cos(ph_hadr) + sigma_sf_final*sin(ph_hadr));


//multiply sigma_total by virtual photon flux
/*Float_t V_flux;

V_flux = alph_const/4./M_PI;
V_flux = V_flux/E_beam/E_beam/MP/MP;
V_flux = V_flux/(1.-eps_t)/Q2;
V_flux = V_flux*W*(W*W-MP*MP);
*/
//sigma_total = sigma_total*V_flux;
//--------

if ((W>=1.8125)&&(W<=2.5)&&(Q2>0.000001)&&(Q2<0.65)) {
get_xsect_scale_gol_18_25(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total);

};

// if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.000001)&&(Q2<0.65)) {
//get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_1);
//get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_2); 
//sigma_total = 1./0.025;
//sigma_total = sigma_total*(sigma_total_1*fabs(1.5875-W)+sigma_total_2*fabs(1.6125-W));
//};
 
 
 
    
    //FOR 1-PIM, 2-PIP, 3-P
//    anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr,  MPIM, MPIP, MP, P4_PIM, P4_PIP,  P4_Pfin);

   //FOR 1-P, 2-PIP, 3-PIM
     // anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MP,MPIP, MPIM, P4_Pfin,P4_PIP, P4_PIM);
   
     //FOR 1-PIP, 2-PIM, 3-P 
     //   anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MPIP,MPIM, MP, P4_PIP,P4_PIM, P4_Pfin);
	
//-------------------------------------------------------------------------------------  
    

    //FOR 1-PIM, 2-PIP, 3-P
//     rot(Q2, E_beam, P4_E_prime,P4_PIM, P4_PIP,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
   
     //FOR 1-P, 2-PIP, 3-PIM
  //   rot(Q2, E_beam, P4_E_prime,P4_Pfin, P4_PIP,  P4_PIM, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
  //FOR 1-PIP, 2-PIM, 3-P
     //  rot(Q2, E_beam, P4_E_prime,P4_PIP, P4_PIM,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
   


        //writing generated events into desired input file
//   out_file_write(i, P4_E_prime,P4_Pfin, P4_PIP,P4_PIM);
   


    

//      h_empty2[qq2][ww]->Fill(Var_2,sigma_total*2.*Var_2[1]*2.*Var_2[0]*((cos((Var_2[2]-(xmax[2] - xmin[2])/thth_max/2.)*M_PI/180.))-(cos((Var_2[2]+(xmax[2] - xmin[2])/thth_max/2.)*M_PI/180.))));

if (!(sigma_total > 0.)) {
if (!(sigma_total < 0.)) {
cout << sigma_total  << "  " << sigma_t_final << "  " << sigma_l_final  << " " << W << "  " << Q2 << "  " << alph_hadr*180./M_PI << "  " << ph_hadr*180./M_PI << "  "<<  th_hadr*180./M_PI   << endl;
};
};
 
  h_empty1[qq2][ww]->Fill(Var_2,sigma_total*2.*Var_2[1]*2.*Var_2[0]*(xmax[0] - xmin[0])*(xmax[1] - xmin[1])*(sin(Var_2[2]*M_PI/180.)));

  
//  cout << Var_2[0] << "  " <<  Var_2[1] << "  " <<  Var_2[2] << "  " <<  Var_2[3] << "  " <<  Var_2[4] << endl;

  }; // end of Bycling function check  
   
  }; // enf of Q2 vs W tiangle check 
   
    }; //end of event PHI HADRON
    }; //end of event ALPHA HADRON
    };//end of event THETA HADRON
    };//end of event M23  
    };//end of event M12 
     h_empty1[qq2][ww]->Scale(1./m12m12_max/m23m23_max/thth_max/alphalph_max/phph_max); 
//       if (h_empty2[qq2][ww]->GetNbins()>0) {
//    if (qq2 == 3)h_W->Fill(W,((h_empty2[qq2][ww]->Projection(3))->Integral())/h_empty2[qq2][ww]->GetNbins());
//    if (qq2 == 7)h_W->Fill(W,((h_empty1[qq2][ww]->Projection(3))->Integral())); 
//    };























   // END with VAR SET 1



    for (Int_t m12m12=0; m12m12<m12m12_max; m12m12++) {   
    Var_2[0] = xmin[0] + (xmax[0] - xmin[0])/m12m12_max/2.+ ( xmax[0] - xmin[0])/m12m12_max*m12m12;
    for (Int_t m23m23=0; m23m23<m23m23_max; m23m23++) { 
    Var_2[1] = xmin[1] + (xmax[1] - xmin[1])/m23m23_max/2.+ ( xmax[1] - xmin[1])/m23m23_max*m23m23;
    for (Int_t thth=0; thth<thth_max; thth++) {    
     Var_2[2] = xmin[2] + (xmax[2] - xmin[2])/thth_max/2.+ ( xmax[2] - xmin[2])/thth_max*thth; 
   for (Int_t phph=0; phph<phph_max; phph++) {     
      Var_2[3] = xmin[3] + (xmax[3] - xmin[3])/phph_max/2.+ ( xmax[3] - xmin[3])/phph_max*phph;  
    for (Int_t alphalph=0; alphalph<alphalph_max; alphalph++) {   
     Var_2[4] = xmin[4] + (xmax[4] - xmin[4])/alphalph_max/2.+ ( xmax[4] - xmin[4])/alphalph_max*alphalph; 
      

phi_e =ph_e_rndm.Uniform(phi_e_min,phi_e_max);

z_EL = z_EL_rndm.Uniform(Targ_off  - Targ_len/2.,Targ_off  + Targ_len/2.);
ph_hadr = Var_2[3]*M_PI/180.;

alph_hadr =  Var_2[4]*M_PI/180.; 

    k++;

    W = W_bin[ww];

    Q2 = Q2_bin[qq2];
    
    nu=(W*W+Q2-MP*MP)/2./MP;
    if (nu <= (E_beam-E_eprime_min))  {
    
   
    E_E_prime=E_beam-nu;
    Theta_e_prime = acos(1.-Q2/E_beam/E_E_prime/2.);
    P4_E_prime.SetXYZT(E_E_prime*cos(phi_e)*sin(Theta_e_prime),E_E_prime*sin(phi_e)*sin(Theta_e_prime),E_E_prime*cos(Theta_e_prime),E_E_prime);
  //  W_tmp = W;
   // W = 1.6125;
    dummy_int = 0;
   s12 = Var_2[1]*Var_2[1];
   s23 = Var_2[0]*Var_2[0];
   inv_m12 = sqrt(s12);
   inv_m23 = sqrt(s23); 
   dummy_int++;

    // cout << W << " " << Q2 << " " << inv_m12 << "  " << inv_m23 << "  " << sigma_total << "\n";
//   } while (G_BYCKLING(inv_m12*inv_m12,inv_m23*inv_m23,W*W,M2*M2,M1*M1,M3*M3) > -5e-07);

   th_hadr = Var_2[2]*M_PI/180.;

   
        //FOR 1-PIM, 2-PIP, 3-P
//     M1 = MPIM;
//     M2 = MPIP;
//     M3 = MP; 
     //FOR 1-P, 2-PIP, 3-PIM
//    M1 = MP;
//    M2 = MPIP;
//    M3 = MPIM;
     //FOR 1-PIP, 2-PIM, 3-P
     M1 = MPIP;
     M2 = MPIM;
     M3 = MP;

   if (G_BYCKLING(inv_m12*inv_m12,inv_m23*inv_m23,W*W,M2*M2,M1*M1,M3*M3) <= -5e-07) {

   
          //FOR 1-PIM, 2-PIP, 3-P
//    anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr,  MPIM, MPIP, MP, P4_PIM, P4_PIP,  P4_Pfin);

   //FOR 1-P, 2-PIP, 3-PIM
//     anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MP,MPIP, MPIM, P4_Pfin,P4_PIP, P4_PIM);
     
//   cout  << "M12_2 = " << sqrt((P4_PIP+P4_Pfin)*(P4_PIP+P4_Pfin)) <<"\n";     
   
     //FOR 1-PIP, 2-PIM, 3-P 
        anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MPIP,MPIM, MP, P4_PIP,P4_PIM, P4_Pfin);
	
//------------------------------------------------------------------------------------- 

 
       
          //FOR 1-PIM, 2-PIP, 3-P
         M1 = MPIM;
         M2 = MPIP;
         M3 = MP; 
     //FOR 1-P, 2-PIP, 3-PIM
 //   M1 = MP;
 //   M2 = MPIP;
 //   M3 = MPIM;
     //FOR 1-PIP, 2-PIM, 3-P
 //    M1 = MPIP;
 //    M2 = MPIM;
 //    M3 = MP; 

    //FOR 1-PIM, 2-PIP, 3-P
     rot(Q2, E_beam, P4_E_prime,P4_PIM, P4_PIP,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
//      cout  << "M12_2 = " << inv_m23 <<"\n";   
          
   
   
     //FOR 1-P, 2-PIP, 3-PIM
  //   rot(Q2, E_beam, P4_E_prime,P4_Pfin, P4_PIP,  P4_PIM, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
  //FOR 1-PIP, 2-PIM, 3-P
     //  rot(Q2, E_beam, P4_E_prime,P4_PIP, P4_PIM,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
   
   	//cout  << "M12  = " << inv_m12 << " M23 = " << inv_m23 << " th_hadr = " << th_hadr <<"\n";   
   
//   h_inv_m12 ->Fill(s12,1.);
   
   s12 = inv_m12*inv_m12;
   s23 = inv_m23*inv_m23;
	
	
 
//Getting cross section in given generated (W, Q2, s12, s23, theta, alpha)-point

 if ((W>=1.4125)&&(W<=1.8125)&&(Q2>=0.65)&&(Q2<=1.3)) {
 get_xsect_ripani(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);

};







 if ((W>=1.6125)&&(W<=1.8125)&&(Q2>0.000001)&&(Q2<0.65)) {

 get_xsect_14_18_lowq2_fit(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);
};



if (((W>=1.3125)&&(W<=1.4375)&&(Q2>0.275)&&(Q2<0.575))||((W>=1.4125)&&(W<=1.4875)&&(Q2>0.275)&&(Q2<0.525))||((W>=1.4125)&&(W<=1.5125)&&(Q2>0.275)&&(Q2<0.425))||((W>=1.4125)&&(W<=1.5375)&&(Q2>0.225)&&(Q2<0.375))||((W>=1.4125)&&(W<=1.5625)&&(Q2>0.225)&&(Q2<0.275))){

 get_xsect_fedotov(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);


};

 if (((W>=1.4125)&&(W<=1.4375)&&(Q2>=0.575)&&(Q2<=0.65))||((W>=1.4375)&&(W<=1.4875)&&(Q2>=0.525)&&(Q2<=0.65))||((W>=1.4875)&&(W<=1.5125)&&(Q2>=0.425)&&(Q2<=0.65))||((W>=1.5125)&&(W<=1.5375)&&(Q2>=0.325)&&(Q2<=0.65))||((W>=1.5375)&&(W<=1.5625)&&(Q2>=0.275)&&(Q2<=0.65))||((W>=1.5625)&&(W<=1.5875)&&(Q2>=0.225)&&(Q2<=0.65))||((W>=1.3125)&&(W<1.4125)&&(Q2>=0.575)&&(Q2<=1.3))||((W>=1.3125)&&(W<1.5125)&&(Q2>=0.002)&&(Q2<=0.275))||((W>=1.5125)&&(W<=1.5875)&&(Q2>=0.002)&&(Q2<=0.225))){
 
 get_xsect_rip_fed_join(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final); 

 };
 
  

 if ((W>=1.2375)&&(W<1.3125)&&(Q2>0.002)&&(Q2<1.3)) {


get_xsect_near_threshold(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );


};


// if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.000001)&&(Q2<0.65)) {
//get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_1);
//get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_2); 
//sigma_total = 1./0.025;
//sigma_total = sigma_total*(sigma_total_1*fabs(1.5875-W)+sigma_total_2*fabs(1.6125-W));
//};

if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.000001)&&(Q2<0.65)) {
get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

 get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1); 

sigma_t_final = 1./0.025;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(1.5875-W)+sigma_t_final_1*fabs(1.6125-W));


sigma_l_final = 1./0.025;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(1.5875-W)+sigma_l_final_1*fabs(1.6125-W));



sigma_c2f_final = 1./0.025;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_1*fabs(1.5875-W)+sigma_c2f_final_2*fabs(1.6125-W));

sigma_s2f_final = 1./0.025;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_1*fabs(1.5875-W)+sigma_s2f_final_2*fabs(1.6125-W));

sigma_cf_final = 1./0.025;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_1*fabs(1.5875-W)+sigma_cf_final_2*fabs(1.6125-W));

sigma_sf_final = 1./0.025;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_1*fabs(1.5875-W)+sigma_sf_final_2*fabs(1.6125-W));

};

//sigma_t_final = 0.;
//sigma_l_final = 0.;

/*sigma_c2f_final = 0.;
sigma_s2f_final = 0.;
sigma_cf_final = 0.;
sigma_sf_final = 0.;*/


//calculating sigma_total from different sigmas, eps_l and eps_t
Float_t eps_l,eps_t,nu_g,theta_el;

nu_g = (W*W + Q2 - MP*MP)/2./MP;
theta_el = acos(1.- Q2/E_beam/(E_beam - nu_g)/2.);

eps_t = 1./(1.+ 2.*(1. + nu_g*nu_g/Q2)*tan(theta_el/2.)*tan(theta_el/2.));
eps_l = Q2*eps_t/nu_g/nu_g;
sigma_total =0.;

sigma_total = sigma_t_final;
sigma_total = sigma_total + eps_l*sigma_l_final;
sigma_total = sigma_total + eps_t*(sigma_c2f_final*cos(2.*ph_hadr) + sigma_s2f_final*sin(2.*ph_hadr));
sigma_total = sigma_total + sqrt(2.*eps_l*(eps_t+1))*(sigma_cf_final*cos(ph_hadr) + sigma_sf_final*sin(ph_hadr));


//multiply sigma_total by virtual photon flux
/*Float_t V_flux;

V_flux = alph_const/4./M_PI;
V_flux = V_flux/E_beam/E_beam/MP/MP;
V_flux = V_flux/(1.-eps_t)/Q2;
V_flux = V_flux*W*(W*W-MP*MP);
*/
//sigma_total = sigma_total*V_flux;
//--------

if ((W>=1.8125)&&(W<=2.5)&&(Q2>0.000001)&&(Q2<0.65)) {
get_xsect_scale_gol_18_25(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total);

};

// if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.000001)&&(Q2<0.65)) {
//get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_1);
//get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_total_2); 
//sigma_total = 1./0.025;
//sigma_total = sigma_total*(sigma_total_1*fabs(1.5875-W)+sigma_total_2*fabs(1.6125-W));
//};
 
 
 
    
    //FOR 1-PIM, 2-PIP, 3-P
//    anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr,  MPIM, MPIP, MP, P4_PIM, P4_PIP,  P4_Pfin);

   //FOR 1-P, 2-PIP, 3-PIM
     // anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MP,MPIP, MPIM, P4_Pfin,P4_PIP, P4_PIM);
   
     //FOR 1-PIP, 2-PIM, 3-P 
     //   anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MPIP,MPIM, MP, P4_PIP,P4_PIM, P4_Pfin);
	
//-------------------------------------------------------------------------------------  
    

    //FOR 1-PIM, 2-PIP, 3-P
//     rot(Q2, E_beam, P4_E_prime,P4_PIM, P4_PIP,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
   
     //FOR 1-P, 2-PIP, 3-PIM
  //   rot(Q2, E_beam, P4_E_prime,P4_Pfin, P4_PIP,  P4_PIM, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
  //FOR 1-PIP, 2-PIM, 3-P
     //  rot(Q2, E_beam, P4_E_prime,P4_PIP, P4_PIM,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
   


        //writing generated events into desired input file
//   out_file_write(i, P4_E_prime,P4_Pfin, P4_PIP,P4_PIM);
   


    

//      h_empty2[qq2][ww]->Fill(Var_2,sigma_total*2.*Var_2[1]*2.*Var_2[0]*((cos((Var_2[2]-(xmax[2] - xmin[2])/thth_max/2.)*M_PI/180.))-(cos((Var_2[2]+(xmax[2] - xmin[2])/thth_max/2.)*M_PI/180.))));

 
  h_empty3[qq2][ww]->Fill(Var_2,sigma_total*2.*Var_2[1]*2.*Var_2[0]*(xmax[0] - xmin[0])*(xmax[1] - xmin[1])*(sin(Var_2[2]*M_PI/180.)));
  }; // end of Bycling function check  
   
  }; // enf of Q2 vs W tiangle check 
   
    }; //end of event PHI HADRON
    }; //end of event ALPHA HADRON
    };//end of event THETA HADRON
    };//end of event M23  
    };//end of event M12 
     h_empty3[qq2][ww]->Scale(1./m12m12_max/m23m23_max/thth_max/alphalph_max/phph_max); 
//       if (h_empty2[qq2][ww]->GetNbins()>0) {
//    if (qq2 == 3)h_W->Fill(W,((h_empty2[qq2][ww]->Projection(3))->Integral())/h_empty2[qq2][ww]->GetNbins());
    if (qq2 == 7)h_W->Fill(W,((h_empty3[qq2][ww]->Projection(3))->Integral())); 
//    };
















    };//end of event W 
    };//end of event Q2 PHI HADRON
    
  
    
    //closing output file
    out_file_close();
    
    
     TFile *file1 = TFile::Open("empty_cells_new_bin.root","RECREATE");



file1->cd();
for (Int_t qq2=0; qq2<22;qq2++) { 
for (Int_t ww=0; ww<22;ww++) {
qqq.str("");
qqq << "q2_" << Q2_bin[qq2] << "/w_" << W_bin[ww];
file1->mkdir(qqq.str().c_str());
file1->cd(qqq.str().c_str());
h_empty1[qq2][ww]->Write();
h_empty2[qq2][ww]->Write();
h_empty3[qq2][ww]->Write();
qqq.str("");
};
};


file1->cd();


h_W->SetMinimum(0.);
h_Q2->SetMinimum(0.);
h_Q2vsW->SetMinimum(0.); 




h_W->Write("", TObject::kOverwrite);
h_Q2->Write("", TObject::kOverwrite);
h_Q2vsW->Write("", TObject::kOverwrite);

//for (Short_t ii=0; ii<15; ii++) {
//h_odn_q2_dep[ii]->Write("", TObject::kOverwrite); 
//};
file1->Write(); 
    
    
    

 
  
    TApplication *theApp = new TApplication("App", &argc, argv);
//    MyMainFrame jopa;
//    jopa.MainFrame(flag,E0,number_of_files,file,outfile_inp);
    
 //   TCanvas *cresolv = new TCanvas("cresolv", "Resolution", 500, 800);
//    TGMainFrame *fMain = new TGMainFrame(gClient->GetRoot(),200,200);
//    fMain->SetWindowName("Simple Example");
 //   TRootEmbeddedCanvas *fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,200,200); 
//    fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY,10,10,10,1));
   
 

   
//   theApp.Run();
  


  delete [] file;
  file = NULL;
  
  return 0;
}
