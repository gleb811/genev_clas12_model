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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <TLorentzVector.h>
#include "global.h"
 using namespace std;

int out_file_open() {

//Create BOS output if needed
 if ((flag_bos == 1)||(flag_bos == 2)){
char mess[256];
      initbos();
      bankList( &bcs_, "E=", "HEADPARTMCTKMCVXMCEV" );
     sprintf( mess, "OPEN BOSOUTPUT UNIT=1 FILE=\"%s\" WRITE STATUS=NEW RECL=3600", out_bos_file.c_str() );
     fparm_c( mess );
     
     
     };//end of BOS output flag check
};
