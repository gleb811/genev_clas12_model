
CXX           = g++ -O2 -Wno-write-strings -Wno-pragmas -I$(CLAS_BUILD)/packages/include
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINCLUDE  := -I$(shell root-config --incdir)



all:
	make Dictationarys
	make genev_clas12_2

genev_clas12_2:  main_prog.$(ObjSuf) global.$(ObjSuf) inp_file_read.$(ObjSuf) read_xsect_files.$(ObjSuf)  read_fit_param_files.$(ObjSuf)  out_file_open.$(ObjSuf) out_file_write.$(ObjSuf)  out_file_close.$(ObjSuf) get_xsect_ripani.$(ObjSuf) get_xsect_fedotov.$(ObjSuf)  get_xsect_rip_fed_join.$(ObjSuf) get_xsect_near_threshold.$(ObjSuf) get_xsect_golovach.$(ObjSuf) get_xsect_14_18_lowq2_fit.$(ObjSuf)  get_xsect_scale_gol_18_25.$(ObjSuf) anti_rot.$(ObjSuf) rot.$(ObjSuf) interpol.$(ObjSuf) interpol_golovach.$(ObjSuf) interpol_fedotov.$(ObjSuf) interpol_fedotov_thresh.$(ObjSuf)
	$(CXX) -g -o $@ $^ -L/usr/lib64 -L/home/fedotov/GLEB -L/u/home/gleb/lib/LinuxRHFC8 -lpid -ltag -llac -lseb -lst -lclasutil -lsc -lc_bos_io -ldc -lec -lcc -ltrk -ldc -lc_bos_io -lsc -lmapmanager -lfputil -lfpack -lrecutl -lonline_dummy -lc_cern -lclasutil -lbos -lfpack -lbankdefs -L/u/home/gleb/cern/2005/lib -lpacklib -lkernlib -lnsl -lgfortran -lmathlib -lpacklib -lkernlib -lpawlib $(ROOTGLIBS)  -lEG 

	
Dictationarys:	
#	rootcint -f MyMainFrameDict.cxx -c -I`root-config --incdir` MyMainFrame.h 
	
#	rootcint -f $(TARGETCINT) -c -I$(ROOTSYS)/include $(TARGET).h

%.$(ObjSuf): %.$(SrcSuf)
	$(CXX) -g -c $(ROOTINCLUDE) -c $<

clean:
	rm -f *.o
	rm -f *Dict.*
	rm -f G__*
	rm genev_clas12_2
	

