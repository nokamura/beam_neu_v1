#F77      = g77
F77      = gfortran
FFLAGS   = -O -ffixed-line-length-200 -ffixed-form
LIBDIR   = -L./lib
LIBS_MYLIB = -lmylib -lminuit

LIBS          = $(LIBS_MYLIB)

OBJS_main       = main.o hfunc1D.o minfunc.o smearing.o gran.o funcs.o\
		  MakeHisto1D.o hsimpson1D.o chi2.o \
		  hmffrd.o bsffrd.o comsep.o msgout.o tolow.o leng.o \
		  interpolate.o CEigensystem.o zheevj3.o sqrabs.o \
		  func_smear.o ret_bins.o convert_th.o \
		  get_nudist.o get_xsecfrac.o fErec_dist.o \
		  read_histo.o get_1pi0dist.o fErec_1pi0dist.o fpi0mom.o \
		  fpi0mom_old.o fpi0_Erec.o fpi0_Erec_old.o write_input.o \
		  init_minuit.o write_minuit_results.o write_results_mode1.o \
		  bining.o get_event.o get_event_nu.o get_event_pi0.o\
		  get_event_allmode.o make_basic_data2.o init_minfunc.o \
		  get_chisq_stat.o display_results.o prepare_events.o \
		  get_event_all.o get_recevent.o \
		  classify_recevent.f apply_Ereccut.o write_chi2event.o \
		  write_params.f write_params_unit.o set_params.o \
		  write_data_oscsum.o write_data.o make_prob_data.o write_pulls.o

OBJS_smeardist_mono  = smeardist_mono.o get_nudist_mode.o MakeHisto1D.o \
		       hsimpson1D.o fErec_dist.o get_xsecfrac.o read_histo.o

OBJS_flux       = interpolate_flux.o make_flux.o

OBJS_test       = test.o apply_Ereccut.o bining.o

.f.o:
	$(F77) $(FFLAGS) $(INCLUDES) -c $<

main: $(OBJS_main)
	 $(F77) $(FFLAGS) $(LIBDIR) $(OBJS_main) $(LIBS) -o $@

smeardist_mono: $(OBJS_smeardist_mono)
	 $(F77) $(FFLAGS) $(LIBDIR) $(OBJS_smeardist_mono) $(LIBS) -o $@

flux: $(OBJS_flux)
	 $(F77) $(FFLAGS) $(LIBDIR) $(OBJS_flux) $(LIBS) -o $@

test: $(OBJS_test)
	 $(F77) $(FFLAGS) $(LIBDIR) $(OBJS_test) $(LIBS) -o $@

clean: 
	@rm *.o *~ *# fort* main test

clean_all:
	@rm *.o *~ *# fort* *.dat *.eps main test