ifndef F90C
F90C=gfortran
endif
ifndef CC
CC=gcc
endif
RTSUPP=w2f__types OAD_active OAD_cp OAD_tape OAD_rev   
driver: $(addsuffix .o, $(RTSUPP) iaddr) lorenz63_passive.o numCore.pre.xb.x2w.w2f.post.o driver.o
	${F90C} -o $@ $^
numCore.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90, $(RTSUPP)) iaddr.c : toolChain 
toolChain : numCore.f90
	openad -c -m rj $<

numCore.f90: Lorenz63.f90 head.prepped.f90
	cat $^ > $@

lorenz63_passive.f90: Lorenz63.f90
	cat $^ | sed 's/Lorenz63/lorenz63_passive/' > $@
  

ad_inline.f:toolChain

%.o : %.f90
	${F90C} -o $@ -c $< 

%.o : %.f
	${F90C} ${F90FLAGS} -g -O -o $@ -c $< 

%.o : %.c
	${CC} -o $@ -c $< 


clean: 
	rm -f ad_template* ad_inline.f OAD_* w2f__*  iaddr* 
	rm -f head.prepped.pre.* *.B *.xaif *.o *.mod driver driverE *~ 
.PHONY: clean toolChain
# the following include has explicit rules that could replace the openad script
include MakeExplRules.inc
