ifndef F90C
F90C=gfortran
endif
ifndef CC
CC=gcc
endif
RTSUPP=w2f__types OAD_active OAD_cp OAD_tape OAD_rev   

solver_h: $(addsuffix .o, $(RTSUPP) iaddr) lorenz63_passive.o numCoreh.pre.xb.x2w.w2f.post.o driver_homogeneous.o
	${F90C} -o $@ $^

solver_i: $(addsuffix .o, $(RTSUPP) iaddr) lorenz63_passive.o numCorei.pre.xb.x2w.w2f.post.o driver_inhomogeneous.o
	${F90C} -o $@ $^


driver: $(addsuffix .o, $(RTSUPP) iaddr) lorenz63_passive.o numCore.pre.xb.x2w.w2f.post.o driver.o
	${F90C} -o $@ $^

numCore.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90, $(RTSUPP)) iaddr.c : toolChain 

numCoreh.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90, $(RTSUPP)) iaddr.c : toolChainh 

numCorei.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90, $(RTSUPP)) iaddr.c : toolChaini 

toolChain : numCore.f90
	openad -c -m rj $<

numCore.f90: Lorenz63.f90 head.prepped.f90
	cat $^ > $@

toolChainh : numCoreh.f90
	openad -c -m rj $<

numCoreh.f90: Lorenz63.f90 head.prepped.homo.f90
	cat $^ > $@

toolChaini : numCorei.f90
	openad -c -m rj $<

numCorei.f90: Lorenz63.f90 head.prepped.inhomo.f90
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
