ifndef F90C
F90C=gfortran
endif
ifndef CC
CC=gcc
endif
RTSUPP=w2f__types OAD_active_custom OAD_cp OAD_tape OAD_rev   
#RTSUPP=w2f__types OAD_active OAD_cp OAD_tape OAD_rev   
#RTSUPP=w2f__types OAD_active  


#Interface with nilss
solver_ht: $(addsuffix .o, $(RTSUPP)) \
cheb.o combustion_qiqi_passive.o numCoreh.pre.xb.x2w.w2f.post.o driver_ht_nilss.o
	${F90C} -o $@ $^

solver_iht: $(addsuffix .o, $(RTSUPP)) \ \
cheb.o combustion_qiqi_passive.o numCorei.pre.xb.x2w.w2f.post.o driver_iht_nilss.o
	${F90C} -o $@ $^


solver_primal: \
cheb.o combustion_qiqi_passive.o driver_primal.o 
	${F90C} -o $@ $^ 

solver_attractor: \
cheb.o combustion_qiqi_passive.o driver_attractor.o 
	${F90C} -o $@ $^ 

solver_phaseportrait: \
cheb.o combustion_qiqi_passive.o phase_portrait.o 
	${F90C} -o $@ $^ 

testcheb: \
cheb.o test_cheb.o
	${F90C} -o $@ $^ 


numCoreh.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90, $(RTSUPP)) iaddr.c : toolChainh 

numCorei.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90,  $(RTSUPP)) iaddr.c : toolChaini 
	


toolChainh : numCoreh.f90
	openad -c -m rj $<
	rsync -av ${OPENADROOT}/runTimeSupport/vector/OAD_active.f90 ./
	cat OAD_active.f90 | sed 's/100/33/' > OAD_active_custom.f90		

numCoreh.f90: combustion_qiqi.f90 head.prepped.homo.f90
	cat $^ > $@

toolChaini : numCorei.f90
	openad -c -m rj $<
	rsync -av ${OPENADROOT}/runTimeSupport/vector/OAD_active.f90 ./
	cat OAD_active.f90 | sed 's/100/33/' > OAD_active_custom.f90	

numCorei.f90: cheb.f90 combustion_qiqi.f90 head.prepped.inhomo.f90
	cat $^ > $@


combustion_qiqi_passive.f90: combustion_qiqi.f90
	cat $^ | sed 's/combustion_qiqi/combustion_qiqi_passive/' > $@
  

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
	rm -f numCore* *.mod*
	rm -f combustion_qiqi_passive.f90 
.PHONY: clean toolChain
# the following include has explicit rules that could replace the openad script
include MakeExplRules.inc
