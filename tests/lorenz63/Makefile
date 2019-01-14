ifndef F90C
F90C=gfortran
endif
ifndef CC
CC=gcc
endif
F90CC=gfortran-debug
RTSUPP=w2f__types OAD_active_custom OAD_cp OAD_tape OAD_rev   
#RTSUPP=w2f__types OAD_active OAD_cp OAD_tape OAD_rev   
#RTSUPP=w2f__types OAD_active  

#Tests

solver_h_test: $(addsuffix .o, $(RTSUPP)) \
lorenz63_passive.o numCoreh.pre.xb.x2w.w2f.post.o driver_homogeneous.o
	${F90C} -o $@ $^

solver_i_test: $(addsuffix .o, $(RTSUPP)) \ \
lorenz63_passive.o numCorei.pre.xb.x2w.w2f.post.o driver_inhomogeneous.o
	${F90C} -o $@ $^

driver: $(addsuffix .o, $(RTSUPP)) \  \
lorenz63_passive.o numCore.pre.xb.x2w.w2f.post.o driver.o
	${F90C} -o $@ $^


#Interface with nilss

solver_ht: $(addsuffix .o, $(RTSUPP)) \
lorenz63_passive.o numCoreh.pre.xb.x2w.w2f.post.o driver_ht_nilss.o
	${F90C} -o $@ $^

#solver_iht: $(addsuffix .o, $(RTSUPP)) \ \
lorenz63_passive.o numCorei.pre.xb.x2w.w2f.post.o driver_iht_nilss.o
#	${F90C} -o $@ $^


solver_primal: \
lorenz63_passive.o driver_primal.o 
	${F90C} -o $@ $^ 

test: $(addsuffix .o, $(RTSUPP)) lorenz63_passive.o  numCoreh.pre.xb.x2w.w2f.post.o \
driver_ht_nilss_c.o driver_fortran_ht.o 
	${F90C} -o $@ $^ 

driver_fortran_ht.o: driver_fortran_ht.f90 
	${F90C} -ffree-form -c $^

for-sub.o: lorenz63_passive.f90 for-sub.f90  
	${F90C} -ffree-form -c $^

driver_primal_c.o: driver_primal_c.c
	${CC} -c -fPIC $^ 

driver_fortran_primal.o : driver_fortran_primal.f90
	${F90C} -c -fPIC $^

lorenz63_passive.o : lorenz63_passive.f90
	${F90C} -c -fPIC $^


test2: lorenz63_passive.f90 \
driver_fortran_primal.f90 
	f2py -c -m libprimal $^ 



test1: driver_primal_c.o lorenz63_passive.o \
driver_fortran_primal.o 
	${F90C} -o $@ $^
	${CC} -shared -lgfortran -Wl,-soname,libfoo.so -o libfoo.so driver_primal_c.o lorenz63_passive.o driver_fortran_primal.o  
	#ld -shared -o libfoo.so $^
numCore.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90, $(RTSUPP)) iaddr.c : toolChain 

numCoreh.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90, $(RTSUPP)) iaddr.c : toolChainh 


#numCorei.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90,  $(RTSUPP)) iaddr.c : toolChaini 
	

#toolChain : numCore.f90
#	openad -c -m rj $<
#	#rsync -av ${OPENADROOT}/runTimeSupport/vector/OAD_active.f90 ./

#numCore.f90: Lorenz63.f90 head.prepped.f90
#	cat $^ > $@

toolChainh : numCoreh.f90
	openad -c -m rj $<
	rsync -av ${OPENADROOT}/runTimeSupport/vector/OAD_active.f90 ./
	cat OAD_active.f90 | sed 's/100/3/' > OAD_active_custom.f90		

numCoreh.f90: Lorenz63.f90 head.prepped.homo.f90
	cat $^ > $@

#toolChaini : numCorei.f90
#	openad -c -m rj $<
#	rsync -av ${OPENADROOT}/runTimeSupport/vector/OAD_active.f90 ./
#	cat OAD_active.f90 | sed 's/100/3/' > OAD_active_custom.f90	
#numCorei.f90: Lorenz63.f90 head.prepped.inhomo.f90
#	cat $^ > $@




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
	rm -f numCore* *.so 
.PHONY: clean toolChain
# the following include has explicit rules that could replace the openad script
include MakeExplRules.inc
