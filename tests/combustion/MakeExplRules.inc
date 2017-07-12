# explicit make rules
# preprocess
head.prepped.pre.f90: head.prepped.f90
	${OPENADFORTTK_BASE}/tools/SourceProcessing/preProcess.py -m r --inputFormat=free -o $@ $<
# fortran -> whirl
head.prepped.pre.B: head.prepped.pre.f90 
	${OPEN64ROOT}/crayf90/sgi/mfef90 -z -F -N132 $<
# whirl -> xaif
head.prepped.pre.xaif : head.prepped.pre.B 
	${OPENADFORTTKROOT}/bin/whirl2xaif -n -o $@ $<
# xaif -> xaif'
head.prepped.pre.xb.xaif : head.prepped.pre.xaif
	${XAIFBOOSTERROOT}/xaifBooster/algorithms/BasicBlockPreaccumulationReverse/driver/oadDriver \
	-c ${XAIFSCHEMAROOT}/schema/examples/inlinable_intrinsics.xaif \
	-s ${XAIFSCHEMAROOT}/schema -i $< -o $@ 
# xaif' -> whirl'
head.prepped.pre.xb.x2w.B : head.prepped.pre.xb.xaif  
	${OPENADFORTTKROOT}/bin/xaif2whirl head.prepped.pre.B $<
# whirl' -> fortran'
head.prepped.pre.xb.x2w.w2f.f: head.prepped.pre.xb.x2w.B
	${OPEN64ROOT}/whirl2f/whirl2f -openad $<
# postprocess 
head.prepped.pre.xb.x2w.w2f.post.E.f90: head.prepped.pre.xb.x2w.w2f.f
	${OPENADFORTTK_BASE}/tools/SourceProcessing/postProcess.py -m r \
	--infoUnitFile=w2f__types.E.f90 \
	-i ${OPENADROOT}/runTimeSupport/simple/ad_inline.f		\
	-t ${OPENADROOT}/runTimeSupport/simple/ad_template.joint.f      \
	--outputFormat=free -o $@ $<
# we add the .E extension  here to distinguish the targets from the 
# rules using the openad script
driverE: $(addsuffix .E.o, $(RTSUPP)) driver.o head.prepped.pre.xb.x2w.w2f.post.E.o
	${F90C} -o $@ $^
w2f__types.E.f90: ${OPENADROOT}/runTimeSupport/all/w2f__types.f90
	cp -f $< $@
%.E.f90: ${OPENADROOT}/runTimeSupport/scalar/%.f90
	cp -f $< $@
%.E.f90: ${OPENADROOT}/runTimeSupport/simple/%.f90
	cp -f $< $@
ad_template.f: ${OPENADROOT}/runTimeSupport/simple/ad_template.joint.f
	cp -f $< $@

