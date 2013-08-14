
bin_PROGRAMS=chenopis\
	hcat \
	unhcat
lib_LTLIBRARIES=libchenopis.la

libchenopis_la_SOURCES=\
	swmod1.F \
	swmod2.F \
	m_constants.F90 \
	m_fileio.F90 \
	serv_xnl4v5.F90 \
	mod_xnl4v5.F90 \
	swancom1.F \
	swancom2.F \
	swancom3.F \
	swancom4.F \
	swancom5.F \
	SwanGriddata.F90 \
	SwanGridobjects.F90 \
	SwanCompdata.F90 \
	swanpre1.F \
	swanpre2.F \
	swanout1.F \
	swanout2.F \
	swanser.F \
	swanparll.F \
	SwanReadGrid.F90 \
	SwanReadADCGrid.F90 \
	SwanReadTriangleGrid.F90 \
	SwanReadEasymeshGrid.F90 \
	SwanInitCompGrid.F90 \
	SwanCheckGrid.F90 \
	SwanCreateEdges.F90 \
	SwanGridTopology.F90 \
	SwanGridVert.F90 \
	SwanGridCell.F90 \
	SwanGridFace.F90 \
	SwanPrintGridInfo.F90 \
	SwanFindPoint.F90 \
	SwanPointinMesh.F90 \
	SwanBpntlist.F90 \
	SwanPrepComp.F90 \
	SwanVertlist.F90 \
	SwanCompUnstruc.F90 \
	SwanDispParm.F90 \
	SwanPropvelX.F90 \
	SwanSweepSel.F90 \
	SwanPropvelS.F90 \
	SwanTranspAc.F90 \
	SwanTranspX.F90 \
	SwanDiffPar.F90 \
	SwanGSECorr.F90 \
	SwanInterpolatePoint.F90 \
	SwanInterpolateAc.F90 \
	SwanInterpolateOutput.F90 \
	SwanConvAccur.F90 \
	SwanConvStopc.F90 \
	SwanThreadBounds.F90 \
	SwanFindObstacles.F90 \
	SwanCrossObstacle.F90 \
	SwanComputeForce.F90 \
	SwanIntgratSpc.F90 \
	SwanBndStruc.F90 \
	SwanReadfort18.F90 \
	SwanPunCollect.F90 \
	SwanSumOverNodes.F90 \
	SwanMinOverNodes.F90 \
	SwanMaxOverNodes.F90 \
	ocpids.F \
	ocpcre.F \
	ocpmix.F \
	swanapi.F

if HAVE_BMI
	libchenopis_la_SOURCES+=SwanBMI.F90
endif

BMIINCLUDES=\
	bmi_get_nd.inc\
	bmi_nvar.inc\
	bmi_set_nd.inc\
	bmi_type.inc\
	bmi_var_name.inc\
	bmi_module.inc\
	bmi_rank.inc\
	bmi_shape.inc\
	bmi_unit.inc

SwanBMI.F90: $(BMIINCLUDES)

$(BMIINCLUDES): variables.json
	python $(top_srcdir)/scripts/codegen.py

libchenopis_la_DIST=variables.json

# Add the library
chenopis_LDADD=\
	libchenopis.la

chenopis_SOURCES=\
	swanmain.F



hcat_SOURCES=swanhcat.F
hcat_DIST=hcat.nml
unhcat_SOURCES=HottifySWAN.F90

clean-local:
	rm -f *.mod





