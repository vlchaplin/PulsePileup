


AUTOMAKE_OPTIONS = foreign

srcbase=src2/

defaultH5=$(top_srcdir)/PPUProb_Volts.h5

EXTRA_DIST = README.txt PPUExample.cpp idl_code/ppu_example.pro idl_code/ppu_fit_example.pro idl_code/link_ppu_library.pro idl_code/drmdata__define.pro ppu_IDL_interface.cpp TestingData/

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} @fits_inc@ @gsl_inc@ @h5_inc@

AM_CPPFLAGS= 	-Wall \
				-I$(top_srcdir) \
				-I$(top_srcdir)/$(srcbase) \
                -DDEFAULT_H5_PPU="\"$(defaultH5)\"" \
				${ACLOCAL_AMFLAGS}
                
basesrcs = \
	$(srcbase)/DBStringUtilities.cpp \
	$(srcbase)/DBStringUtilities.hh \
	$(srcbase)/DynMatrix.h \
	$(srcbase)/EdgeSet.hh \
	$(srcbase)/PHAElements.hh \
	$(srcbase)/PHAStructures.cpp \
	$(srcbase)/PHAStructures.hh \
	$(srcbase)/PHA_IO.cpp \
	$(srcbase)/PHA_IO.hh \
	$(srcbase)/VectorSpan.h \
	$(srcbase)/spoccExeUtilities.cpp \
	$(srcbase)/spoccExeUtilities.h
    
ppusources = \
    ADC.h \
    AnalyticalPPU.h \
    ChannelEnergy.h \
    GapStatistics.h \
    MeasuredEnergies.h \
    PrRspsHDF5.h \
    PulseShape.h \
    SpectralShapes.h \
    SpectrumHDF5.h \
    fits_static_interface.h \
    ppu_static_interface.h \
    ADC.cpp \
    AnalyticalPPU.cpp \
    GapStatistics.cpp \
    MeasuredEnergies.cpp \
    PrRspsHDF5.cpp \
    SpectrumHDF5.cpp \
    ppu_static_interface.cpp
    
    
bin_PROGRAMS = ppu_pha_spec ppu_MakeH5Response

ppu_pha_spec_LDADD = @fits_lib@ @gsl_lib@ @h5_lib@ -lgsl -lcfitsio @HDF5_LDFLAGS@ 
ppu_MakeH5Response_LDADD = @fits_lib@ @gsl_lib@ @h5_lib@ -lgsl -lcfitsio @HDF5_LDFLAGS@ 

ppu_pha_spec_SOURCES = $(basesrcs) $(ppusources) ppu_pha_spec.cpp
ppu_MakeH5Response_SOURCES = $(basesrcs) $(ppusources) PPU_MakeH5Response.cpp    