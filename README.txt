$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$  Pulse-pileup code
$
$  BUILD / INSTALLATION INSTRUCTIONS
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

Requires g++

There are several libraries that will be needed: 
CFITSIO
GSL (Gnu Scientific Library)
HDF5 (see below)

Installing HDF5:
--------------------------------
* Get the source code
The source code can be obtained via:
http://www.hdfgroup.org/HDF5/release/obtain5.html

I clicked the ALL PLATFORMS link about halfway down the page. From the next page I downloaded the tar.gz distribution (under "has UNIX line endings"). 


* Configure HDF5.
Follow the installation instructions, but MAKE SURE TO USE '--enable-fortran' and '--enable-cxx' when configuring HDF5.
The prefix is just an example I used and could be anywhere on your system.

./configure --prefix=/usr/local/hdf5 --enable-fortran --enable-cxx CFLAGS="-arch x86_64" CXXFLAGS="-arch x86_64"

Then make, make install.

OF course, the library will have be linkable when the PPU code is run, so I suggest adding lines like this to your profile

export HDF5=/usr/local/hdf5
export PATH=$PATH:$HDF5/bin
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$HDF5/lib


--------------------------------




BUILDING PPU CODE
--------------------------------

cd to the PPU code directory.  There are several source files here, and more in a subdirectory called src2/. 
The location of the libraries needs to be given before configuring (unless they are in a default path)

export GSL=/path/to/GSL      (contains, for example, include/gsl/gsl_rng.h   and   lib/libgsl.*, among others)
export HDF5=/path/to/HDF5    (contains include/hdf5.h   and   lib/libhdf5.*.  From the example above, I have HDF5=/usr/local/hdf5  )
export HEADAS=/path/to/cfitsio

If any of these were installed into default paths, those steps aren't necessary.

The final configuration option to consider is the environment variable HDF5_LDFLAGS. Since HDF5 has its own library dependencies, it has to be compiled with several flags.  This variables provides that option, but you probably do not need to do anything.  The defaults are set to "-lz -lm -lhdf5 -lhdf5_cpp".  If however you compiled HDF5 without ZLIB, the "-lz" would need to be removed. 

Now, configure and build the PPU code. Use CXXFLAGS to pass the other g++ compilation options needed so it can later be called from fortran:

./configure --prefix=<prefix>  CXXFLAGS="-fPIC"
make

Doing `make install` installs two programs into <prefix>/bin, one for computing distortion of PHA file, the other for remaking the HDF5 pileup response.  Each one has a '--help' option that explains them.


OBJECT FILES (.o) LEFT OVER FROM BUILD:
--------------------------------
In the main directory will be several .o files from the build.  To figure out your GFORTRAN statement, I would just look at the last command of `make` (before `make install`). It should look something like:

g++  -g -O2   -o ppu_MakeH5Response DBStringUtilities.o PHAStructures.o PHA_IO.o spoccExeUtilities.o ADC.o AnalyticalPPU.o GapStatistics.o MeasuredEnergies.o PrRspsHDF5.o SpectrumHDF5.o ppu_static_interface.o PPU_MakeH5Response.o -L/usr/local/heasoft-6.12/x86_64-apple-darwin11.3.0/lib  -L/usr/local/hdf5/lib -lgsl -lcfitsio -lz -lm -lhdf5 -lhdf5_cpp


Only these object files are needed for the PPU interface:
ADC.o AnalyticalPPU.o GapStatistics.o MeasuredEnergies.o PrRspsHDF5.o SpectrumHDF5.o ppu_static_interface.o

And they are linked against GSL and HDF5 (and dependencies), not fitsio.  For example linking with gfotran might look like

gfortran -L/usr/local/hdf5/lib -lgsl -lz -lm -lhdf5 -lhdf5_cpp \
	ADC.o AnalyticalPPU.o GapStatistics.o MeasuredEnergies.o PrRspsHDF5.o SpectrumHDF5.o ppu_static_interface.o \
	…your objects and flags…





IDL Interface
---------------

The file Make_PPU_IDL.sh compiles the idl-c++ interface.  Open it for instructions.

Then put the idl_code/ subdirectory in your IDL_PATH. This aren't "installed" by 'make install' and can be put anywhere.

Edit the file link_ppu_library. At the top is a line to specify the library path.  Replace it with the one build by Make_PPU_IDL.sh.

Upon calling IDL, THE FIRST THING YOU MUST DO is run link_ppu_library.  If it's not run before any routine containing the ppu functions, they will not compile properly.  IDL is brilliant.


















