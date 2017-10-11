molBLOCKS v0.1 -- 01/23/14
Copyright (c) 2014  Dario Ghersi.

Thank you for downloading molBLOCKS. molBLOCKS is
described in:

molBLOCKS: decomposing small molecule sets and uncovering enriched fragments
Dario Ghersi and Mona Singh
Bioinformatics, 2014 Jul 15;30(14):2081-3. doi: 10.1093/bioinformatics/btu173

If you use this program, please cite the paper.

-----------------------------------------------------------------------------

INSTALLATION -- SOURCE CODE

molBLOCKS is available as source code and as a ready-to-go Virtual Machine.

molBLOCKS has one external dependency, the Open Babel library.
Open Babel needs to be compiled and installed first.

In order to compile and install Open Babel, make sure your system is
capable of building C/C++ programs. Compilers can be readily obtained
for both Mac OS X (Xcode development environment) and Linux.
Open Babel also requires CMake, available for download here
(http://www.cmake.org/cmake/resources/software.html).

Then, please type the following:

# tar xzvf openbabel-2.3.2.tar.gz (if the file has not been decompressed
already)
# cd openbabel-2.3.2
# mkdir build
# cd build
# cmake ../
# make -j2
# sudo make install

Mac OS X users:

To get OpenBabel to compile, you might find it helpful to get homebrew (http://brew.sh).
Just type:

# ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)"

at a terminal prompt.

Then this program can be used to obtain missing software.
You may need to type the following at terminal prompts to get cmake and pkg_config,
if you do not have them already:

# brew install cmake

# brew install Pkg_Config

For several Linux distributions it is also possible to install
Open Babel using the built-in packaging system, e.g. apt under
Ubuntu Linux. Both the libopenbabel and libopenbabel-dev packages
need to be installed in this case.

The next and final step is the compilation of the molBLOCKS suite.
This is simply accomplished by entering the molblocks directory and
typing:

# make.

In case of errors, it might be necessary to edit the path to the
openbabel library in the Makefile by modifying the following line:

INCLUDES := −Iboost −I/usr/local/include/openbabel−2.0

with the correct location of Open Babel on your system.

-----------------------------------------------------------------------------

INSTALLATION -- VIRTUAL MACHINE

For users who do not wish to or cannot compile molBLOCKS, we prepared an
image of Linux Debian with a pre-installed copy of molBLOCKS 
(http://compbio.princeton.edu/molblocks/download.html).
Right-click (or control-click on Mac OS X) and download the .ova file
containing the virtual machine image.
The image should run out of the box on any virtualization environment,
but we recommend VirtualBox (https://www.virtualbox.org/wiki/Downloads),
which is freely available for Windows, Linux and Mac OS X.

After installing VirtualBox, double-click on the Linux image and 
import the Virtual Machine with standard settings. Alternatively,
choose File->Import Appliance from Virtual Box menu.
More information on importing a Virtual Machine can be found at
https://www.virtualbox.org/manual/ch01.html\#ovf.

After successfully importing the Virtual Machine, start it by pushing the
play button. Once booted, the molBLOCKS program will be in the
molblocks directory, ready for use.

-----------------------------------------------------------------------------

EXAMPLES:

1. The 'exampleCephalosp' folder contains the example described in 
Section 4.2 of the User's Guide.

To run it, please type the following:

# cd molblocks/exampleCephalosp
#
# ../fragment -i cephalosp.smi -r RECAP.txt -n 4 -o cephalosp.frag -e
# ../fragment -i background.smi -r RECAP.txt -n 4 -o background.frag -e
#
# ../analyze -i cephalosp.frag -c 0.7 -o cephalosp.enrichCl07
             -e background.frag

The output directory already contains the output files that should be
obtained after following the steps outlined above.


2. The 'exampleAntiNeopl' folder contains the application outlined in
the Supplementary Material of the accompanying paper.

To run this example, please type the following:

# cd molblocks/exampleAntiNeopl
#
# ../fragment -i antineoplastic.smi -r RECAP.txt -n 4
              -o antineoplastic.frag -e
# ../fragment -i background.smi -r RECAP.txt -n 4 -o background.frag -e
#
# ../analyze -i antineoplastic.frag -c 0.8 -o antineoplastic.enrichCl08
             -c 0.8 -e background.frag

-----------------------------------------------------------------------------

VISUALIZING FRAMENTS:

To visualize the structure of the fragments produced by molBLOCKS,
you can use the freely available Marvin Sketch
(http://www.chemaxon.com/download/marvin-suite).
Fragments can be directly copied and pasted into the main window of
Marvin Sketch.


Please send an email to Dario Ghersi (dghersi [at] princeton.edu) if you
have any questions.
