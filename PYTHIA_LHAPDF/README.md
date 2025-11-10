### How to install PYTHIA and LHAPDF together

### (1) Installing LHAPDF
    wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.5.tar.gz
    mv index.html?f=LHAPDF-6.5.5.tar.gz LHAPDF-6.5.5.tar.gz
    tar -xvzf LHAPDF-6.5.5.tar.gz
    mkdir LHAPDF-6.5.5Install
    cd LHAPDF-6.5.5
    ./configure --prefix=/nfs/zfs2/Software/LHAPDF-6.5.5Install --enable-shared PYTHON=/usr/bin/python3.10
    make
    make install
    
### (2) Now download NuclearPDF sets for proton and Lead Ion
    cd LHAPDF-6.5.5Install/share/LHAPDF/

    wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/nCTEQ15_1_1.tar.gz
    wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/nCTEQ15_208_82.tar.gz

    tar -xvzf nCTEQ15_1_1.tar.gz
    tar -xvzf nCTEQ15_208_82.tar.gz
    
### (3) Download and Install PYTHIA with LHAPDF8.309 for x-scape
        #4th argument of the function LHAGrid1 is different between pythia 8.309 and 8.316
        #https://pythia.org/doxygen/pythia8309/classPythia8_1_1LHAGrid1.html
	#https://pythia.org/latest-doxygen/classPythia8_1_1LHAGrid1.html
        #The X-SCAPE-COMP compiles and runs using pythia 8.309
    wget https://pythia.org/download/pythia83/pythia8309.tar
    tar -xvf pythia8309.tar
    cd pythia8309
    ./configure --prefix=/nfs/zfs2/Software/pythia8309  --with-lhapdf6=/nfs/zfs2/Software/LHAPDF-6.5.5Install
    make
    
### (4) Goto examples directory and check sample codes
    cd examples/
    emacs -nw main01.cc
    add pythia.readString("PDF:pSetB = LHAPDF6:nCTEQ15_208_82"); inside main01.cc file
    make main01
    ./main01


    
