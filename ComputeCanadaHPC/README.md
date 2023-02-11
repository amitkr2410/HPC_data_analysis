

### To run JETSCAPE software package in Compute Canada cluster(Beluga,Narval) high-performance computing facilities, follow the command outlined below:

(1) Edit the directory structure and provide proper path of data folders, input XML file and output DIR: \
    Open and Edit "python run_jetscape.py" 
    
(2) Set appropriate flags in the reference input XML file: \
    Open and Edit "jetscape_user_pbpb_grid.xml"

(3) Download JETSCAPE code package from git-hub and compile the software: \
    git clone -b main https://github.com/JETSCAPE/JETSCAPE.git \
    cd JETSCAPE \
    mkdir build \
    cd external_packages \
    ./get_music.sh	 \
    ./get_iSS.sh \
    ./get_lbtTab.sh \
    cd ../build \
    cmake ..   \
    make -j4 \
    
(4) To Launch the job, edit file "run_jetscape.py" and set directory and run the command below: \
    python run_jetscape.py \
    
