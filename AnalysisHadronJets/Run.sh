
PYTHIA8=/nfs/zfs2/Software/pythia8309

echo g++ -std=c++17 Test.C -o Test \
$(root-config --cflags) \
-I${PYTHIA8}/include \
-L${PYTHIA8}/lib -Wl,-rpath,${PYTHIA8}/lib \
-lpythia8 -ldl \
$(root-config --libs)


#g++ -std=c++17 AnalysisSpectraSingleHadron.C -o Test $(root-config --cflags) -I${PYTHIA8}/include -L${PYTHIA8}/lib -Wl,-rpath,${PYTHIA8}/lib -lpythia8 -ldl $(root-config --libs) 


g++ -std=c++17 AnalysisSpectraSingleHadron.C -o Test -pthread -std=c++17 -m64 -I/snap/root-framework/950/usr/local/include -I/nfs/zfs2/Software/pythia8309/include -L/nfs/zfs2/Software/pythia8309/lib -Wl,-rpath,/nfs/zfs2/Software/pythia8309/lib -lpythia8 -ldl -L/snap/root-framework/950/usr/local/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/snap/root-framework/950/usr/local/lib -pthread -lm -ldl -rdynamic


# g++  -std=c++17 AnalysisSpectraSingleHadron.C /nfs/zfs2/Software/pythia8309/lib/libpythia8.a -o AnalysisSpectraSingleHadron  -I/nfs/zfs2/Software/pythia8309/include -O2  -pedantic -W -Wall -Wshadow -fPIC -L/nfs/zfs2/Software/pythia8309/lib -Wl,-rpath,/nfs/zfs2/Software/pythia8309/lib -lpythia8 -ldl \
#-Wl,-rpath, -w  -pthread -std=c++17 -m64 -I/snap/root-framework/950/usr/local/include -I/nfs/zfs2/Software/pythia8309/include \
#-L/nfs/zfs2/Software/pythia8309/lib -Wl,-rpath,/nfs/zfs2/Software/pythia8309/lib -lpythia8 -ldl -L/snap/root-framework/950/usr/local/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/snap/root-framework/950/usr/local/lib -pthread -lm -ldl -rdynamic
