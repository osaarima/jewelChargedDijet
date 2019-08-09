# need to have correct $ROOTSYS, $FASTJET and $PYTHIA8 are needed in Makefile
# In case you have installed all from Aliroot installation including fastjet packages based on https://dberzano.github.io/alice/install-aliroot/
# no need to setup $ROOTSYS, $FASTJET need to have only $PYTHIA8
#source ~/alice/root/alice_v5-34-30/inst/bin/thisroot.sh  from AliROOT
#export FASTJET=$HOME/alice/fastjet/3.1.3/inst  from AliROOT
#export PYTHIA8=$HOME/alice/PythiaMy/pythia8226

#export HEPMC=$HOME/alice/HepMCMy/hepmc2.06.09-install/
#export HEPPDT=$HOME/alice/HepPDTMy/HepPDT-3.04.01-install

# In puck:
export FASTJET=/n/work01/osanmasa/fastjet/fastjet-install
source /n/work01/osanmasa/root/root-6.10.08-build/bin/thisroot.sh
export HEPMC=/n/work01/osanmasa/hepmc/hepmc2.06.09-install
export HEPPDT=/n/work01/osanmasa/HepPDT/HepPDT-3.04.01-install
