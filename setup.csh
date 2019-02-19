# need to have correct $ROOTSYS, $FASTJET and $PYTHIA8 are needed in Makefile
# In case you have installed all from Aliroot installation including fastjet packages based on https://dberzano.github.io/alice/install-aliroot/
# no need to setup $ROOTSYS, $FASTJET need to have only $PYTHIA8
#source ~/alice/root/alice_v5-34-30/inst/bin/thisroot.sh  from AliROOT
#export FASTJET=$HOME/alice/fastjet/3.1.3/inst  from AliROOT
setenv FASTJET /n/work00/osanmasa/fastjet/fastjet-install
source /n/work00/osanmasa/root/root/bin/thisroot.csh
setenv HEPMC /n/work00/osanmasa/hepmc/hepmc2.06.09-install
setenv HEPPDT /n/work00/osanmasa/HepPDT/HepPDT-3.04.01-install
