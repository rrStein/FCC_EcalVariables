import ROOT, sys
from bin.tools import train

#### input variables to be used for training ###
variables = [
    'e2max',
    'eocore',
    'edmax',
    'w3st',
    'w21st',
    #'Jet_trk02_SD_Corr_m',
    #'Jet_trk04_SD_Corr_m',
    #'Jet_trk08_SD_Corr_m',
    # test
    #'Jet_Flow15',
    #'Jet_Flow25',
    #'Jet_Flow35',
    #'Jet_Flow45',
    #'Jet_Flow55',
]

#### define signal and background trees + selection to be applied ###

eospath = '/afs/cern.ch/user/r/rastein/TMVAFacility/testfiles/'
#eospath = '/afs/cern.ch/user/d/djamin/fcc_work/heppy/FCChhAnalyses/output/tagger/W_top_vs_QCD_tagger/'
# to make these special trees ,execute scripts/do_jet_list.py
bkgTree = eospath + 'pi0_20k.root'
sigTree = eospath + 'photon_20k.root'

# take only hadronic decays
SIGcuts = ''
#SIGcuts = 'fullhad_fullhadsemilep_lep_decays<2 && Jet_trk02_tau21>0 && Jet_trk02_tau31>0 && Jet_trk02_tau32>0'
BKGcuts = SIGcuts
cuts    = [SIGcuts,BKGcuts]

#### define MVA method ###

method = 'BDT'
#nTrain = 50000
nTrain = 0
label = 'thad_vs_QCD'


#### train the MVA ####
train(bkgTree, sigTree, variables, method, nTrain, label, cuts)


#### check the output with the TMVAGUI ####

#outFileName = 'TMVA_{}_{}.root'.format(method, label)
#ROOT.TMVA.TMVAGui(outFileName, label)
#raw_input('Press Enter to continue...')
