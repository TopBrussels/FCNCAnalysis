from ROOT import TChain
from glob import glob

path = '/pnfs/iihe/cms/store/user/kderoove/FCNC_Run2_Samples/tHToBB_1L_Kappa_hct/tHFCNC13TeV/crab_tHToBB_1L_Kappa_hct/151014_120023/0000/TOPTREE_*.root'
files = glob(path)
root_files = []
for f in files:
	root_files.append('dcap://maite.iihe.ac.be' + f)
print root_files
chain = TChain('eventTree')
for rf in root_files:
	chain.Add(rf)
print 'added files'
nEntries = chain.GetEntries();
print nEntries

