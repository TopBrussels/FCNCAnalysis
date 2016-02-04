from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TopTree_FCNC_tHWW_1' # this will be the task name so it is better to be adapted for every new submission
config.General.workArea = 'crab_projects_TopTree_tHWW_1'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'TOPTREE_fromminiAOD.py' #modify if needed according to which configartion file you need to be used 
#config.JobType.inputFiles = ['MCRUN2_74_V7::All']
config.JobType.outputFiles = ['TOPTREE.root']

config.section_("Data")
config.Data.inputDataset = '/tHFCNC13TeV/kskovpen-kskovpen_tHToWW_1L_Kappa_hut_AODFASTSIM_4d8f2bc2e6b46f86d9c660f8ab1d998e_USER-881b16c8bee650855968d7f6ca374946/USER' #modify if needed according to which dataset you want to run over
#config.Data.inputDBS = 'global'
config.Data.outLFNDirBase = '/store/user/sabuzeid/FCNC_MCSamples_Kirill/tHFCNC13TeV/kskovpen-kskovpen_tHToWW_1L_Kappa_hut_AODFASTSIM_4d8f2bc2e6b46f86d9c660f8ab1d998e_USER-881b16c8bee650855968d7f6ca374946/USER/' #modify according to your stoarge element 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5 #this is number of events per job 
#config.Data.totalUnits = 100
config.Data.publication = False
#config.Data.Data.publishDBS = 'phys03'
#config.Data.publishDataName = 'TEST'
#config.Data.ignoreLocality = True #activate this line if you want to run remotely (i.e xroot)

config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE' # to be adapted according to your network  

config.section_("User") # to be adapted according to your network
config.User.voGroup = 'becms' # to be adapted according to your network
