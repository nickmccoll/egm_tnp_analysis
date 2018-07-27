from libPython.tnpClassUtils import tnpSample

### qll stat
recoDir = 'reco_skim/'

Moriond17_80X_reco = {
    'DY_madgraph' : tnpSample('DY_madgraph', recoDir + 'TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_DYToLL_madgraph.root', 
                              isMC = True),
    'DY_amcatnlo' : tnpSample('DY_amcatnlo', recoDir + 'TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root', 
                              isMC = True),
    'data'  : tnpSample('data'  , recoDir + 'TnPTree_data.root' , lumi = -1 )
    }


idDir = 'id_skim/'

Moriond17_80X_id = {
    'DY_madgraph' : tnpSample('DY_madgraph', idDir + 'TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_DYToLL_madgraph.root', 
                              isMC = True),
    'DY_amcatnlo' : tnpSample('DY_amcatnlo', idDir + 'TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root', 
                              isMC = True),
    'data'  : tnpSample('data'  , idDir + 'TnPTree_data.root' , lumi = -1 )
    }
