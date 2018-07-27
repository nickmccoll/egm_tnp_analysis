#############################################################
########## General settings --> RECO
#############################################################
# flag to be Tested
flags = {
    'passingRECO'       : '(passingRECO   == 1)',
    }

doFunc = True
doRefitFunc=False

baseOutDir = 'results/recoISO/'

#############################################################
########## samples definition  - preparing the samples
#############################################################
### samples are defined in etc/inputs/tnpSampleDef.py
### not: you can setup another sampleDef File in inputs
# import samples as tnpSamples
import etc.inputs.tnpSampleDefComp as tnpSamples
tnpTreeDir = 'tnpEleReco'

samplesDef = {
    'data'   : tnpSamples.Moriond17_80X_reco['data'].clone(),
    'mcNom'  : tnpSamples.Moriond17_80X_reco['DY_madgraph'].clone(),
    'mcAlt'  : tnpSamples.Moriond17_80X_reco['DY_amcatnlo'].clone(),
    'tagSel' : tnpSamples.Moriond17_80X_reco['DY_madgraph'].clone(),
}

samplesDef['data' ].set_tnpTree(tnpTreeDir)
if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_tnpTree(tnpTreeDir)
if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_tnpTree(tnpTreeDir)
if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_tnpTree(tnpTreeDir)

if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_mcTruth()
if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_mcTruth()
if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_mcTruth()
if not samplesDef['tagSel'] is None:
    samplesDef['tagSel'].rename('mcAltSel_DY_madgraph')
    samplesDef['tagSel'].set_cut('tag_Ele_pt > 35')

## set MC weight, simple way (use tree weight) 
weightName = 'totWeight'
if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_weight(weightName)
if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_weight(weightName)
if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_weight(weightName)

#############################################################
########## bining definition  [can be nD bining]
#############################################################
biningDef = [
   { 'var' : 'sc_isoPT/sc_pt-1' , 'type': 'float', 'bins': [-1,0.2,.3,.4,.5,.6,.8,1] },
   { 'var' : 'sc_isoDR'  , 'type': 'float', 'bins': [0,0.01,0.02,0.03,0.05,0.08,.12,.4]},
]

#############################################################
########## Cuts definition for all samples
#############################################################
### cut
# cutBase   = 'sc_pt >30 && abs(sc_eta) < 2.5 && tag_Ele_pt > 30 && abs(tag_sc_eta) < 2.17 && tag_Ele_nonTrigMVA >0.96 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 60'
#running on skim
cutBase   = 'sc_pt >30 && abs(sc_eta) < 2.5'

#### or remove any additional cut (default)
additionalCuts = None


#############################################################
########## fitting params to tune fit by hand if necessary
#############################################################
if doFunc:
    cmsShapePars = [
        "acmsP[70.0]","betaP[0.05,0.01,0.08]","gammaP[0.1, 0, 1]","peakP[90.0]",
        "acmsF[70.0]","betaF[0.05,0.01,0.08]","gammaF[0.1, 0, 1]","peakF[90.0]"
    ]
    expoPars = ["alphaP[0.,-5.,5.]","alphaF[0.,-5.,5.]"]
    zGausPars = ["peakPassMean[90,88,92]","peakPassSigma[3,2,4]","peakFailMean[90,88,92]","peakFailSigma[3,2,4]"]
    radGausPars = ["radMean[70,60,82]","radSigma[4,3,8]","radFracPass[0.5,0.0,1.0]","radFracFail[0.5,0.0,1.0]"]            
    cmsShapePDF = ["RooCMSShape::expoPass(x, acmsP, betaP, gammaP, peakP)","RooCMSShape::expoFail(x, acmsF, betaF, gammaF, peakF)"]
    expoPDF   = ["Exponential::expoPass(x, alphaP)","Exponential::expoFail(x, alphaF)"]
    zGausPDF = ["Gaussian::peakPass(x,peakPassMean,peakPassSigma)","Gaussian::peakFail(x,peakFailMean,peakFailSigma)"]
    radGausPDF = ["Gaussian::radPass(x,radMean,radSigma)","Gaussian::radFail(x,radMean,radSigma)"]
    stdSigPDF = ["Gaussian::sigPass(x,peakPassMean,peakPassSigma)","Gaussian::sigFail(x,peakFailMean,peakFailSigma)"]
    altSigPDF = ["SUM::sigPass(radFracPass*radPass,peakPass)","SUM::sigFail(radFracFail*radFail,peakFail)"]
    stdBkgPDF = ["SUM::bkgPass(radFracPass*radPass,expoPass)","SUM::bkgFail(radFracFail*radFail,expoFail)"]
    altSigBkgPDF = ["RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)","RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)"]
    altBkgPDF = ["SUM::bkgPass(radFracPass*radPass,expoPass)","SUM::bkgFail(radFracFail*radFail,expoFail)"]
    
    if doRefitFunc :
        zGausPars = ["peakPassMean[90,80,92]","peakPassSigma[3,2,5.5]","peakFailMean[90,80,92]","peakFailSigma[3,2,5.5]"]
        radGausPars = ["radMean[70,70,80]","radSigma[4,3,8]","radFracPass[0.5,0.0,1.0]","radFracFail[0.5,0.0,1.0]"]            
                
    stdPars = []
    stdPars.extend(cmsShapePars)
    stdPars.extend(expoPars)
    stdPars.extend(zGausPars)
    stdPars.extend(radGausPars)
    
    tnpParNomFit = []
    tnpParNomFit.extend(stdPars)
    tnpParNomFit.extend(cmsShapePDF)
    tnpParNomFit.extend(radGausPDF)
    tnpParNomFit.extend(stdBkgPDF)
    tnpParNomFit.extend(stdSigPDF)
    
    tnpParAltSigFit = []
    tnpParAltSigFit.extend(stdPars)
    tnpParAltSigFit.extend(altSigBkgPDF)
    tnpParAltSigFit.extend(radGausPDF)
    tnpParAltSigFit.extend(zGausPDF)
    tnpParAltSigFit.extend(altSigPDF)
    
    tnpParAltBkgFit = []
    tnpParAltBkgFit.extend(stdPars)
    tnpParAltBkgFit.extend(expoPDF)
    tnpParAltBkgFit.extend(radGausPDF)
    tnpParAltBkgFit.extend(altBkgPDF)    
    tnpParAltBkgFit.extend(stdSigPDF)
        
else :             
    tnpParNomFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,5.0]",
        "meanF[-0.0,-5.0,5.0]","sigmaF[0.5,0.1,5.0]",
        "acmsP[60.,50.,80.]","betaP[0.05,0.01,0.08]","gammaP[0.1, 0, 1]","peakP[90.0]",
        "acmsF[60.,50.,80.]","betaF[0.05,0.01,0.08]","gammaF[0.1, 0, 1]","peakF[90.0]",
        ]
     
    tnpParAltSigFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[1,0.7,6.0]","alphaP[2.0,1.2,3.5]" ,'nP[3,-5,5]',"sigmaP_2[1.5,0.5,6.0]","sosP[1,0.5,5.0]",
        "meanF[-0.0,-5.0,5.0]","sigmaF[2,0.7,15.0]","alphaF[2.0,1.2,3.5]",'nF[3,-5,5]',"sigmaF_2[2.0,0.5,6.0]","sosF[1,0.5,5.0]",
        "acmsP[60.,50.,75.]","betaP[0.04,0.01,0.06]","gammaP[0.1, 0.005, 1]","peakP[90.0]",
        "acmsF[60.,50.,75.]","betaF[0.04,0.01,0.06]","gammaF[0.1, 0.005, 1]","peakF[90.0]",
        ]
          
    tnpParAltBkgFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,5.0]",
        "meanF[-0.0,-5.0,5.0]","sigmaF[0.5,0.1,5.0]",
        "alphaP[0.,-5.,5.]",
        "alphaF[0.,-5.,5.]",
        ]