#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "TH1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"

/// include pdfs
#include "RooCBExGaussShape.h"
#include "RooCMSShape.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"

#include <vector>
#include <string>
#ifdef __CINT__
#pragma link C++ class std::vector<std::string>+;
#endif

using namespace RooFit;
using namespace std;

class FuncFitter {
public:
  FuncFitter( TFile *file, std::string histname  );
  FuncFitter( TH1 *hPass, TH1 *hFail, std::string histname  );
  ~FuncFitter(void) {if( _work != 0 ) delete _work; }
  void setWorkspace(std::vector<std::string>);
  void setOutputFile(TFile *fOut ) {_fOut = fOut;}
  void fits(bool mcTruth,std::string title = "");
  void textParForCanvas(RooFitResult *resP, TPad *p);
  

  void setFitRange(double xMin,double xMax) { _xFitMin = xMin; _xFitMax = xMax; }
private:
  RooWorkspace *_work;
  std::string _histname_base;
  TFile *_fOut =0;
  double _nTotP, _nTotF;
  double _xFitMin,_xFitMax;
};

FuncFitter::FuncFitter(TFile *filein, std::string histname   ) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  _histname_base = histname;  

  TH1 *hPass = (TH1*) filein->Get(TString::Format("%s_Pass",histname.c_str()).Data());
  TH1 *hFail = (TH1*) filein->Get(TString::Format("%s_Fail",histname.c_str()).Data());
  _nTotP = hPass->Integral();
  _nTotF = hFail->Integral();
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
   if(  hPass->GetXaxis()->GetBinCenter(ib) <= 60 || hPass->GetXaxis()->GetBinCenter(ib) >= 120 ) {
     hPass->SetBinContent(ib,0);
     hFail->SetBinContent(ib,0);
   }
  
  _work = new RooWorkspace("w") ;
  _work->factory("x[50,130]");
  _work->factory("index[pass,fail]");

  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  RooDataHist rooFail("hFail","hFail",*_work->var("x"),hFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;
  RooDataHist combData("combData","combined data",*_work->var("x"),RooFit::Index(*_work->cat("index")),
          RooFit::Import("pass",*(RooDataHist*)_work->data("hPass")),RooFit::Import("fail",*(RooDataHist*)_work->data("hFail")));
  _work->import(combData);
  _xFitMin = 60;
  _xFitMax = 120;
}

FuncFitter::FuncFitter(TH1 *hPass, TH1 *hFail, std::string histname  )  {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  _histname_base = histname;
  
  _nTotP = hPass->Integral();
  _nTotF = hFail->Integral();
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
    if(  hPass->GetXaxis()->GetBinCenter(ib) <= 60 || hPass->GetXaxis()->GetBinCenter(ib) >= 120 ) {
      hPass->SetBinContent(ib,0);
      hFail->SetBinContent(ib,0);
    }
  
  _work = new RooWorkspace("w") ;
  _work->factory("x[50,130]");
  _work->factory("index[pass,fail]");

  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  RooDataHist rooFail("hFail","hFail",*_work->var("x"),hFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;
  RooDataHist combData("combData","combined data",*_work->var("x"),RooFit::Index(*_work->cat("index")),
          RooFit::Import("pass",*(RooDataHist*)_work->data("hPass")),RooFit::Import("fail",*(RooDataHist*)_work->data("hFail")));
  _work->import(combData) ;

  _xFitMin = 60;
  _xFitMax = 120;
  
}




void FuncFitter::setWorkspace(std::vector<std::string> workspace) {
  for( unsigned icom = 0 ; icom < workspace.size(); ++icom ) {
    _work->factory(workspace[icom].c_str());
  }
  _work->factory(TString::Format("expr::nSigP('eff*fSig*nTot',eff[0.9,0,1.0],fSig[.9,0,1],nTot[%f,0,%f])",_nTotF+_nTotP,2*(_nTotF+_nTotP)+10));
  _work->factory("expr::nSigF('(1-eff)*fSig*nTot',eff,fSig,nTot)");
  _work->factory("expr::nBkgP('effBkg*(1-fSig)*nTot',effBkg[.5,0,1],fSig,nTot)");
  _work->factory("expr::nBkgF('(1-effBkg)*(1-fSig)*nTot',effBkg,fSig,nTot)");

  _work->factory("SUM::pdfPass(nSigP*sigPass,nBkgP*bkgPass)");
  _work->factory("SUM::pdfFail(nSigF*sigFail,nBkgF*bkgFail)");
  _work->factory("SIMUL:jointModel(index,pass=pdfPass,fail=pdfFail)");

  _work->Print();
}

void FuncFitter::fits(bool mcTruth,string title) {

  cout << " title : " << title << endl;


  RooAbsPdf *pdfPass = _work->pdf("pdfPass");
  RooAbsPdf *pdfFail = _work->pdf("pdfFail");
  _work->var("x")->setRange(_xFitMin,_xFitMax);
  _work->var("x")->setRange("fitMassRange",_xFitMin,_xFitMax);
  RooFitResult* resPass = _work->pdf("jointModel")->fitTo(*_work->data("combData"),Minos(*_work->var("eff")),SumW2Error(mcTruth),Save(true),Extended(true),Range("fitMassRange"));

  RooPlot *pPass = _work->var("x")->frame(60,120);
  RooPlot *pFail = _work->var("x")->frame(60,120);
  pPass->SetTitle("passing probe");
  pFail->SetTitle("failing probe");

  _work->data("hPass") ->plotOn( pPass );
  _work->pdf("pdfPass")->plotOn( pPass, LineColor(kRed) );
  _work->pdf("pdfPass")->plotOn( pPass, Components("bkgPass"),LineColor(kBlue),LineStyle(kDashed));
  _work->data("hPass") ->plotOn( pPass );

  _work->data("hFail") ->plotOn( pFail );
  _work->pdf("pdfFail")->plotOn( pFail, LineColor(kRed) );
  _work->pdf("pdfFail")->plotOn( pFail, Components("bkgFail"),LineColor(kBlue),LineStyle(kDashed));
  _work->data("hFail") ->plotOn( pFail );

  TCanvas c("c","c",1100,450);
  c.Divide(3,1);
  TPad *padText = (TPad*)c.GetPad(1);
  textParForCanvas( resPass, padText );
  c.cd(2); pPass->Draw();
  c.cd(3); pFail->Draw();

  _fOut->cd();
  c.Write(TString::Format("%s_Canv",_histname_base.c_str()),TObject::kOverwrite);
  resPass->Write(TString::Format("%s_resP",_histname_base.c_str()),TObject::kOverwrite);
}





/////// Stupid parameter dumper /////////
void FuncFitter::textParForCanvas(RooFitResult *resP,TPad *p) {

  TPaveText *text1 = new TPaveText(0,0.8,1,1);
  text1->SetFillColor(0);
  text1->SetBorderSize(0);
  text1->SetTextAlign(12);
  text1->AddText(TString::Format("* fit status: %d",resP->status()));
  text1->AddText(TString::Format("* eff = %1.4f #pm %1.4f",_work->var("eff")->getVal(),_work->var("eff")->getError()));

    TPaveText *text = new TPaveText(0,0,1,0.8);
  text->SetFillColor(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->AddText("    --- parmeters " );
  if(resP){
      RooArgList listParFinalP = resP->floatParsFinal();
      for( int ip = 0; ip < listParFinalP.getSize(); ip++ ) {
        TString vName = listParFinalP[ip].GetName();
        if(vName == "eff") continue;
        text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f",
                      vName.Data(),
                      _work->var(vName)->getVal(),
                      _work->var(vName)->getError() ) );
  }

  }
  p->cd();
  text1->Draw();
  text->Draw();
}
