
#include <vector>
#include <string>
#include <utility>
#include <fstream>

#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#ifdef __CINT__
#pragma link C++ class std::vector<std::string>+;
#endif

using namespace std;

class Uncertainty {
public:
    enum UType {MC, DATA, CORR, STAT};
    std::string name="";
    UType type=STAT;
    TH2* mcEff = 0;
    TH2* dataEff = 0;
    TH2* unc  = 0;
    Uncertainty(){}
    Uncertainty(const std::string& name, const UType type) :
        name(name),type(type){
    }
    ~Uncertainty(){
//        if(mcEff) delete mcEff;
//        if(dataEff) delete dataEff;
//        if(unc) delete unc;
    }

};

class SFPlotter {
public:

    SFPlotter(const std::string& xTitle,const std::vector<double>& xBins,
            const std::string& yTitle, const std::vector<double>& yBins ) :
                xTitle(xTitle),yTitle(yTitle),xBins(xBins),yBins(yBins),nBinsX(xBins.size()-1)
,nBinsY(yBins.size()-1)
{
        bins.reserve(xBins.size()*yBins.size());
        dataEff = makeHistogram("EGamma_EffData2D");
        mcEff   = makeHistogram("EGamma_EffMC2D");
        sf = makeHistogram("EGamma_SF2D");
        totUnc = makeHistogram("EGamma_UncTot2D");
}
    ~SFPlotter() {
//        if( dataEff) delete dataEff;
//        if( mcEff  ) delete mcEff  ;
//        if( sf     ) delete sf     ;
//        if( totUnc ) delete totUnc ;
    }

    void addBin(const std::string& binName, const float xV, const float yV){
        bins.emplace_back(binName,dataEff->FindFixBin(xV,yV));

    }

    /// Must add all bins before adding nominal or systematics
    void addNominal(const std::string& dataFileName,const std::string& mcFileName){
        auto dataFile = new TFile( dataFileName.c_str(),"READ" );
        fillResults(dataFile,dataEff);
        delete dataFile;
        auto mcFile = new TFile( mcFileName.c_str(),"READ" );
        fillResults(mcFile,mcEff);
        delete mcFile;
    }

    void addDataOnlySyst(const std::string& name, const std::string& dataFileName){
        addSyst(name,Uncertainty::DATA,dataFileName,"");
    }
    void addMCOnlySyst(const std::string& name, const std::string& mcFileName){
        addSyst(name,Uncertainty::MC,"",mcFileName);
    }
    void addCorrSyst(const std::string& name, const std::string& dataFileName,const std::string& mcFileName){
        addSyst(name,Uncertainty::CORR,dataFileName,mcFileName);
    }
    void process(){
        //first get the sf values and stat uncertainties
        auto dataStat = makeHistogram("EGamma_UncStatData2D");
        auto mcStat   = makeHistogram("EGamma_UncStatMC2D");
        for(unsigned int iX = 1; iX <=nBinsX;++iX)
            for(unsigned int iY = 1; iY <= nBinsY;++iY){
                const auto eM    = mcEff->GetBinContent(iX,iY);
                const auto eMErr = mcEff->GetBinError(iX,iY);
                const auto eD    = dataEff->GetBinContent(iX,iY);
                const auto eDErr = dataEff->GetBinError(iX,iY);
                if(eM == 0) continue;
                sf->SetBinContent(iX,iY, eD/eM  );
                dataStat->SetBinContent(iX,iY,eDErr/eM);
                mcStat->SetBinContent(iX,iY,eMErr*eD/(eM*eM));
            }
        //process systematics
        for(auto& unc : uncs){ procSyst(unc);}
        //add in statistical
        uncs.emplace_back("dataStat",Uncertainty::STAT);
        uncs.back().unc = dataStat;
        uncs.emplace_back("mcStat",Uncertainty::STAT);
        uncs.back().unc = mcStat;
        //add them up

        for(unsigned int iX = 1; iX <=nBinsX;++iX)
            for(unsigned int iY = 1; iY <= nBinsY;++iY){
                double total = 0;
                for(auto& unc : uncs) total += unc.unc->GetBinContent(iX,iY)*unc.unc->GetBinContent(iX,iY);
                total = std::sqrt(total);
                totUnc->SetBinContent(iX,iY,total);
                sf->SetBinError(iX,iY,total);
            }
    }
    void makeOutput(const std::string& outFileName){
        makeTxtFile(outFileName);
        makePDF(outFileName+"_egammaPlots.pdf");
        makeMainRootFile(outFileName +"_EGM2D.root");
    }

private:
    void makePDF(const TString& outFileName){
        auto cDummy = new TCanvas();
        cDummy->Print(outFileName+"[");
        gStyle->SetPalette(1);
        gStyle->SetPaintTextFormat("1.3f");
        gStyle->SetOptTitle(1);
        gStyle->SetOptStat(0);
        auto mkCan = [&](const std::string& name,const std::string& title) ->TCanvas* {
            auto c2D = new TCanvas(name.c_str(),title.c_str(),900,600);
            c2D->Divide(2,1)                         ;
            c2D->GetPad(1)->SetRightMargin(0.15)      ;
            c2D->GetPad(1)->SetLeftMargin( 0.15)      ;
            c2D->GetPad(1)->SetTopMargin(  0.10)      ;
            c2D->GetPad(2)->SetRightMargin(0.15)      ;
            c2D->GetPad(2)->SetLeftMargin( 0.15)      ;
            c2D->GetPad(2)->SetTopMargin(  0.10)      ;
            return c2D;
        };

        sf->SetTitle("e/#gamma scale factors")        ;
        totUnc->SetTitle ("e/#gamma uncertainties")        ;
        auto c2D  = mkCan("canScaleFactor","canScaleFactor");
        c2D->cd(1)                                 ;
        double dmin = 1.0 - sf->GetMinimum()      ;
        double dmax = sf->GetMaximum() - 1.0             ;
        double dall = std::max(dmin,dmax)                      ;
        sf->SetMinimum(1-dall)                    ;
        sf->SetMaximum(1+dall)                    ;
        sf->DrawCopy("colz TEXT45")               ;

        c2D->cd(2);
        totUnc->SetMinimum(0)                                 ;
        totUnc->SetMaximum(std::min(totUnc->GetMaximum(),0.2))     ;

        totUnc->DrawCopy("colz TEXT45")                       ;
        c2D->Print( outFileName );

        for(auto& unc:uncs){
            std::string  nm = "canScaleFactor";
            auto c2D_Err  = mkCan(nm+"_"+unc.name,nm+": "+unc.name);
            auto uncRel = (TH2*)unc.unc->Clone();
            uncRel->SetDirectory(0);
            uncRel->Divide(totUnc);
            unc.unc->SetMinimum(0)                                                           ;
            unc.unc->SetMaximum(min(unc.unc->GetMaximum(),0.2))                         ;
            uncRel->SetMinimum(0)                                                           ;
            uncRel->SetMaximum(1)                                                           ;
            unc.unc->SetTitle((std::string("e/#gamma absolute SF syst: ")+unc.name).c_str())        ;
            uncRel->SetTitle ((std::string("e/#gamma relative SF syst: ")+unc.name).c_str())        ;
            c2D_Err->cd(1)                                                                         ;
            unc.unc->DrawCopy("colz TEXT45")                                                 ;
            c2D_Err->cd(2)                                                                         ;
            uncRel->DrawCopy("colz TEXT45")                                                 ;
            c2D_Err->Print(outFileName)                                                                ;
            delete  uncRel;
        }
        cDummy->Print( outFileName + "]" );
    }
    void makeMainRootFile(const std::string& outFileName){
        auto oF = new TFile(outFileName.c_str(),"recreate");
        oF->cd();
        dataEff->Write();
        mcEff->Write();
        sf->Write();
        oF->Close();

    }
    void makeTxtFile(const std::string& outFileName){
        std::ofstream outF(outFileName.c_str(),std::ios::out|std::ios::trunc);
        outF <<"#"<< xTitle <<" min\t"<< xTitle <<" max"
                << yTitle <<" min\t"<< yTitle <<" max";
        outF<<"\tDataEff\tDataStat\tMCEff\tMCStat";
        for(const auto& unc :uncs){
            if(unc.type == Uncertainty::DATA) outF <<"\t"<<unc.name;
            if(unc.type == Uncertainty::MC) outF <<"\t"<<unc.name;
            if(unc.type == Uncertainty::CORR)
                outF <<"\t"<<unc.name+"Data\t"<<unc.name+"MC";
        }
        outF <<"\n";

        for(unsigned int iX = 1; iX <=nBinsX;++iX)
            for(unsigned int iY = 1; iY <= nBinsY;++iY){
            outF << TString::Format("%+8.3f\t%+8.3f\t%+8.3f\t%+8.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f",
                    sf->GetXaxis()->GetBinLowEdge(iX),sf->GetXaxis()->GetBinLowEdge(iX)+sf->GetXaxis()->GetBinWidth(iX),
                    sf->GetYaxis()->GetBinLowEdge(iY),sf->GetYaxis()->GetBinLowEdge(iY)+sf->GetYaxis()->GetBinWidth(iY),
                    dataEff->GetBinContent(iX,iY),dataEff->GetBinError(iX,iY),mcEff->GetBinContent(iX,iY),mcEff->GetBinError(iX,iY));
            for(const auto& unc :uncs){
                if(unc.type == Uncertainty::DATA) outF <<TString::Format("\t%5.3f",unc.dataEff->GetBinContent(iX,iY));
                if(unc.type == Uncertainty::MC) outF <<TString::Format("\t%5.3f",unc.mcEff->GetBinContent(iX,iY));
                if(unc.type == Uncertainty::CORR)
                    outF<<TString::Format("\t%5.3f\t%5.3f",unc.dataEff->GetBinContent(iX,iY),unc.mcEff->GetBinContent(iX,iY));
            }
                outF <<"\n";
            }
        outF.close();
    }
    void procSyst(Uncertainty& unc){
        for(unsigned int iX = 1; iX <=nBinsX;++iX)
            for(unsigned int iY = 1; iY <= nBinsY;++iY){
                const double sfv = sf->GetBinContent(iX,iY);
                const double eMC = mcEff->GetBinContent(iX,iY);
                const double eData = dataEff->GetBinContent(iX,iY);
                double uncVal = 0;
                if(unc.type == Uncertainty::MC){
                    if(eMC==0)continue;
                    uncVal =  sfv*std::fabs(unc.mcEff->GetBinContent(iX,iY) -eMC)/eMC;
                } else if(unc.type == Uncertainty::DATA){
                    if(eMC==0)continue;
                    uncVal =  std::fabs(unc.dataEff->GetBinContent(iX,iY) -eData)/eMC;
                } else{
                    if(eMC == 0 || unc.mcEff->GetBinContent(iX,iY) ==0 ) continue;
                    uncVal =  std::fabs(unc.dataEff->GetBinContent(iX,iY)/unc.mcEff->GetBinContent(iX,iY) -sfv);
                }
                unc.unc->SetBinContent(iX,iY,uncVal );
            }
    }
    void addSyst(const std::string& name, const Uncertainty::UType type, const std::string& dataFileName,const std::string& mcFileName){
        std::string hSystName = "EGamma_Unc" + name + "2D";
        uncs.emplace_back(name,type);
        uncs.back().unc = makeHistogram(hSystName);

        if(type == Uncertainty::DATA ||type == Uncertainty::CORR){
            std::string hDataName = "EGamma_Eff" + name + "Data2D";
            uncs.back().dataEff = makeHistogram(hDataName);
            auto dataFile = new TFile( dataFileName.c_str(),"READ" );
            fillResults(dataFile,uncs.back().dataEff);
            delete dataFile;
        }
        if(type == Uncertainty::MC ||type == Uncertainty::CORR){
            std::string hMCName = "EGamma_Eff" + name + "MC2D";
            uncs.back().mcEff = makeHistogram(hMCName);
            auto mcFile = new TFile( mcFileName.c_str(),"READ" );
            fillResults(mcFile,uncs.back().mcEff);
            delete mcFile;
        }



    }
    void fillResults(TFile* inFile, TH2* h){
        for(const auto& b : bins){
            RooFitResult * fitRes = 0;
            inFile->GetObject((b.first +"_resP").c_str(),fitRes);
            if(fitRes==0) continue;
            auto fitP = (RooRealVar*)fitRes->floatParsFinal().find("eff");
            h->SetBinContent(b.second,fitP->getVal());
            h->SetBinError(b.second,fitP->getError());
            delete fitRes;
        }
    }

    TH2 * makeHistogram(const std::string& name, const std::string& title=""){
        auto h2 = new TH2F(name.c_str(),(title+";"+xTitle+";"+yTitle).c_str(),
                nBinsX, &xBins[0],nBinsY,&yBins[0]);
        h2->Sumw2(true);
        h2->SetDirectory(0);
        return h2;
    }

    const std::string xTitle;
    const std::string yTitle;
    const std::vector<double> xBins;
    const std::vector<double> yBins;
    const unsigned int nBinsX;
    const unsigned int nBinsY;
    std::vector<std::pair<std::string,int>> bins;
    std::vector<Uncertainty> uncs;
    TH2* dataEff;
    TH2* mcEff;
    TH2* sf;
    TH2* totUnc;

};
