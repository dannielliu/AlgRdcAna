#include "DmpAlgRdcAna.h"
#include "DmpDataBuffer.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "DmpBgoBase.h"
//#include "mysql.h"
using std::string;

//-------------------------------------------------------------------
DmpAlgRdcAna::DmpAlgRdcAna()
 :DmpVAlg("RdcAna"),
 eventID(0),
 fHeader(0),
 fBgo(0),
 h1(0)
{
  time_low=-1;
	time_up=-1;
	adc_low=-300;
	adc_up=300;
	OptMap.insert(std::make_pair("RunMode",0));
	OptMap.insert(std::make_pair("SetRange/Low",1));
	OptMap.insert(std::make_pair("SetRange/Up",2));
}

//-------------------------------------------------------------------
DmpAlgRdcAna::~DmpAlgRdcAna(){
}

//-------------------------------------------------------------------
bool DmpAlgRdcAna::Initialize(){
  fHeader = new DmpEvtHeader();
  gDataBuffer->ReadObject("Event/Rdc/EventHeader", fHeader);
	fBgo = new DmpEvtBgoRaw();
  gDataBuffer->ReadObject("Event/Rdc/Bgo",fBgo);

  h1 = new TH1D("h1","time",100,time_low,time_up);
  // for pedestal mode
	if (runMode=="Pedestal"||runMode=="Ped"||runMode=="pedestal"){
	  InitializePedMode();
	}
	if (runMode == "Dynode"||runMode=="DynodeRatio"){
	  InitializeDynodeRatioMode();
	}

	std::cout<<"run mode is "<<runMode<<std::endl;
	return true;
}

//-------------------------------------------------------------------
bool DmpAlgRdcAna::InitializePedMode()
{
  for(int layer=0;layer<14;layer++){
    for(int barid=0;barid<24;barid++){
	    for(int sidedy=0;sidedy<6;sidedy++){
		    char tmp[100];
			  sprintf(tmp,"layer_%d_Bar_%d_side_%d_dynode_%d",layer,barid,sidedy/3,(sidedy%3)*3+2);
		    histped[layer][barid][sidedy]=new TH1D(tmp,tmp,200,-300,300);
		  }
	  }
  }
	return true;
}

//-------------------------------------------------------------------
bool DmpAlgRdcAna::InitializeDynodeRatioMode()
{
  for(int layer=0;layer<14;layer++){
    for(int barid=0;barid<22;barid++){
	    for(int side=0;side<2;side++){
		    char tmp[100];
			  sprintf(tmp,"layer_%d_Bar_%d_side_%d",layer,barid,side);
		    dynodeRatio[layer][barid][side]=new TH2D(tmp,tmp,200,adc_low,adc_up/10,200,adc_low,adc_up);
				sprintf(tmp,"dynode %d",dy1);
				dynodeRatio[layer][barid][side]->GetXaxis()->SetTitle(tmp);
				sprintf(tmp,"dynode %d",dy2);
				dynodeRatio[layer][barid][side]->GetYaxis()->SetTitle(tmp);
		  }
	  }
  }
	return true;
}

//-------------------------------------------------------------------
bool DmpAlgRdcAna::ProcessThisEvent(){
  if(time_low==-1)
	  time_low=fHeader->GetSecond();
	time_up=fHeader->GetSecond();
	time_cur=time_up;

  h1->Fill(time_cur);
  //std::cout<<"debuging... current event id"<<eventID++<<std::endl;
	if (runMode=="Pedestal"||runMode=="Ped"||runMode=="pedestal"){
    ProcessPedMode();
  }
	if (runMode == "Dynode"||runMode=="DynodeRatio"){
	  ProcessDynodeRatioMode(dy1,dy2);
	}

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgRdcAna::ProcessPedMode()
{
  int layer,barid,sidedy;
  int evtsize;
  short gid,signal;
  evtsize = fBgo->GetSignalSize();
  for(short i=0;i<evtsize;i++){
    fBgo->GetSignal(i,gid,signal);
	  layer = DmpBgoBase::GetLayerID(gid);
		barid = DmpBgoBase::GetBarID(gid);
		// care for dynode ID
		sidedy = DmpBgoBase::GetSideID(gid)*3+(DmpBgoBase::GetDynodeID(gid)+1)/3-1;
    histped[layer][barid][sidedy]->Fill(signal);
	}
}

//-------------------------------------------------------------------
bool DmpAlgRdcAna::ProcessDynodeRatioMode(short dy1,short dy2)
{
  //std::cout<<"dy1 is "<<dy1<<" dy2 is "<<dy2<<std::endl;
  int layer,barid,side,dynode;
  int bar[14][22][2][2];
	memset(bar,-300,sizeof(bar));
  int evtsize;
  short gid,signal;
  evtsize = fBgo->GetSignalSize();
	// reading an event to a buffer
  for(short i=0;i<evtsize;i++){
    fBgo->GetSignal(i,gid,signal);
	  layer = DmpBgoBase::GetLayerID(gid);
		barid = DmpBgoBase::GetBarID(gid);
		// care for dynode ID
		side = DmpBgoBase::GetSideID(gid);
		dynode = DmpBgoBase::GetDynodeID(gid);
  //std::cout<<layer<<'\t'<<barid<<'\t'<<side<<'\t'<<dynode<<'\t'<<signal<<std::endl;
		if(barid<22){
		  if(dynode==dy1){
		    bar[layer][barid][side][0]=signal;
		  }
		  if(dynode==dy2){
		    bar[layer][barid][side][1]=signal;
		  }
		}
	}
	for(int i=0;i<14;i++){
	for(int j=0;j<22;j++)
	for(int k=0;k<2;k++)
	  if(bar[i][j][k][0]!=-300 && bar[i][j][k][1]!=-300){
      dynodeRatio[i][j][k]->Fill(bar[i][j][k][0],bar[i][j][k][1]);
    }
	}
	
}

//-------------------------------------------------------------------
bool DmpAlgRdcAna::Finalize(){
  //draw time distribution
  TCanvas *c1=new TCanvas("c1","time",800,600);
  h1->Draw();
	c1->Print("TimeDisplay.pdf");
	delete c1;
	delete h1;
  
	// draw pedestals
	if (runMode=="Pedestal"||runMode=="Ped"||runMode=="pedestal"){
    FinalizePedMode();
	}
	if (runMode == "Dynode"||runMode=="DynodeRatio"){
	  FinalizeDynodeRatioMode();
	}

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgRdcAna::FinalizePedMode()
{	
	TSpectrum *spectrum=new TSpectrum();
	int npeak;
	float *peakValue;
	TF1 *fit;
	double pedestalV,pedestalE;
  string peds="pedestal.pdf";
  string peds_start=peds+"[";
  string peds_stop=peds+"]";
  TCanvas *c2 = new TCanvas("c2","peds",0,0,800,600);
	c2->Divide(3,2);
  c2->Print(peds_start.c_str());//,"Portrait");
	for(int layer=0;layer<14;layer++){
	  for(int barid=0;barid<24;barid++){
		  for(int sidedy=0;sidedy<6;sidedy++){
				c2->cd(sidedy+1);
				histped[layer][barid][sidedy]->Draw();
				npeak = spectrum->Search(histped[layer][barid][sidedy]);
				if(npeak>0){
				  peakValue = spectrum->GetPositionX();
					histped[layer][barid][sidedy]->Fit("gaus","","",peakValue[0]-50,peakValue[0]+50);
					fit = histped[layer][barid][sidedy]->GetFunction("gaus");
          pedestalV = fit->GetParameter(1);
          //pedestalE = fit->GetParError(1);
          pedestalE = fit->GetParameter(2);
				}
				else{
				  std::cout<<"peak number is "<<npeak<<std::endl;
				}
			}
			c2->Print(peds.c_str());
		}
	}
  c2->Print(peds_stop.c_str());//,"Portrait");

}

//-------------------------------------------------------------------
bool DmpAlgRdcAna::FinalizeDynodeRatioMode()
{	
	//TSpectrum *spectrum=new TSpectrum();
	//int npeak;
	//float *peakValue;
	//TF1 *fit;
	//double pedestalV,pedestalE;
  string dyr="dynoderatio.pdf";
  string dyr_start=dyr+"[";
  string dyr_stop=dyr+"]";
  TCanvas *c2 = new TCanvas("c2","dyRatio",0,0,600,800);
	c2->Divide(1,2);
  c2->Print(dyr_start.c_str());//,"Portrait");
	for(int layer=0;layer<14;layer++){
	  for(int barid=0;barid<22;barid++){
		  for(int side=0;side<2;side++){
				c2->cd(side+1);
				dynodeRatio[layer][barid][side]->Draw();
				/*
				npeak = spectrum->Search(histped[layer][barid][sidedy]);
				if(npeak>0){
				  peakValue = spectrum->GetPositionX();
					histped[layer][barid][sidedy]->Fit("gaus","","",peakValue[0]-50,peakValue[0]+50);
					fit = histped[layer][barid][sidedy]->GetFunction("gaus");
          pedestalV = fit->GetParameter(1);
          //pedestalE = fit->GetParError(1);
          pedestalE = fit->GetParameter(2);
		 		}
		 		else{
				  std::cout<<"peak number is "<<npeak<<std::endl;
				}
				*/
			}
			c2->Print(dyr.c_str());
		}
	}
  c2->Print(dyr_stop.c_str());//,"Portrait");

}

//-------------------------------------------------------------------
void DmpAlgRdcAna::Set(const std::string &type,const std::string &value){
  if(OptMap.find(type) == OptMap.end()){
    DmpLogError<<"No argument type: "<<type<<DmpLogEndl;
		DmpLogInfo<<"Possible option are "<<DmpLogEndl;
		for (std::map<string,short>::iterator anOpt= OptMap.begin();anOpt!=OptMap.end();anOpt++){
		  DmpLogInfo<<anOpt->first<<DmpLogEndl;
		} 
		throw;
	}
	
	switch(OptMap[type])
	{
	  case 0:
		  runMode = value;
			break;
		case 1:
		  adc_low=atoi(value.c_str());
			break;
		case 2:
		  adc_up=atoi(value.c_str());
			break;
	}
	if(value=="Dynode25"){
	  runMode="Dynode";
		dy1=2;
		dy2=5;
	}
	else if(value=="Dynode28"){
	  runMode="Dynode";
		dy1=2;
		dy2=8;
	}
  else if(value=="Dynode58"||value=="Dynode"){
	  runMode="Dynode";
		dy1=5;
		dy2=8;
	}

}

