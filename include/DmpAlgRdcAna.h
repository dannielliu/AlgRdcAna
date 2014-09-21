#ifndef DmpAlgRdcAna_H
#define DmpAlgRdcAna_H

#include "DmpVAlg.h"
#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"

class TH1D;
class TH2D;

class DmpAlgRdcAna : public DmpVAlg{
/*
 *  DmpAlgRdcAna
 *
 */
public:
  DmpAlgRdcAna();
  ~DmpAlgRdcAna();

  void Set(const std::string &type,const std::string &value);
  // if you need to set some options for your algorithm at run time. Overload Set()
  bool Initialize();
  bool ProcessThisEvent();    // only for algorithm
  bool Finalize();

private:
  unsigned long eventID;
	DmpEvtHeader* fHeader;
	DmpEvtBgoRaw* fBgo;
	std::string runMode; //can be mips, pedestal, dynode ratio

  TH1D *h1;
	long time_low,time_up;
	long time_cur;
	int adc_low,adc_up;
	TH1D *histped[14][24][6];
	short dy1,dy2;
	TH2D *dynodeRatio[14][22][2];

public:
  bool InitializePedMode();
  bool ProcessPedMode();
	bool FinalizePedMode();
	bool InitializeDynodeRatioMode();
	bool ProcessDynodeRatioMode(short dy1=5,short dy2=8);
  bool FinalizeDynodeRatioMode();
};

#endif
