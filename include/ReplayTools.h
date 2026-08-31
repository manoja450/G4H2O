#ifndef __ReplayTools_H__
#define __ReplayTools_H__

#include "TObject.h"
#include "TTimer.h"
#include "TStopwatch.h"
#include "TDatime.h"
#include "TObjArray.h"

#include <iostream>

class ReplayTools : public TObject
{    
public:
	ReplayTools();
	virtual ~ReplayTools() {;}
	
    void PrepareETFTimer(Double_t timeInterval, Int_t tEvents);
    void UpdateStatus(); 
    void StartUpdateTimer();
    void DisplayProgress(Int_t numEvent);
    void SetCurrentEvent(Int_t theEventNum) {currentEvent = theEventNum;}
    
private:
	
    TTimer *updateTimer;
    Int_t numPrimaryEvents;
    Int_t currentEvent;
    Double_t timeOutTime;
    
    void GetDHMS( Double_t timeInSecs, Char_t d[], Char_t h[], Char_t m[], Char_t s[]);

    TStopwatch *timer;
    struct Date_t { UInt_t month, day, year; };
    Date_t runDate;
    TDatime *dateTime;
    struct Time_t { UInt_t hour, minute, second; };
    Time_t currentTime;
    
};

#endif
