
#include "ReplayTools.h"
#include "TMath.h"
#include "TObjString.h"

#include <stdio.h>
#include <unistd.h>

using namespace std;


ReplayTools::ReplayTools(){

    updateTimer = 0;
    numPrimaryEvents = 0;
    currentEvent = 0;

    dateTime = new TDatime();
    timer = new TStopwatch();

}

void ReplayTools::PrepareETFTimer(Double_t timeInterval, Int_t tEvents){
 
    timeOutTime = timeInterval;
    numPrimaryEvents = tEvents;
    updateTimer = new TTimer(timeInterval,kTRUE);
    updateTimer->Connect("Timeout()", "ReplayTools", this, "UpdateStatus()");
    
}

void ReplayTools::StartUpdateTimer(){

    timer->Start();
    updateTimer->TurnOn();

}

void ReplayTools::UpdateStatus(){
    
    DisplayProgress(currentEvent);

}

void ReplayTools::DisplayProgress( Int_t numEvent )
{
    
    Char_t elapsedSecs[100];
    Char_t elapsedMins[100];
    Char_t elapsedHrs[100];
    Char_t elapsedDays[100];
    Char_t eRate[100];
    Char_t thisSec[100];
    Char_t thisMin[100];
    Char_t thisHr[100];
    Char_t finishSecs[100];
    Char_t finishMins[100];
    Char_t finishHrs[100];
    Char_t finishDays[100];
    Char_t month[100];
    Char_t day[100];
    Char_t year[100];
    
    Double_t eventRate;
    Double_t etf;
    Double_t totalRealTime = timer->RealTime();
    timer->Continue();
    
    if(numEvent == -1){
        eventRate = -0.01;
        etf = 0.0;
    }
    else if(numEvent == -2){
        eventRate = ((Double_t)numPrimaryEvents)/totalRealTime;
        etf = 0.0;
    }
    else{
        eventRate = ((Double_t)numEvent)/totalRealTime;
        etf = (Double_t)(numPrimaryEvents-numEvent)/eventRate;
    }
    
    GetDHMS(totalRealTime, elapsedDays,elapsedHrs,elapsedMins,elapsedSecs);
    GetDHMS(etf,           finishDays, finishHrs, finishMins, finishSecs);
    
    dateTime->Set();
    currentTime.hour   = dateTime->GetHour();
    currentTime.minute = dateTime->GetMinute();
    currentTime.second = dateTime->GetSecond();
    runDate.month      = dateTime->GetMonth();
    runDate.day        = dateTime->GetDay();
    runDate.year       = dateTime->GetYear();
    
    //
    if( runDate.month <= 9 ) sprintf(month, "0%d", runDate.month);
    else                     sprintf(month, "%d", runDate.month);
    if( runDate.day   <= 9 ) sprintf(day, "0%d", runDate.day);
    else                     sprintf(day, "%d", runDate.day);
    sprintf(year, "%d", runDate.year);
    
    //
    sprintf(eRate, "%8.2f", eventRate);
    
    //
    if( currentTime.second <= 9 ) sprintf(thisSec, "0%d", currentTime.second);
    else                          sprintf(thisSec, "%d", currentTime.second);
    if( currentTime.minute <= 9 ) sprintf(thisMin, "0%d", currentTime.minute);
    else                          sprintf(thisMin, "%d", currentTime.minute);
    if( currentTime.hour   <= 9 ) sprintf(thisHr, "0%d", currentTime.hour);
    else                          sprintf(thisHr, "%d", currentTime.hour);
    
    //
    if( numEvent==-1 )
    {
        cout << "\tRun Began On ";
        cout << month << "/";
        cout << day << "/";
        cout << year << "\n";
        cout << "\tNumber of Primary Events To Be Generated: ";
        cout << numPrimaryEvents << "\n";
        cout << endl;
        printf("\t Clock   -   Elapsed   -  To Finish  -      Event Number   -       Event Rate   \n");
        printf("\tHH:MM:SS - DD:HH:MM:SS - DD:HH:MM:SS - \n");
    }
    
    printf("\t%2s:%2s:%2s - %2s:%2s:%2s:%2s - %2s:%2s:%2s:%2s - ",thisHr,thisMin,thisSec,elapsedDays,elapsedHrs,elapsedMins,elapsedSecs,finishDays,finishHrs,finishMins,finishSecs);
    if(numEvent==-1) printf("%12d events",0);
    if(numEvent==-2) printf("%12d events",numPrimaryEvents);
    if(numEvent>0) printf("%12d events",numEvent);
    printf(" - ");
    if( totalRealTime <= 0 || eventRate < 0.0 ) printf("%8s events/sec","----");
    else printf("%8s events/sec",eRate);
    
    cout<<endl;
    
    if( numEvent==-2 )
    {
        cout << endl;
        cout << "\tRun Ended On ";
        cout << month << "/";
        cout << day << "/";
        cout << year << "\n";
        cout << "\tNumber of Primary Events Generated: ";
        cout << numPrimaryEvents << "\n\n";
    }
    
}//END of DisplayProgress()

void ReplayTools::GetDHMS( Double_t timeInSecs, Char_t d[], Char_t h[], Char_t m[], Char_t s[])
{
    Int_t t = TMath::Nint(timeInSecs);
    
    Int_t days = 0;
    Int_t hrs  = 0;
    Int_t mins = 0;
    Int_t secs = 0;
    
    Int_t remaining_s = 0;
    
    days = t/86400;
    remaining_s = t%86400;
    hrs = remaining_s/3600;
    remaining_s = remaining_s%3600;
    mins = remaining_s/60;
    remaining_s = remaining_s%60;
    secs = remaining_s;
    
    if( days <=9 ) sprintf(d,"0%d",days);
    else           sprintf(d, "%d",days);
    if( hrs  <=9 ) sprintf(h,"0%d",hrs );
    else           sprintf(h, "%d",hrs );
    if( mins <=9 ) sprintf(m,"0%d",mins);
    else           sprintf(m, "%d",mins);
    if( secs <=9 ) sprintf(s,"0%d",secs);
    else           sprintf(s, "%d",secs);
    
}//END of GetDHMS()

