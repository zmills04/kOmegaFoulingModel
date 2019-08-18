// clDisplay.cpp: implementation of the clDisplay class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "clDisplay.h"
#include "clVariablesLS.h"
#include "clVariablesLB.h"
#include "clVariablesFD.h"
#include "clVariablesTR.h"
#include "clVariablesFL.h"
#include "oclEnvironment.h"

void clDisplay::finish()
{
	ShowTime();
	printf("\n");
}

void clDisplay::ini()
{
	counter = 0;
	time(&StartTimeSys);
	strcpy(cStartTimeSys,&(ctime(&StartTimeSys)[11]));
};

void clDisplay::ShowTime()
{
	time_t CurrTimeSys, ElapsTimeSys, FinishTimeSys, RemTimeSys;
	char cFinishTimeSys[26];

	time(&CurrTimeSys);

	ElapsTimeSys = CurrTimeSys - StartTimeSys;
	FinishTimeSys = StartTimeSys + (time_t)((double)ElapsTimeSys * (p.StopTime - p.RestartTime) / (p.TimeN - p.RestartTime));
	int ElapsH = 0, ElapsM = 0, ElapsS = 0;
	if (ElapsTimeSys > 0)
	{
		ElapsH = (int)(ElapsTimeSys / (60 * 60));
		ElapsM = (int)((ElapsTimeSys - ElapsH * 60 * 60) / 60);
		ElapsS = (int)((ElapsTimeSys - ElapsH * 60 * 60 - ElapsM * 60));
	}
	int RemH = -1, RemM = 0, RemS = 0;
	if (FinishTimeSys > 0)
	{
		RemTimeSys = FinishTimeSys - CurrTimeSys;
		RemH = (int)(RemTimeSys / (60 * 60));
		RemM = (int)((RemTimeSys - RemH * 60 * 60) / 60);
		RemS = (int)((RemTimeSys - RemH * 60 * 60 - RemM * 60));
	}

	if (FinishTimeSys > 0)
	{
		strcpy(cFinishTimeSys, &(ctime(&FinishTimeSys)[11]));
	}
	else
	{
		sprintf(cFinishTimeSys, "00:00:00");
	}
	
	double usum_ = vlb.sumUx.reduceSingle();
	double rhosum_ = vlb.sumRo.reduceSingle();
	double tempsum_ = vfd.sumTemp.reduceSingle();

	printf("T%f U%f Ro%g E%d:%.2d:%.2d R%d:%.2d:%.2d P%5.2f%%\r",
		tempsum_ / rhosum_, usum_/rhosum_, rhosum_, ElapsH, ElapsM, ElapsS, RemH, RemM, RemS, 100. * p.Time / p.StopTime);

};

void clDisplay::step()
{
	counter++;
	if (++counter < p.displaySignalFreq)
		return;

	ShowTime();
	counter = 0;
}




