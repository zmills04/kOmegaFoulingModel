// clDisplay.h: interface for the clDisplay class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLDISPLAY_H__904010C9_8CE3_405E_9266_71D4D0D77458__INCLUDED_)
#define AFX_CLDISPLAY_H__904010C9_8CE3_405E_9266_71D4D0D77458__INCLUDED_


//TODO: Clean this up and implement openGL calls in this cpp/h file


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "clProblem.h"
class clDisplay  
{
public:
	time_t StartTimeSys;
	char cStartTimeSys[26];
	unsigned int counter;


	clDisplay(){ counter = 0; };
	virtual ~clDisplay(){};

	void finish();
	
	void ini(void);      

	void ShowTime(void);

	void start(){ShowTime();}

	void step();
};

#endif // !defined(AFX_CLDISPLAY_H__904010C9_8CE3_405E_9266_71D4D0D77458__INCLUDED_)
