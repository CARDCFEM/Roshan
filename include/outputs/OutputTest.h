
#pragma once

#include "Console.h"

class OutputTest : public Console
{
public:
	OutputTest(const InputParameters & parameters);
	virtual ~OutputTest(){};

	void initialSetup();
	void outputSystemInformation();
};

template<>
InputParameters validParams<OutputTest> ();
