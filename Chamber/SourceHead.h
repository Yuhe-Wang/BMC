#pragma once
#ifndef _SOURCEHEAD_H_
#define _SOURCEHEAD_H_

#include "../Tools/Tools.h"
// interface functions provided by this Source Head module

bool SourceHead_Init(ConfigFile* config, int NT);

void SourceHead_Delete();

int SourceHead_Sample(int it, double& E, double& weight, MonteCarlo::Vector& x, MonteCarlo::Vector& v, ParticleType& type);

double SourceHead_Emax();

double SourceHead_NormFactor(); // return the real history number at given time

#endif