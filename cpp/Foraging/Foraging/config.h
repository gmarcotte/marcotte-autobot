#ifndef CONFIG_H
#define CONFIG_H

#include <windows.h>

// Window Setup
const char* g_szWindowClassName = "forager";
const char* g_szApplicationName = "Foraging Simulation";
const UINT nWindow_Width = 600;
const UINT nWindow_Height = 600;

// Simulation Parameters
const double sim_dT = 1.0;

// Field parameters
const UINT fieldSide = 60;
const double neutralPct = 0.2;

const double redPct = 0.4;
const double redRwd = 0.7;
const double redPrb = 1.0;

const double bluePct = 0.4;
const double blueRwd = 1.0;
const double bluePrb = 0.2;

// Forager parameters
const double foragerSpeed = 1.5;
const double visualConeAngle = 0.1745;  //radians, 10 degrees
const double projectionRadius = 1.0;
const double projection_dR = 1.0 / 50.0;
const double projection_dTheta = 2*3.1415926 / 50.0;

#endif // CONFIG_H