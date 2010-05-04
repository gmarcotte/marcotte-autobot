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
const double neutralPct = 0.1;

const double redPct = 0.6;
const double redRwd = 8.0;
const double redPrb = 0.3;

const double bluePct = 0.3;
const double blueRwd = 3.0;
const double bluePrb = 0.8;

// Forager parameters
const double foragerSpeed = 1.0;
const double visualConeAngle = 0.1745;  //radians, 10 degrees
const double projectionRadius = 1.0;
const double projection_dR = 1.0 / 200.0;
const double projection_dTheta = 2*3.1415926 / 200.0;

#endif // CONFIG_H