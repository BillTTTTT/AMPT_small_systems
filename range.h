#ifndef RANGE_H
#define RANGE_H

#include <iostream>
#include <cmath>

using namespace std;

// FVTXS
// -3.0 < eta < -1.0
bool ifFVTXS(double eta)
{
	if (eta > -3.0 && eta < -1.0) return true;
	else return false;
}

// FVTXN
// 1.0 < eta < 3.0
bool ifFVTXN(double eta)
{
	if (eta > 1.0 && eta < 3.0) return true;
	else return false;
}

// FVTX
// 1.0 < eta < 3.0  or  -3.0 < eta < -1.0
bool ifFVTX(double eta)
{
	if (ifFVTXS(eta) || ifFVTXN(eta)) return true;
	else return false;
}

// BBCS
// -3.9 < eta < -3.1
bool ifBBCS(double eta)
{
	if (eta > -3.9 && eta < -3.1) return true;
	else return false;
}

// BBCN
// 3.1 < eta < 3.9
bool ifBBCN(double eta)
{
	if (eta > 3.1 && eta < 3.9) return true;
	else return false;
}

// CNT
// -0.35 < eta < 0.35
bool ifCNT(double eta)
{
	if (eta > -0.35 && eta < 0.35) return true;
	else return false;
}


#endif