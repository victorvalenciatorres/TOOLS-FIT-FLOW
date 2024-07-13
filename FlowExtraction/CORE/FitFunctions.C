// Author: Victor Valencia
#include <vector>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include "TF1.h"

Bool_t reject;

// enumerators used to select fit functions and tails parameters
enum Efunction
{
	kVWG = 100,
	kVWGQuadratic = 200,
	kPol2Pol3 = 300,
	kDoubleExp = 400,
	kExp = 500,
	kv2Background = 600,
	kv2BackgroundPol3 = 700,
	kv2BackgroundPolExp = 800,
	kCB = 10,
	kCBExtended = 20,
	kNA60 = 30,
	kv2Pol2vsMassCB2VWG = 40,
	kv2Pol2vsMassNA60QVWG = 45,
	kv2Pol2vsMassCB2Pol2Pol3 = 50,
	kv2Pol2vsMassCB2QVWG = 60,
	kv2Pol2vsMassNA60VWG = 55,
	kv2Pol2vsMassNA60Pol2Pol3 = 65,
	kv2Pol3vsMassCB2VWG = 22,
	kv2Pol3vsMassNA60QVWG = 32,
	kv2Pol3vsMassCB2Pol2Pol3 = 42,
	kv2Pol3vsMassCB2QVWG = 52,
	kv2Pol3vsMassNA60VWG = 62,
	kv2Pol3vsMassNA60Pol2Pol3 = 72,
	kv2PolExpvsMassCB2VWG = 21,
	kv2PolExpvsMassNA60QVWG = 31,
	kv2PolExpvsMassCB2Pol2Pol3 = 41,
	kv2PolExpvsMassCB2QVWG = 51,
	kv2PolExpvsMassNA60VWG = 61,
	kv2PolExpvsMassNA60Pol2Pol3 = 71

};

enum Etails
{
	kEMB = 0,
	kPP = 1,
	kSTARLIGHTcoh = 2,
	kSTARLIGHTincoh = 3
};

enum Epart
{
	kJPsi = 0,
	kPsi2S = 1
};

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
Double_t VWG(Double_t *x, Double_t *par) // Variable width Gaussian
{
	// par[0]=1/(sigma*sqrt(2*PI)  Normalization
	// par[1]=mu  Mean
	// par[2]=alpha
	// par[3]=beta

	Double_t sigma = par[2] + par[3] * ((x[0] - par[1]) / par[1]);
	return par[0] * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / (2. * sigma * sigma));
}

Double_t Pol3(Double_t *x, Double_t *par)
{
	return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
}

//___________________________________________________________________________________________________________
Double_t VWGQuadratic(Double_t *x, Double_t *par) // Quadratic Variable width Gaussian
{
	// par[0]=1/(sigma*sqrt(2*PI)  Normalization
	// par[1]=mu  Mean
	// par[2]=alpha
	// par[3]=beta
	// par[4]=gamma

	if (reject && x[0] > 2.8 && x[0] < 3.8)
	{
		TF1::RejectPoint();
		return 0;
	}

	Double_t sigma = par[2] + par[3] * ((x[0] - par[1]) / par[1]) + par[4] * ((x[0] - par[1]) / par[1]) * ((x[0] - par[1]) / par[1]);
	return par[0] * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / (2. * sigma * sigma));
}

//___________________________________________________________________________________________________________
Double_t Pol2Pol3(Double_t *x, Double_t *par) // Ratio of polynomials
{
	// par[0]=Normalization
	// par[1]=a1
	// par[2]=a2
	// par[3]=b1
	// par[4]=b2
	// par[5]=b3

	if (reject && x[0] > 2.8 && x[0] < 3.8)
	{
		TF1::RejectPoint();
		return 0;
	}

	return par[0] * (1 + par[1] * x[0] + par[2] * x[0] * x[0]) / (1 + par[3] * x[0] + par[4] * x[0] * x[0] + par[5] * x[0] * x[0] * x[0]);
}

//___________________________________________________________________________________________________________
Double_t DoubleExp(Double_t *x, Double_t *par) // sum of exponentials
{
	// par[0]=Norm1
	// par[1]=alpha1
	// par[2]=Norm2
	// par[3]=alpha2

	if (reject && x[0] > 2.8 && x[0] < 3.8)
	{
		TF1::RejectPoint();
		return 0;
	}

	return par[0] * TMath::Exp(par[1] * x[0]) + par[2] * TMath::Exp(par[3] * x[0]);
}

//___________________________________________________________________________________________________________
Double_t Exp(Double_t *x, Double_t *par) // exponential
{
	// par[0]=Norm
	// par[1]=alpha

	if (reject && x[0] > 2.8 && x[0] < 3.8)
	{
		TF1::RejectPoint();
		return 0;
	}

	return par[0] * TMath::Exp(par[1] * x[0]);
}

//___________________________________________________________________________________________________________
Double_t CrystalBall(Double_t *x, Double_t *par)
{
	// par[0]=N Normalization
	// par[1]=mu Mean
	// par[2]=sigma Width
	// par[3]=alphaL alpha for left tail
	// par[4]=nL for left tail

	Double_t t = (x[0] - par[1]) / par[2];
	Double_t absAlphaL = TMath::Abs(par[3]);

	if (par[2] < 0)
		t = -t;

	if (t > -absAlphaL)
	{
		return par[0] * (TMath::Exp(-0.5 * t * t));
	}

	if (t <= -absAlphaL)
	{
		Double_t A = TMath::Power(par[4] / absAlphaL, par[4]) * TMath::Exp(-0.5 * absAlphaL * absAlphaL);
		Double_t B = par[4] / absAlphaL - absAlphaL;

		return par[0] * A * TMath::Power(B - t, -par[4]);
	}

	return 0.;
}

//___________________________________________________________________________________________________________
Double_t CrystalBallExtended(Double_t *x, Double_t *par)
{
	// par[0]=N  Normalization
	// par[1]=mu  Mean
	// par[2]=sigma  Width
	// par[3]=alphaL  Alpha for the left tail
	// par[4]=nL for the left tail
	// par[5]=alphaR  Alpha for the right tail
	// par[6]=nR for the right tail

	Double_t t = (x[0] - par[1]) / par[2];
	if (par[2] < 0)
		t = -t;

	Double_t absAlphaL = TMath::Abs(par[3]);
	Double_t absAlphaR = TMath::Abs(par[5]);

	if (t > -absAlphaL && t < absAlphaR) // gaussian core
	{
		return par[0] * (TMath::Exp(-0.5 * t * t));
	}

	if (t <= -absAlphaL) // left tail
	{
		Double_t A = TMath::Power(par[4] / absAlphaL, par[4]) * TMath::Exp(-0.5 * absAlphaL * absAlphaL);
		Double_t B = par[4] / absAlphaL - absAlphaL;

		return par[0] * (A / TMath::Power(B - t, par[4]));
	}

	if (t >= absAlphaR) // right tail
	{
		Double_t C = TMath::Power(par[6] / absAlphaR, par[6]) * TMath::Exp(-0.5 * absAlphaR * absAlphaR);
		Double_t D = par[6] / absAlphaR - absAlphaR;

		return par[0] * (C / TMath::Power(D + t, par[6]));
	}

	return 0.;
}

//___________________________________________________________________________________________________________
Double_t DoubleCrystalBallExtended(Double_t *x, Double_t *par)
{
	// par[0]=N  Normalization
	// par[1]=mu  Mean
	// par[2]=sigma  Width
	// par[3]=alphaL  Alpha for the left tail
	// par[4]=nL for the left tail
	// par[5]=alphaR  Alpha for the right tail
	// par[6]=nR for the right tail
	// par[7]=Npsi'

	Double_t absAlphaL = TMath::Abs(par[3]);
	Double_t absAlphaR = TMath::Abs(par[5]);
	Int_t signAlphaL = TMath::Sign(1, par[3]);
	Int_t signAlphaR = TMath::Sign(1, par[5]);

	Double_t t1 = signAlphaL * signAlphaR * (x[0] - par[1]) / par[2];
	// if (par[2] < 0) t1 = -t1;

	Double_t jpsi = 0.0;
	Double_t psi2S = 0.0;
	if (t1 > -absAlphaL && t1 < absAlphaR) // gaussian core
	{
		jpsi = par[0] * (TMath::Exp(-0.5 * t1 * t1));
	}

	if (t1 <= -absAlphaL) // left tail
	{
		Double_t A = TMath::Power(par[4] / absAlphaL, par[4]) * TMath::Exp(-0.5 * absAlphaL * absAlphaL);
		Double_t B = par[4] / absAlphaL - absAlphaL;

		jpsi = par[0] * (A / TMath::Power(B - t1, par[4]));
	}

	if (t1 >= absAlphaR) // right tail
	{
		Double_t C = TMath::Power(par[6] / absAlphaR, par[6]) * TMath::Exp(-0.5 * absAlphaR * absAlphaR);
		Double_t D = par[6] / absAlphaR - absAlphaR;

		jpsi = par[0] * (C / TMath::Power(D + t1, par[6]));
	}

	// psi2S
	Double_t t2 = signAlphaL * signAlphaR * (x[0] - (par[1] + 3.686109 - 3.096916)) / (par[2] * 1.05);
	if (t2 > -absAlphaL && t2 < absAlphaR) // gaussian core
	{
		psi2S = par[7] * (TMath::Exp(-0.5 * t2 * t2));
	}

	if (t2 <= -absAlphaL) // left tail
	{
		Double_t A = TMath::Power(par[4] / absAlphaL, par[4]) * TMath::Exp(-0.5 * absAlphaL * absAlphaL);
		Double_t B = par[4] / absAlphaL - absAlphaL;

		psi2S = par[7] * (A / TMath::Power(B - t2, par[4]));
	}

	if (t2 >= absAlphaR) // right tail
	{
		Double_t C = TMath::Power(par[6] / absAlphaR, par[6]) * TMath::Exp(-0.5 * absAlphaR * absAlphaR);
		Double_t D = par[6] / absAlphaR - absAlphaR;

		psi2S = par[7] * (C / TMath::Power(D + t2, par[6]));
	}
	return jpsi + psi2S;
}

//___________________________________________________________________________________________________________
Double_t NA60(Double_t *x, Double_t *par)
{
	// par[0]=N  Normalization
	// par[1]=mu  Mean
	// par[2]=sigma  Width
	// par[3]=alphaL  Alpha for the left tail
	// par[4]=p1 for the left tail
	// par[5]=p2 for the left tail
	// par[6]=p3 for the left tail
	// par[7]=alphaR  Alpha for the right tail
	// par[8]=p1 for the right tail
	// par[9]=p2 for the right tail
	// par[10]=p3 for the right tail

	Double_t t = (x[0] - par[1]) / par[2];
	if (par[2] < 0)
		t = -t;
	Double_t t0 = 0.0;

	if (t < par[3])
	{
		t0 = 1 + TMath::Power(par[4] * (par[3] - t), par[5] - par[6] * TMath::Sqrt(par[3] - t));
	}

	if (t >= par[3] && t < par[7])
	{
		t0 = 1;
	}

	if (t >= par[7])
	{
		t0 = 1 + TMath::Power(par[8] * (t - par[7]), par[9] - par[10] * TMath::Sqrt(t - par[7]));
	}

	return par[0] * TMath::Exp(-0.5 * (t / t0) * (t / t0));
}

//___________________________________________________________________________________________________________
Double_t DoubleNA60(Double_t *x, Double_t *par)
{
	// par[0]=N  Normalization
	// par[1]=mu  Mean
	// par[2]=sigma  Width
	// par[3]=alphaL  Alpha for the left tail
	// par[4]=p1 for the left tail
	// par[5]=p2 for the left tail
	// par[6]=p3 for the left tail
	// par[7]=alphaR  Alpha for the right tail
	// par[8]=p1 for the right tail
	// par[9]=p2 for the right tail
	// par[10]=p3 for the right tail
	// par[11]=N psi2s

	Double_t t1 = (x[0] - par[1]) / par[2];
	if (par[2] < 0)
		t1 = -t1;
	Double_t tjpsi = 0.0;
	Double_t tpsi2S = 0.0;
	Double_t jpsi = 0.0;
	Double_t psi2S = 0.0;

	if (t1 < par[3])
	{
		tjpsi = 1 + TMath::Power(par[4] * (par[3] - t1), par[5] - par[6] * TMath::Sqrt(par[3] - t1));
	}

	if (t1 >= par[3] && t1 < par[7])
	{
		tjpsi = 1;
	}

	if (t1 >= par[7])
	{
		tjpsi = 1 + TMath::Power(par[8] * (t1 - par[7]), par[9] - par[10] * TMath::Sqrt(t1 - par[7]));
	}

	// psi2S

	Double_t t2 = (x[0] - (par[1] + 3.686109 - 3.096916)) / (par[2] * 1.05);
	if (par[2] < 0)
		t2 = -t2;

	if (t2 < par[3])
	{
		tpsi2S = 1 + TMath::Power(par[4] * (par[3] - t2), par[5] - par[6] * TMath::Sqrt(par[3] - t2));
	}

	if (t2 >= par[3] && t2 < par[7])
	{
		tpsi2S = 1;
	}

	if (t2 >= par[7])
	{
		tpsi2S = 1 + TMath::Power(par[8] * (t2 - par[7]), par[9] - par[10] * TMath::Sqrt(t2 - par[7]));
	}

	return (par[0] * TMath::Exp(-0.5 * (t1 / tjpsi) * (t1 / tjpsi))) + (par[11] * TMath::Exp(-0.5 * (t2 / tpsi2S) * (t2 / tpsi2S)));
}

//___________________________________________________________________________________________________________
Double_t VWG_CBext(Double_t *x, Double_t *par) // 4 parmeters for BG and 7 parameters for signal
{
	return VWG(x, par) + CrystalBallExtended(x, &(par[4]));
}

//___________________________________________________________________________________________________________
Double_t VWGquad_CBext(Double_t *x, Double_t *par) // 5 parmeters for BG and 7 parameters for signal
{
	return VWGQuadratic(x, par) + CrystalBallExtended(x, &(par[5]));
}

//___________________________________________________________________________________________________________
Double_t VWGquad_DoubleCBext(Double_t *x, Double_t *par) // 5 parmeters for BG and 7+1 parameters for signal
{
	// BackGround
	// par[0]=1/(sigma*sqrt(2*PI)  Normalization
	// par[1]=mu  Mean
	// par[2]=alpha
	// par[3]=beta
	// par[4]=gamma

	// J/psi signal
	// par[5]=N  Normalization
	// par[6]=mu  Mean
	// par[7]=sigma  Width
	// par[8]=alphaL  Alpha for the left tail
	// par[9]=nL  n for the left tail
	// par[10]=alphaR  Alpha for the right tail
	// par[11]=nR  n for the right tail

	// psi' signal
	// par[12]=N  Normalization

	return VWGQuadratic(x, par) + DoubleCrystalBallExtended(x, &(par[5]));
}

//___________________________________________________________________________________________________________
Double_t VWGquad_DoubleNA60(Double_t *x, Double_t *par) // 5 parmeters for BG and 11+1 parameters for signal
{
	// BackGround
	// par[0]=1/(sigma*sqrt(2*PI)  Normalization
	// par[1]=mu  Mean
	// par[2]=alpha
	// par[3]=beta
	// par[4]=gamma

	// J/psi signal
	// par[5]=N  Normalization
	// par[6]=mu  Mean
	// par[7]=sigma  Width
	// par[8]=alphaL  Alpha for the left tail
	// par[9]=p1 for the left tail
	// par[10]=p2 for the left tail
	// par[11]=p3 for the left tail
	// par[12]=alphaR  Alpha for the right tail
	// par[13]=p1 for the right tail
	// par[14]=p2 for the right tail
	// par[15]=p3 for the right tail

	// psi' signal
	// par[16]=N  Normalization

	return VWGQuadratic(x, par) + DoubleNA60(x, &(par[5]));
}

Double_t VWGquad_NA60(Double_t *x, Double_t *par) // 5 parmeters for BG and 11+1 parameters for signal
{
	// BackGround
	// par[0]=1/(sigma*sqrt(2*PI)  Normalization
	// par[1]=mu  Mean
	// par[2]=alpha
	// par[3]=beta
	// par[4]=gamma

	// J/psi signal
	// par[5]=N  Normalization
	// par[6]=mu  Mean
	// par[7]=sigma  Width
	// par[8]=alphaL  Alpha for the left tail
	// par[9]=p1 for the left tail
	// par[10]=p2 for the left tail
	// par[11]=p3 for the left tail
	// par[12]=alphaR  Alpha for the right tail
	// par[13]=p1 for the right tail
	// par[14]=p2 for the right tail
	// par[15]=p3 for the right tail

	// psi' signal
	// par[16]=N  Normalization

	return VWGQuadratic(x, par) + NA60(x, &(par[5]));
}

Double_t VWG_NA60(Double_t *x, Double_t *par) // 4 parmeters for BG and 11 parameters for signal
{
	// BackGround
	// par[0]=1/(sigma*sqrt(2*PI)  Normalization
	// par[1]=mu  Mean
	// par[2]=alpha
	// par[3]=beta
	// par[4]=gamma

	// J/psi signal
	// par[5]=N  Normalization
	// par[6]=mu  Mean
	// par[7]=sigma  Width
	// par[8]=alphaL  Alpha for the left tail
	// par[9]=p1 for the left tail
	// par[10]=p2 for the left tail
	// par[11]=p3 for the left tail
	// par[12]=alphaR  Alpha for the right tail
	// par[13]=p1 for the right tail
	// par[14]=p2 for the right tail
	// par[15]=p3 for the right tail

	// psi' signal
	// par[16]=N  Normalization

	return VWG(x, par) + NA60(x, &(par[4]));
}

//___________________________________________________________________________________________________________
Double_t Pol2Pol3_DoubleCBext(Double_t *x, Double_t *par) // 6 parmeters for BG and 7+1 parameters for signal
{
	// BackGround
	// par[0]=Normalization
	// par[1]=a1
	// par[2]=a2
	// par[3]=b1
	// par[4]=b2
	// par[5]=b3

	// J/psi signal
	// par[6]=N  Normalization
	// par[7]=mu  Mean
	// par[8]=sigma  Width
	// par[9]=alphaL  Alpha for the left tail
	// par[10]=nL  n for the left tail
	// par[11]=alphaR  Alpha for the right tail
	// par[12]=nR  n for the right tail

	// psi' signal
	// par[13]=N  Normalization

	return Pol2Pol3(x, par) + DoubleCrystalBallExtended(x, &(par[6]));
}

//___________________________________________________________________________________________________________
Double_t Pol2Pol3_DoubleNA60(Double_t *x, Double_t *par) // 6 parmeters for BG and 11+1 parameters for signal
{
	// BackGround
	// par[0]=Normalization
	// par[1]=a1
	// par[2]=a2
	// par[3]=b1
	// par[4]=b2
	// par[5]=b3

	// J/psi signal
	// par[6]=N  Normalization
	// par[7]=mu  Mean
	// par[8]=sigma  Width
	// par[9]=alphaL  Alpha for the left tail
	// par[10]=p1 for the left tail
	// par[11]=p2 for the left tail
	// par[12]=p3 for the left tail
	// par[13]=alphaR  Alpha for the right tail
	// par[14]=p1 for the right tail
	// par[15]=p2 for the right tail
	// par[16]=p3 for the right tail

	// psi' signal
	// par[17]=N  Normalization

	return Pol2Pol3(x, par) + DoubleNA60(x, &(par[6]));
}

Double_t Pol2Pol3_NA60(Double_t *x, Double_t *par) // 6 parmeters for BG and 11+1 parameters for signal
{
	// BackGround
	// par[0]=Normalization
	// par[1]=a1
	// par[2]=a2
	// par[3]=b1
	// par[4]=b2
	// par[5]=b3

	// J/psi signal
	// par[6]=N  Normalization
	// par[7]=mu  Mean
	// par[8]=sigma  Width
	// par[9]=alphaL  Alpha for the left tail
	// par[10]=p1 for the left tail
	// par[11]=p2 for the left tail
	// par[12]=p3 for the left tail
	// par[13]=alphaR  Alpha for the right tail
	// par[14]=p1 for the right tail
	// par[15]=p2 for the right tail
	// par[16]=p3 for the right tail

	// psi' signal
	// par[17]=N  Normalization

	return Pol2Pol3(x, par) + NA60(x, &(par[6]));
}

//___________________________________________________________________________________________________________
Double_t DoubleExp_DoubleCBext(Double_t *x, Double_t *par) // 4 parameters for BG ans 7+1 for signals
{
	// BackGround
	// par[0]=Norm1
	// par[1]=alpha1
	// par[2]=Norm2
	// par[3]=alpha2

	// J/psi signal
	// par[4]=N  Normalization
	// par[5]=mu  Mean
	// par[6]=sigma  Width
	// par[7]=alphaL  Alpha for the left tail
	// par[8]=nL  n for the left tail
	// par[9]=alphaR  Alpha for the right tail
	// par[10]=nR  n for the right tail

	// psi' signal
	// par[11]=N  Normalization

	return DoubleExp(x, par) + DoubleCrystalBallExtended(x, &(par[4]));
}

//___________________________________________________________________________________________________________
Double_t DoubleExp_DoubleNA60(Double_t *x, Double_t *par) // 4 parameters for BG ans 11+1 for signals
{
	// BackGround
	// par[0]=Norm1
	// par[1]=alpha1
	// par[2]=Norm2
	// par[3]=alpha2

	// J/psi signal
	// par[4]=N  Normalization
	// par[5]=mu  Mean
	// par[6]=sigma  Width
	// par[7]=alphaL  Alpha for the left tail
	// par[8]=p1 for the left tail
	// par[9]=p2 for the left tail
	// par[10]=p3 for the left tail
	// par[11]=alphaR  Alpha for the right tail
	// par[12]=p1 for the right tail
	// par[13]=p2 for the right tail
	// par[14]=p3 for the right tail

	// psi' signal
	// par[15]=N  Normalization

	return DoubleExp(x, par) + DoubleNA60(x, &(par[4]));
}
//___________________________________________________________________________________________________________
Double_t Exp_DoubleCBext(Double_t *x, Double_t *par) // 4 parameters for BG ans 7+1 for signals
{
	// BackGround
	// par[0]=Norm
	// par[1]=alpha

	// J/psi signal
	// par[2]=N  Normalization
	// par[3]=mu  Mean
	// par[4]=sigma  Width
	// par[5]=alphaL  Alpha for the left tail
	// par[6]=nL  n for the left tail
	// par[7]=alphaR  Alpha for the right tail
	// par[8]=nR  n for the right tail

	// psi' signal
	// par[9]=N  Normalization

	return DoubleExp(x, par) + DoubleCrystalBallExtended(x, &(par[2]));
}

//___________________________________________________________________________________________________________
Double_t Exp_DoubleNA60(Double_t *x, Double_t *par) // 4 parameters for BG ans 11+1 for signals
{
	// BackGround
	// par[0]=Norm
	// par[1]=alpha

	// J/psi signal
	// par[2]=N  Normalization
	// par[3]=mu  Mean
	// par[4]=sigma  Width
	// par[5]=alphaL  Alpha for the left tail
	// par[6]=p1 for the left tail
	// par[7]=p2 for the left tail
	// par[8]=p3 for the left tail
	// par[9]=alphaR  Alpha for the right tail
	// par[10]=p1 for the right tail
	// par[11]=p2 for the right tail
	// par[12]=p3 for the right tail

	// psi' signal
	// par[13]=N  Normalization

	return DoubleExp(x, par) + DoubleNA60(x, &(par[2]));
}

Double_t v2Background(Double_t *x, Double_t *par)
{
	return par[0] * ((x[0] - 3.1) * (x[0] - 3.1)) + par[1] * (x[0] - 3.1) + par[2];
}

Double_t v2BackgroundPol2(Double_t *x, Double_t *par)
{
	return par[0] + par[1] * (x[0] - 3.1) + par[2] * ((x[0] - 3.1) * (x[0] - 3.1));
}

// v2 background functions pol3
//_______________________________________________________________________________
Double_t v2BackgroundPol3(Double_t *x, Double_t *par)
{
	return par[0] + par[1] * (x[0] - 3.1) + par[2] * ((x[0] - 3.1) * (x[0] - 3.1)) + par[3] * ((x[0] - 3.1) * (x[0] - 3.1) * (x[0] - 3.1));
}

Double_t v2Bkg(Double_t *x, Double_t *par)
{ // bkg 4 + sig 7 + flow bkg 3
	return (VWG(x, par) * v2Background(x, &par[11])) / (VWG(x, par) + CrystalBallExtended(x, &par[4]));
}

Double_t v2BackgroundCheb(Double_t *x, Double_t *par)
{
	Double_t xmin = 2.2;
	Double_t xmax = 4.7;
	// Double_t xmin = 8;
	// Double_t xmax = 11;
	double xx = (2.0 * x[0] - xmin - xmax) / (xmax - xmin);
	const int order = 4;
	Double_t fT[order + 1] = {0};
	if (order == 1)
		return par[0];
	if (order == 2)
		return par[0] + xx * par[1];
	// build the polynomials
	fT[0] = 1;
	fT[1] = xx;
	for (int i = 1; i < order; ++i)
	{
		fT[i + 1] = 2 * xx * fT[i] - fT[i - 1];
	}
	double sum = par[0] * fT[0];
	for (int i = 1; i <= order; ++i)
	{
		sum += par[i] * fT[i];
	}
	if (reject && x[0] > 9.1 && x[0] < 9.8)
	{
		TF1::RejectPoint();
		return 0;
	}
	return sum;
}

Double_t v2BackgroundPolExp(Double_t *x, Double_t *par)
{
	//	return (par[2] + par[3]*(x[0]-3.1) + par[4]*(x[0]-3.1)*(x[0]-3.1))*exp(par[0]+par[1]*(x[0]-3.1)); //par[0]*x[0] + par[1];
	return par[0] * (1. / exp(par[1]) + par[3] * (x[0] - 3.1) + par[4] * (x[0] - 3.1) * (x[0] - 3.1)) * exp(par[1] + par[2] * (x[0] - 3.1)); // par[0]*x[0] + par[1];
}

////////////////////////////////////////////////////////////////////////////////////

Double_t v2Pol2vsMassCB2VWG(Double_t *x, Double_t *par)
{ // bkg 4 + sig 7 + flow bkg 3 + flow sig 1= 15
	return (VWG(x, par) * v2Background(x, &par[11]) + par[14] * CrystalBallExtended(x, &par[4])) / (VWG(x, par) + CrystalBallExtended(x, &par[4]));
}

Double_t v2Pol2vsMassNA60QVWG(Double_t *x, Double_t *par)
{ // bkg 5 + sig 11 + flow bkg 3 + flow sig 1 = 20
	return (VWGQuadratic(x, par) * v2Background(x, &par[16]) + par[19] * NA60(x, &par[5])) / (VWGQuadratic(x, par) + NA60(x, &par[5]));
}

Double_t v2Pol2vsMassNA60Pol2Pol3(Double_t *x, Double_t *par)
{ // bkg 6 + sig 11 + flow bkg 3 + flow sig 1 = 21
	return (Pol2Pol3(x, par) * v2Background(x, &par[17]) + par[20] * NA60(x, &par[6])) / (Pol2Pol3(x, par) + NA60(x, &par[6]));
}

Double_t v2Pol2vsMassNA60VWG(Double_t *x, Double_t *par)
{ // bkg 4 + sig 11 + flow bkg 3 + flow sig 1 = 19
	return (VWG(x, par) * v2Background(x, &par[15]) + par[18] * NA60(x, &par[4])) / (VWG(x, par) + NA60(x, &par[4]));
}

Double_t v2Pol2vsMassCB2Pol2Pol3(Double_t *x, Double_t *par)
{ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
	return (Pol2Pol3(x, par) * v2Background(x, &par[13]) + par[16] * CrystalBallExtended(x, &par[6])) / (Pol2Pol3(x, par) + CrystalBallExtended(x, &par[6]));
}

Double_t v2Pol2vsMassCB2QVWG(Double_t *x, Double_t *par)
{ // bkg 5 + sig 7 + flow bkg 3 + flow sig 1 = 16
	return (VWGQuadratic(x, par) * v2Background(x, &par[12]) + par[15] * CrystalBallExtended(x, &par[5])) / (VWGQuadratic(x, par) + CrystalBallExtended(x, &par[5]));
}

Double_t v2Pol3vsMassCB2VWG(Double_t *x, Double_t *par)
{ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1= 16

	return (VWG(x, par) * v2BackgroundPol3(x, &par[11]) + par[15] * CrystalBallExtended(x, &par[4])) / (VWG(x, par) + CrystalBallExtended(x, &par[4]));
}

Double_t v2Pol3vsMassNA60QVWG(Double_t *x, Double_t *par)
{ // bkg 5 + sig 11 + flow bkg 4 + flow sig 1 = 21
	return (VWGQuadratic(x, par) * v2BackgroundPol3(x, &par[16]) + par[20] * NA60(x, &par[5])) / (VWGQuadratic(x, par) + NA60(x, &par[5]));
}

Double_t v2Pol3vsMassNA60Pol2Pol3(Double_t *x, Double_t *par)
{ // bkg 6 + sig 11 + flow bkg 4 + flow sig 1 = 22
	return (Pol2Pol3(x, par) * v2BackgroundPol3(x, &par[17]) + par[21] * NA60(x, &par[6])) / (Pol2Pol3(x, par) + NA60(x, &par[6]));
}

Double_t v2Pol3vsMassCB2Pol2Pol3(Double_t *x, Double_t *par)
{ // bkg 6 + sig 7 + flow bkg 4 + flow sig 1 = 18
	return (Pol2Pol3(x, par) * v2BackgroundPol3(x, &par[13]) + par[17] * CrystalBallExtended(x, &par[6])) / (Pol2Pol3(x, par) + CrystalBallExtended(x, &par[6]));
}

Double_t v2Pol3vsMassCB2QVWG(Double_t *x, Double_t *par)
{ // bkg 5 + sig 7 + flow bkg 4 + flow sig 1 = 17
	return (VWGQuadratic(x, par) * v2BackgroundPol3(x, &par[12]) + par[16] * CrystalBallExtended(x, &par[5])) / (VWGQuadratic(x, par) + CrystalBallExtended(x, &par[5]));
}

Double_t v2Pol3vsMassNA60VWG(Double_t *x, Double_t *par)
{ // bkg 4 + sig 11 + flow bkg 4 + flow sig 1 = 20
	return (VWG(x, par) * v2BackgroundPol3(x, &par[15]) + par[19] * NA60(x, &par[4])) / (VWG(x, par) + NA60(x, &par[4]));
}

////////////////////////////////

Double_t v2PolExpvsMassCB2VWG(Double_t *x, Double_t *par)
{ // bkg 4 + sig 7 + flow bkg 5 + flow sig 1= 17

	return (VWG(x, par) * v2BackgroundPolExp(x, &par[11]) + par[16] * CrystalBallExtended(x, &par[4])) / (VWG(x, par) + CrystalBallExtended(x, &par[4]));
}

Double_t v2PolExpvsMassNA60QVWG(Double_t *x, Double_t *par)
{ // bkg 5 + sig 11 + flow bkg 5 + flow sig 1 = 22
	return (VWGQuadratic(x, par) * v2BackgroundPolExp(x, &par[16]) + par[21] * NA60(x, &par[5])) / (VWGQuadratic(x, par) + NA60(x, &par[5]));
}

Double_t v2PolExpvsMassNA60Pol2Pol3(Double_t *x, Double_t *par)
{ // bkg 6 + sig 11 + flow bkg 5 + flow sig 1 = 23
	return (Pol2Pol3(x, par) * v2BackgroundPolExp(x, &par[17]) + par[22] * NA60(x, &par[6])) / (Pol2Pol3(x, par) + NA60(x, &par[6]));
}

Double_t v2PolExpvsMassCB2Pol2Pol3(Double_t *x, Double_t *par)
{ // bkg 6 + sig 7 + flow bkg 5 + flow sig 1 = 19
	return (Pol2Pol3(x, par) * v2BackgroundPolExp(x, &par[13]) + par[18] * CrystalBallExtended(x, &par[6])) / (Pol2Pol3(x, par) + CrystalBallExtended(x, &par[6]));
}

Double_t v2PolExpvsMassCB2QVWG(Double_t *x, Double_t *par)
{ // bkg 5 + sig 7 + flow bkg 5 + flow sig 1 = 18
	return (VWGQuadratic(x, par) * v2BackgroundPolExp(x, &par[12]) + par[17] * CrystalBallExtended(x, &par[5])) / (VWGQuadratic(x, par) + CrystalBallExtended(x, &par[5]));
}

Double_t v2PolExpvsMassNA60VWG(Double_t *x, Double_t *par)
{ // bkg 4 + sig 11 + flow bkg 5 + flow sig 1 = 21
	return (VWG(x, par) * v2BackgroundPolExp(x, &par[15]) + par[20] * NA60(x, &par[4])) / (VWG(x, par) + NA60(x, &par[4]));
}

////////////////

//___________________________________________________________________________________________________________
Int_t GetNPar(Efunction fName)
{
	switch (fName)
	{
	case kVWG:
	case kDoubleExp:
	case kv2BackgroundPol3:
		return 4;
		break;

	case kExp:
		return 2;
		break;

	case kVWGQuadratic:
	case kCB:
	case kv2BackgroundPolExp:
		return 5;
		break;

	case kPol2Pol3:
		return 6;
		break;

	case kCBExtended:
		return 7;
		break;

	case kNA60:
		return 11;
		break;

	case kv2Background:
		return 3;
		break;

	// Group 1
	case kv2Pol2vsMassCB2VWG:
		return 15; // Number of parameters for v2Pol2
		break;
	case kv2Pol3vsMassCB2VWG:
		return 16; // One more parameter than v2Pol2
		break;

	// Group 2
	case kv2Pol2vsMassNA60VWG:
		return 19; // Number of parameters for v2Pol2
		break;
	case kv2Pol3vsMassNA60VWG:
		return 20; // One more parameter than v2Pol2
		break;

	// Group 3
	case kv2Pol2vsMassNA60QVWG:
		return 20; // Number of parameters for v2Pol2
		break;
	case kv2Pol3vsMassNA60QVWG:
		return 21; // One more parameter than v2Pol2
		break;

	// Group 4
	case kv2Pol2vsMassNA60Pol2Pol3:
	case kv2PolExpvsMassNA60VWG:
		return 21; // Number of parameters for v2Pol2
		break;
	case kv2Pol3vsMassNA60Pol2Pol3:
	case kv2PolExpvsMassNA60QVWG:
		return 22; // One more parameter than v2Pol2
		break;

	case kv2PolExpvsMassNA60Pol2Pol3:
		return 23; // One more parameter than v2Pol2
		break;

	// Group 5
	case kv2Pol2vsMassCB2Pol2Pol3:
	case kv2PolExpvsMassCB2VWG:
		return 17; // Number of parameters for v2Pol2
		break;

	case kv2Pol3vsMassCB2Pol2Pol3:
	case kv2PolExpvsMassCB2QVWG:
		return 18; // One more parameter than v2Pol2
		break;

	case kv2PolExpvsMassCB2Pol2Pol3:
		return 19; // One more parameter than v2Pol2
		break;

	// Group 6
	case kv2Pol2vsMassCB2QVWG:
		return 16; // Number of parameters for v2Pol2
		break;
	case kv2Pol3vsMassCB2QVWG:
		return 17; // One more parameter than v2Pol2
		break;
	}
}
