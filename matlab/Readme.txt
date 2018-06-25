17 April 2018

Alexandros Kostakis, Tassos Magdalinos and Michalis P. Stamatogiannis, 
"Robust econometric inference for stock return predictability", 
The Review of Financial Studies, Vol. 28, No. 5, 2014, pp. 1506-1553.

This file details the data and computer programs used to produce the 
empirical results in the published paper.

The potential predictors dataset updated up to December 2012 is sourced 
from Amit Goyal's website: http://www.hec.unil.ch/agoyal/

We use S&P 500 value-weighted log excess returns to proxy for 
excess market returns.

*****
Data:
*****

monthly.xlsx:  monthly data with 13 series with rows corresponding to
dates and columns to variables ordered as:

1	year and month
2	dividend payout ratio (DE)
3	long-term yield (LTY)
4	dividend yield (DY)
5	dividend-price ratio (DP)
6	T-bill rate (TBL) 
7	earnings-price ratio (EP)
8	book-to-market value ratio (BM)
9	inflation rate (INF)
10	default yield spread (DFY)
11	net equity expansion (NTIS)
12	term spread (TMS)
13	S&P 500 value-weighted log excess returns (LOG EXCESS VW)


quarterly.xlsx:  quarterly data, 14 series with rows corresponding to
dates and columns to variables ordered as:

1	year and quarter
2	dividend payout ratio (DE)
3	long-term yield (LTY)
4	dividend yield (DY)
5	dividend-price ratio (DP)
6	T-bill rate (TBL) 
7	earnings-price ratio (EP)
8	book-to-market value ratio (BM)
9	inflation rate (INF)
10	default yield spread (DFY)
11	net equity expansion (NTIS)
12	term spread (TMS)
13	S&P 500 value-weighted log excess returns (LOG EXCESS VW)
14	consumption-wealth ratio (CAY)


****************
Matlab Programs:
****************

ivxlh.m 
This is the Matlab function used to run the long-horizon IVX
estimation and tests of individual and joint significance of 
the predictors. Details about the input required and output
are provided in comments at the beginning of the file. 
Note, that short-horizon regressions can be run with this 
function by setting horizon K=1.


EmpiricalTables.m
This is the Matlab file that calculates the empirical Tables of 
the paper. This file calls both monthly.xlsx and quarterly.xlsx 
data files and the ivxlh.m function to run the regressions. 
All files need to be placed in the same folder.


Please address questions to:  	Michalis Stamatogiannis
				m.stamatogiannis [AT] liverpool.ac.uk
				https://sites.google.com/site/mpstamatogiannis/


