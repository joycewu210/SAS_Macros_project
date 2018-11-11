/* Linear Mixed Model Macro. Additive model*/

%macro lmm_add_jm;
/* Auxiliary macro creating data for JM */ 

/* data Randm created by PROC MIXED */

/*** Extract random intercept and random slope from LMM and merge it with work.one_sub for survival model in joint model**/
data Randm;                  /* Random effects: one record per subject */
set Randm;
 retain u0 ;
 if effect = "Intercept" then u0   = estimate;
 if effect = "FUTIME"    then u1 = estimate;
 if effect = "FUTIME";
 KEEP _numid u0 u1;
run;

data WORK.one_sub_rand&m;    /* Data set for survival model in JM */
merge WORK.one_sub Randm;
by _numid;
run;

%drop_data (data = Randm);

/***Extract estimates from LMM for JM analysis**/

data work.lmm&m._est1;
 set work.lmm_est;
 length Parameter $ 20;
 nx = _n_ - 1;
 Parameter = cats("bl", nx);
 drop nx;

/*  SKIP THIS PARA
if Effect     = "Intercept"     then Parameter = 'bl0' ;
if Effect     = "SNP&i"         then Parameter = 'bl1' ;
if Effect     = "FUTIME"        then Parameter = 'bl2' ;
if Effect     = "SNP&i*FUTIME"  then Parameter = 'bl3' ;
if Effect     = "BASE_GFR"      then Parameter = 'bl4' ;
*/

run;

data work.lmm&m._cpe1;
 set work.lmm_cpe;
 length Parameter $ 5;
 if CovParm    = "FA(1,1)"       then Parameter = 'a11' ;
 if CovParm    = "FA(2,1)"       then Parameter = 'a12' ;
 if CovParm    = "FA(2,2)"       then Parameter = 'a22' ;
 if CovParm    = "Residual"      then Parameter = 'S2'  ; 
run;

/* 
Note: Datasets created
  work.one_sub_rand&m 
  work.lmm&m._est1
  work.lmm&m._cpe1
*/
%mend lmm_add_jm;


%macro lmm_add_main(xvar = , m = 1, drop = T);

/* ---- PROC MIXED_ADDITIVE MODEL ---- */
*ods select none;
ods output "Solution for Fixed Effects"     = lmm_Est(drop = tValue);
ods output "Covariance Parameter Estimates" = lmm_Cpe(drop = Subject ZValue ProbZ);
ods output "Fit Statistics"                 = lmm_Fit;
ods output "Number of Observations"         = lmm_Num(drop = NObsRead NObsUsed SumFreqsRead SumFreqsUsed);
ods output "Solution for Random Effects"    = Randm;

proc mixed data = work.mul_sub(keep = GFR &xvar FUTIME _numid SNP&i) covtest method = ml &lmmoptions;
 class _numid;
 model gfr = SNP&i FUTIME SNP&i*FUTIME &xvar/ solution ;
 random intercept FUTIME/ subject = _numid type = fa0(2) solution;
run; 



%let syserr1 = &syserr;
data _AUX; 
LENGTH SYSERR1 3 ;
LENGTH SNP_NO 8 ;
LENGTH _SNP $16;
SNP_NO=symget("i");
SYSERR1=symget("syserr1");
_SNP=compress("SNP",snp_no);
run;


%*if &i=2 %then %drop_data (data = lmm_Est);
%if %sysfunc(exist(lmm_Est)) and %sysfunc(exist(lmm_Cpe)) and %sysfunc(exist(lmm_Fit)) and %sysfunc(exist(lmm_Num)) and %sysfunc(exist(Randm))%then %do;
%lmm_add_main_append;
%end;
%mend lmm_add_main;

%macro lmm_add_main_append;

%if (%upcase(&drop) ne T) %then %lmm_add_jm;

data lmm_diag1;
 set lmm_est end=_last_;
 if Estimate = "." then LMM_ADD&m._DIAG1=1;
 else  LMM_ADD&m._DIAG1=0;

 if StdErr="." then LMM_ADD&m._DIAG2=1;
else  LMM_ADD&m._DIAG2=0;

LMM_ADD&m._DIAG_BETA+LMM_ADD&m._DIAG1;
LMM_ADD&m._DIAG_SE+LMM_ADD&m._DIAG2;
if _last_;
KEEP LMM_ADD&m._DIAG_BETA LMM_ADD&m._DIAG_SE;
run;

data lmm_diag2;
 set lmm_cpe end=_last_;
 if Estimate = "." then LMM_ADD&m._DIAG3=1;
 else  LMM_ADD&m._DIAG3=0;

 if StdErr="." then LMM_ADD&m._DIAG4=1;
else  LMM_ADD&m._DIAG4=0;
LMM_ADD&m._DIAG_CP+LMM_ADD&m._DIAG3;
LMM_ADD&m._DIAG_CP_SE+LMM_ADD&m._DIAG4;
if _last_;
KEEP LMM_ADD&m._DIAG_CP LMM_ADD&m._DIAG_CP_SE;
run;

/* Only keep the results from interaction term for final results */
/*
data lmm_est;
length Effect $15;
 set lmm_est (rename=(Effect=Effect1));
 Effect=Effect1;
 if Effect = "SNP&i*FUTIME";
 drop Effect1;
run;

*/

data lmm_est;
set lmm_est;
retain 

LMM_ADD&m._INT_BETA
LMM_ADD&m._INT_SE
LMM_ADD&m._SNP_BETA
LMM_ADD&m._SNP_SE
LMM_ADD&m._SNP_P
LMM_ADD&m._TIME_BETA
LMM_ADD&m._TIME_SE
LMM_ADD&m._TIME_P
LMM_ADD&m._SNP_TIME_BETA
LMM_ADD&m._SNP_TIME_SE
LMM_ADD&m._SNP_TIME_P;

if Effect = "Intercept" then do;
LMM_ADD&m._INT_BETA=Estimate;
LMM_ADD&m._INT_SE=StdErr;
end;

if Effect = "SNP&i" then do;
LMM_ADD&m._SNP_BETA=Estimate;
LMM_ADD&m._SNP_SE=StdErr;
LMM_ADD&m._SNP_P=Probt;
end;

if Effect = "FUTIME" then do;
LMM_ADD&m._TIME_BETA=Estimate;
LMM_ADD&m._TIME_SE=StdErr;
LMM_ADD&m._TIME_P=Probt;
end;

if Effect = "SNP&i*FUTIME" then do;
LMM_ADD&m._SNP_TIME_BETA=Estimate;
LMM_ADD&m._SNP_TIME_SE=StdErr;
LMM_ADD&m._SNP_TIME_P=Probt;
LMM_ADD&m._SNP_TIME_DF=DF;
end;

if Effect = "SNP&i*FUTIME";
drop Effect Estimate StdErr Probt DF;
RUN;





*proc print data=lmm_est;run;
*proc contents data=lmm_est;run;

%append_snp(lmm_add&m._est, lmm_est, &drop);
%append_snp(lmm_add&m._cpe, lmm_cpe, &drop);
%append_snp(lmm_add&m._fit, lmm_fit, &drop);
%append_snp(lmm_add&m._num, lmm_num, &drop);
%append_snp(lmm_add&m._AUX, _AUX,    &drop);
%append_snp(lmm_add&m._diag1, lmm_diag1,    &drop);
%append_snp(lmm_add&m._diag2, lmm_diag2,    &drop);

%mend lmm_add_main_append;

%macro lmm_add_init(xvar = , m = 1, drop = T);
data lmm_add&m._est;

length LMM_ADD&m._INT_BETA 8;
length LMM_ADD&m._INT_SE 8;
length LMM_ADD&m._SNP_BETA 8;
length LMM_ADD&m._SNP_SE 8;
length LMM_ADD&m._SNP_P 8;
length LMM_ADD&m._TIME_BETA 8;
length LMM_ADD&m._TIME_SE 8;
length LMM_ADD&m._TIME_P 8;
length LMM_ADD&m._SNP_TIME_BETA 8;
length LMM_ADD&m._SNP_TIME_SE 8;
length LMM_ADD&m._SNP_TIME_P 8;
length LMM_ADD&m._SNP_TIME_DF 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data lmm_add&m._fit;
length Descr $25;
length Value 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data lmm_add&m._Num;
length Label $31;
length N 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data lmm_add&m._cpe;
length CovParm $9;
length Estimate 8;
length StdErr 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data lmm_add&m._AUX;
length SYSERR1 3;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data lmm_add&m._diag1;
length LMM_ADD&m._DIAG_BETA 8;
length LMM_ADD&m._DIAG_SE 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data lmm_add&m._diag2;
length LMM_ADD&m._DIAG_CP 8;
length LMM_ADD&m._DIAG_CP_SE 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

%mend lmm_add_init;

%macro lmm_add_finish(xvar = , m = 1, drop = T);
data lmm_add&m._fit1;
set lmm_add&m._fit;
retain min_2log;
if  Descr = "AIC (smaller is better)" then AIC = Value;
if  Descr = "-2 Log Likelihood" then min_2log = Value;
if  Descr = "AIC (smaller is better)" ;
drop Descr Value;
run;

data lmm_add&m._Num1;
set lmm_add&m._Num; 
where label="Number of Observations Used";
drop label;
run;

data lmm_add&m._cpe1;
set lmm_add&m._cpe;
retain a11 a12 a22 Residual;
if CovParm="FA(1,1)" then a11=Estimate;
if CovParm="FA(2,1)" then a12=Estimate;
if CovParm="FA(2,2)" then a22=Estimate;
if CovParm="Residual" then Residual=Estimate;
drop CovParm Estimate StdErr;
if CovParm="Residual";
run;

proc print data=lmm_add&m._cpe;run;

data lmm_add&m;
merge lmm_add&m._est lmm_add&m._fit1 lmm_add&m._cpe1 lmm_add&m._Num1 lmm_add&m._AUX lmm_add&m._DIAG1 lmm_add&m._DIAG2;
by snp_no;
run;

data lmm_add&m (rename= ( 

                 AIC       = LMM_ADD&m._AIC 
                 Min_2log  = LMM_ADD&m._MIN2LOG 
                 N         = LMM_ADD&m._N 
                 SYSERR1   = LMM_ADD&m._SYSERR)) ;
set lmm_add&m;
keep 
LMM_ADD&m._INT_BETA
LMM_ADD&m._INT_SE
LMM_ADD&m._SNP_BETA
LMM_ADD&m._SNP_SE
LMM_ADD&m._SNP_P
LMM_ADD&m._TIME_BETA
LMM_ADD&m._TIME_SE
LMM_ADD&m._TIME_P
LMM_ADD&m._SNP_TIME_BETA
LMM_ADD&m._SNP_TIME_SE
LMM_ADD&m._SNP_TIME_P
LMM_ADD&m._SNP_TIME_DF

SNP_NO AIC Min_2log N SYSERR1 LMM_ADD&m._DIAG_BETA LMM_ADD&m._DIAG_SE LMM_ADD&m._DIAG_CP LMM_ADD&m._DIAG_CP_SE;     

LABEL  
LMM_ADD&m._INT_BETA       ="Beta Estimate of intercept in the Linear Mixed Model"
LMM_ADD&m._INT_SE         ="Standard Error of intercept in the Linear Mixed Model" 
LMM_ADD&m._SNP_BETA       ="Beta Estimate of SNP in the Linear Mixed Model"
LMM_ADD&m._SNP_SE         ="Standard Error of SNP in the Linear Mixed Model" 
LMM_ADD&m._SNP_P          ="P-Value of SNP in the Linear Mixed Model" 
LMM_ADD&m._TIME_BETA      ="Beta Estimate of FUTIME in the Linear Mixed Model"
LMM_ADD&m._TIME_SE        ="Standard Error of FUTIME in the Linear Mixed Model"
LMM_ADD&m._TIME_P         ="P-Value of FUTIME in the Linear Mixed Model"
LMM_ADD&m._SNP_TIME_BETA  ="Beta Estimate of the interaction term (SNP*FUTIME) in the Linear Mixed Model"
LMM_ADD&m._SNP_TIME_SE    ="Standard Error of the interaction term (SNP*FUTIME) in the Linear Mixed Model"
LMM_ADD&m._SNP_TIME_P     ="P-Value of the interaction term (SNP*FUTIME) in the Linear Mixed Model"
LMM_ADD&m._SNP_TIME_DF    ="Degree of freedom of the interaction term (SNP*FUTIME) in the Linear Mixed Model"


SNP_NO                    ="Modified SNP order in the final analysis work._geno file"
AIC                       ="AIC (Akaike information criterion) from fit statistics in the full Linear Mixed Model"
Min_2log                  ="-2 Log Likelihood in the full Linear Mixed Model"
N                         ="Number of observations used in the full Linear Mixed Model"
SYSERR1                   ="System error message in the full Linear Mixed Model"
LMM_ADD&m._DIAG_BETA      ="Count of Beta Estimates that are missing"
LMM_ADD&m._DIAG_SE        ="Count of Standard Errors for Beta Estimates that are missing"
LMM_ADD&m._DIAG_CP        ="Count of Covariance Parameter Estimates that are missing"
LMM_ADD&m._DIAG_CP_SE     ="Count of Standard Errors for Covariance Parameter Estimates that are missing";
run;

proc printto print = mylog; run;
proc means data = lmm_add&m (drop=snp_no);
Title "Summary of the results from lmm_add&m for &nsnps SNPs";
Title2 "LMM Model: gfr = SNP FUTIME SNP*FUTIME &xvar" ;
Title3 "lmmoptions:= &lmmoptions" ;
proc printto; run;

%let merges=&merges lmm_add&m;
%mend lmm_add_finish;

