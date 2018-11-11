/* Cox Model Macro*/  

%macro cox_add_main(xvar = BASE_GFR, m = 1, data = one_sub, drop = T);

/*---------------------------------Proc PHReg-Additive Model  ------------------------------------*/ 
*ods select none;
ods output "Maximum Likelihood Estimates of Model Parameters"=cox_Est(drop=ChiSq HazardRatio rename=(Parameter=CoxEffect));
ods output "Model Fit Statistics"=cox_Fit(drop=WithoutCovariates);
ods output "Number of Observations"=cox_Num (drop= NObsRead SumFreqsRead SumFreqsUsed rename=(NObsUsed=N));

proc phreg data = work.&data(keep = ESRDTIME ESRD &xvar SNP&i );
model ESRDTIME*ESRD(0) = SNP&i &xvar;
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

%*if &i=1 %then %drop_data (data = cox_Est);
%if %sysfunc(exist(cox_Est))and %sysfunc(exist(cox_Fit)) and %sysfunc(exist(cox_Num)) %then %do;
%cox_add_main_append;
%end;
%mend cox_add_main;

%macro cox_add_main_append;

data cox_diag1;
set cox_est end=_last_;
if Estimate = "." then cox_ADD&m._DIAG1=1;
else  cox_ADD&m._DIAG1=0;

if StdErr="." then cox_ADD&m._DIAG2=1;
else  cox_ADD&m._DIAG2=0;

COX_ADD&m._DIAG_BETA+cox_ADD&m._DIAG1;
COX_ADD&m._DIAG_SE+cox_ADD&m._DIAG2;
if _last_;
keep COX_ADD&m._DIAG_BETA COX_ADD&m._DIAG_SE;
run;

/*Only keep the results from SNP&i*/
data cox_est2;
set cox_est;
if CoxEffect="SNP&i";
run;

%append_snp(cox_add&m._est, cox_est2, &drop);
%append_snp(cox_add&m._fit, cox_fit, &drop);
%append_snp(cox_add&m._num, cox_num, &drop);
%append_snp(cox_add&m._AUX, _AUX, &drop);
%append_snp(cox_add&m._diag1, cox_diag1, &drop);

%mend cox_add_main_append;

%macro cox_add_init(xvar = BASE_GFR, m = 1, data = one_sub, drop = T);
data cox_add&m._est;
length CoxEffect $ 13;
length Estimate 8;
length StdErr 8;
length ProbChiSq 8;
length df 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data cox_add&m._fit;
length Criterion $35;
length WithCovariates 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;


data cox_add&m._num;
length N 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data cox_add&m._AUX;
length SYSERR1 3;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data cox_add&m._diag1;
length COX_ADD&m._DIAG_BETA 8;
length COX_ADD&m._DIAG_SE 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;
%mend cox_add_init;

%macro cox_add_finish(xvar = BASE_GFR, m = 1, data = one_sub, drop = T);

data cox_add&m._fit1;
set cox_add&m._fit;
retain Min_2LOG;
if Criterion="-2 LOG L" then Min_2LOG=WithCovariates;
if Criterion="AIC" then AIC=WithCovariates;
if Criterion="AIC";
drop Criterion WithCovariates;
run;

data cox_add&m;
merge cox_add&m._est cox_add&m._fit1 cox_add&m._num cox_add&m._AUX cox_add&m._DIAG1;
by SNP_no;
run;

data cox_add&m (rename= ( 
                CoxEffect=COX_ADD&m._VAR
                df=COX_ADD&m._DF
                Estimate=COX_ADD&m._BETA
				AIC=COX_ADD&m._AIC
                Min_2LOG=COX_ADD&m._MIN2LOG
                StdErr=COX_ADD&m._SE
                ProbChiSq=COX_ADD&m._P
                N=COX_ADD&m._N 
                SYSERR1=COX_ADD&m._SYSERR)) ;
set COX_add&m;
keep CoxEffect Estimate StdErr df ProbChiSq  SNP_NO  N SYSERR1 AIC Min_2LOG COX_ADD&m._DIAG_BETA COX_ADD&m._DIAG_se;

LABEL  CoxEffect                ="Variable in the Cox Model. SNP number is the modified SNP order in the final analysis geno file"
       df                       ="Degree of freedom of each SNP in the Cox Model"
       Estimate                 ="Beta Estimate of each SNP in the Cox Model"
       StdErr                   ="Standard Error of each SNP in the Cox Model"
       ProbChiSq                ="P-Value of each SNP in the Cox Model"
	    AIC                     ="AIC (Akaike information criterion) from fit statistics in the full Cox Model"
       Min_2LOG                 ="-2 LOG Likelihood in the full Cox Model"
       SNP_NO                   ="Modified SNP order in the final analysis work._geno file"
       N                        ="Number of observations used in the full Cox Model"
       SYSERR1                  ="System error message in the full Cox Model"
       COX_ADD&m._DIAG_BETA     ="Count of Beta Estimates that are missing."
       COX_ADD&m._DIAG_SE       ="Count of Standard Errors for Beta Estimates that are missing.";
run;

proc printto print = myLOG; run;
proc means data = COX_add&m (drop=snp_no);
Title "Summary of the results from COX_add&m for &nsnps SNPs";
Title2 "PHreg model: esrdtime*esrd(0)= SNP &xvar";
proc printto; run;
       
       
%let merges=&merges COX_add&m;
%mend COX_add_finish;
