/*LOGistic Regression Macro. Additive model.*/

%macro LOG_add_main(xvar=, m=1 ,drop=T);

ods output "Parameter Estimates"=LOG_Est(keep= Variable df Estimate StdErr ProbChiSq);
ods output "Fit Statistics"=LOG_Fit(keep= Criterion InterceptAndCovariates);
ods output "Observations Summary" = LOG_Num (keep= Label N);

proc LOGistic data=one_sub (keep = ESRD SNP&i &xvar) descending &logoptions;
model ESRD=SNP&i &xvar;
run; 

%let syserr1 = &syserr;
data _AUX; 
LENGTH SYSERR1 3 ;
LENGTH SNP_NO 8 ;
LENGTH _SNP $16;
SNP_NO=symget("i");
SYSERR1=symget("syserr1");
_SNP=compress("SNP", snp_no);
run;

%*if &i=1 %then %drop_data (data = LOG_Est);
%if %sysfunc(exist(LOG_Est)) and %sysfunc(exist(LOG_Fit)) and %sysfunc(exist(LOG_Num)) %then %do;
%log_add_main_append
%end;
%mend LOG_add_main;

%macro log_add_main_append;

data LOG_diag1;
 set LOG_est end=_last_;
 if Estimate = "." then LOG_ADD&m._DIAG1=1;
 else  LOG_ADD&m._DIAG1=0;

 if StdErr="." then LOG_ADD&m._DIAG2=1;
else  LOG_ADD&m._DIAG2=0;

LOG_ADD&m._DIAG_BETA+LOG_ADD&m._DIAG1;
LOG_ADD&m._DIAG_SE+LOG_ADD&m._DIAG2;
if _last_;
keep LOG_ADD&m._DIAG_BETA LOG_ADD&m._DIAG_SE;
run;

/*Only keep the results from SNP&i*/
data LOG_est;
set LOG_est;
if Variable="SNP&i";
run;


%append_snp(LOG_add&m._est, LOG_est, &drop);
%append_snp(LOG_add&m._num, LOG_num, &drop);
%append_snp(LOG_add&m._fit, LOG_fit, &drop);
%append_snp(LOG_add&m._AUX, _AUX,    &drop);
%append_snp(LOG_add&m._diag1, LOG_diag1,    &drop);
%mend log_add_main_append;

%macro LOG_add_init(xvar=, m=1, drop = T);
data LOG_add&m._est;
length Variable $ 9;
length Estimate 8;
length StdErr 8;
length ProbChiSq 8;
length df 8;
length _SNP $16;
length SNP_NO 8;
stop;
run;

data LOG_add&m._fit;
length Criterion $9;
length InterceptAndCovariates 8;
length _SNP $16;
length SNP_NO 8;
stop;
run;

data LOG_add&m._Num;
length Label $27;
length N 8;
length _SNP $16;
length SNP_NO 8;
stop;
run;

data LOG_add&m._AUX;
length SYSERR1 3;
length _SNP $16;
length SNP_NO 8;
stop;
run;

data LOG_add&m._DIAG1;
length LOG_ADD&m._DIAG_BETA 8;
length LOG_ADD&m._DIAG_SE 8;
length _SNP $16;
length SNP_NO 8;
stop;
run;
%mend LOG_add_init;

%macro LOG_add_finish(xvar=, m=1, drop = T);

data LOG_add&m._fit1;
set LOG_add&m._fit;
retain AIC;
if  Criterion="AIC" then AIC=InterceptAndCovariates;
if  Criterion="-2 Log L" then Min_2LOG=InterceptAndCovariates;
if  Criterion="-2 Log L" ;
drop Criterion InterceptAndCovariates;
run;

data LOG_add&m._Num1;
set LOG_add&m._Num; 
where label="Number of Observations Used";
keep N SNP_no;
run;

data LOG_add&m;
merge LOG_add&m._est LOG_add&m._fit1 LOG_add&m._Num1 LOG_add&m._aux LOG_add&m._diag1;
by SNP_no;
OR=exp(estimate);
run;

data LOG_add&m (rename= ( 
                  Variable  = LOG_ADD&m._VAR 
                  df        = LOG_ADD&m._DF 
                  Estimate  = LOG_ADD&m._BETA 
                  StdErr    = LOG_ADD&m._SE 
                  ProbChiSq = LOG_ADD&m._P 
                  AIC       = LOG_ADD&m._AIC 
                  Min_2LOG  = LOG_ADD&m._MIN2LOG 
                  N         = LOG_ADD&m._N 
                  OR        = LOG_ADD&m._OR 
                  SYSERR1   = LOG_ADD&m._SYSERR)) ;
set LOG_add&m;

keep Variable Estimate StdErr ProbChiSq df SNP_NO AIC Min_2LOG N OR SYSERR1 LOG_ADD&m._DIAG_BETA LOG_ADD&m._DIAG_SE;

LABEL  Variable                ="Variable in the Logistic Regression Model. SNP number is the modified SNP order in the final analysis geno file"
       df                      ="Degree of freedom of each SNP in the Logistic Regression Model"
       Estimate                ="Beta Estimate of each SNP in the Logistic Regression Model"
       StdErr                  ="Standard Error of each SNP in the Logistic Regression Model"
       ProbChiSq               ="P-Value of each SNP in the Logistic Regression Model"
       SNP_NO                  ="Modified SNP order in the final analysis work._geno file"
       AIC                     ="AIC (Akaike information criterion) from fit statistics in the full Logistic Regression Model"
       Min_2log                ="-2 log Likelihood in the full Logistic Regression Model"
       N                       ="Number of observations used in the full Logistic Regression Model"
       OR                      ="Odds ratio of each SNP in the Logistic Regression Model"
       SYSERR1                 ="System error message in the full Logistic Regression Model"
       LOG_ADD&m._DIAG_BETA    ="Count of Beta Estimates that are missing"
       LOG_ADD&m._DIAG_SE      ="Count of Standard Errors for Beta Estimates that are missing";
run;

proc printto print = myLOG; run;
proc means data = LOG_add&m (drop=snp_no);
Title "Summary of the results from LOG_add&m for &nsnps SNPs";
Title2 "Logistic reg model: ESRD = &xvar +SNP";
Title3 "logoptions:= &logoptions" ;
proc printto; run;
%let merges=&merges LOG_add&m;
%mend LOG_add_finish;