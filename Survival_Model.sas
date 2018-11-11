/* Survival Model Macro*/

%macro srv_add_jm;
/* Auxiliary macro creating data for JM */ 

/* data srv_est created by PROC LIFEREG */

data srv&m._est1;
 set srv_est;
 length Parameter $ 20;
 nx = _n_ - 1;
 Parameter = cats("bs", nx);
 sdist=symget("sdist");
 if sdist = "llogistic" and SurvEffect     = "Scale"       then Parameter = 'sgma';
 if sdist = "lnormal" and SurvEffect       = "Scale"       then Parameter = 'sigma';
 if sdist = "weibull" and SurvEffect       = "Weibull Shape"       then Parameter = 'gamma';
 drop nx sdist;

/* SKIP THIS PARA
if SurvEffect     = "Intercept"      then Parameter = 'bs0' ;
if SurvEffect     = "SNP&i"          then Parameter = 'bs1' ;
if SurvEffect     = "u0"             then Parameter = 'bs2'  ;
if SurvEffect     = "u1"             then Parameter = 'bs3'  ;
if SurvEffect     = "BASE_GFR"       then Parameter = 'bs4' ;
if SurvEffect     = "Scale"          then delete;
if SurvEffect     = "Weibull Shape"  then delete;
*/
run;

%mend srv_add_jm;

%macro srv_add_main(xvar = BASE_GFR, m = 1, data = one_sub, sdist=weibull, drop = T);

/*---------------------------------Proc LifeReg-Additive Model  ------------------------------------*/ 
*ods select none;
ods output "Analysis of Parameter Estimates"=srv_Est(drop=ChiSq rename=(Parameter=SurvEffect));
ods output "Model Fit Statistics"=srv_Fit;
ods output "Observations Summary"=srv_Num (keep= Label N);

proc lifereg data = work.&data(keep = ESRDTIME ESRD &xvar SNP&i ) &srvoptions;
model ESRDTIME*ESRD(0) = SNP&i &xvar / dist = &sdist ;
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

%*if &i=1 %then %drop_data (data = srv_Est);
%if %sysfunc(exist(srv_Est)) and %sysfunc(exist(srv_Fit)) and %sysfunc(exist(srv_Num)) %then %do;
%srv_add_main_append;
%end;
%mend srv_add_main;

%macro srv_add_main_append;

data srv_scale ;
set srv_est ;
if SurvEffect = "Scale" then SRV_ADD&m._SCALE = Estimate;
if SurvEffect = "Scale" then SRV_ADD&m._SCALE_SE = StdErr;
if SurvEffect = "Scale";
keep SRV_ADD&m._SCALE SRV_ADD&m._SCALE_SE;
run;

%if (%upcase(&drop) ne T) %then %srv_add_jm;

data srv_diag1;
 set srv_est end=_last_;
 if Estimate = "." then srv_ADD&m._DIAG1=1;
 else  srv_ADD&m._DIAG1=0;

 if StdErr="." then srv_ADD&m._DIAG2=1;
else  srv_ADD&m._DIAG2=0;

SRV_ADD&m._DIAG_BETA+srv_ADD&m._DIAG1;
SRV_ADD&m._DIAG_SE+srv_ADD&m._DIAG2;
if _last_;
keep SRV_ADD&m._DIAG_BETA SRV_ADD&m._DIAG_SE;
run;

/*Only keep the results from SNP&i*/
data srv_est2;
set srv_est;
if SurvEffect in ("Scale","Weibull Shape")then delete;
if SurvEffect="SNP&i";
run;



%append_snp(srv_add&m._est, srv_est2, &drop);
%append_snp(srv_add&m._scale, srv_scale, &drop);
%append_snp(srv_add&m._fit, srv_fit, &drop);
%append_snp(srv_add&m._num, srv_num, &drop);
%append_snp(srv_add&m._AUX, _AUX, &drop);
%append_snp(srv_add&m._diag1, srv_diag1, &drop);

%mend srv_add_main_append;

%macro srv_add_init(xvar = BASE_GFR, m = 1, data = one_sub, sdist=weilbull, drop = T);
data srv_add&m._est;
length SurvEffect $ 13;
length Estimate 8;
length StdErr 8;
length LowerCL 8;
length UpperCL 8;
length ProbChiSq 8;
length df 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data srv_add&m._scale;
length SRV_ADD&m._SCALE 8;
length SRV_ADD&m._SCALE_SE 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data srv_add&m._fit;
length Criterion $35;
length Value 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data srv_add&m._num;
length Label $27;
length N 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data srv_add&m._AUX;
length SYSERR1 3;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data srv_add&m._diag1;
length SRV_ADD&m._DIAG_BETA 8;
length SRV_ADD&m._DIAG_SE 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;
%mend srv_add_init;

%macro srv_add_finish(xvar = BASE_GFR, m = 1, data = one_sub, sdist=weilbull, drop = T);
data srv_add&m._fit1;
set srv_add&m._fit;
retain Min_2LOG;
if Criterion="-2 Log Likelihood" then Min_2LOG=Value;
if Criterion="AIC (smaller is better)" then AIC=Value;
if Criterion="AIC (smaller is better)";
drop Criterion Value;
run;

data srv_add&m._num1;
set srv_add&m._num; 
where label="Number of Observations Used";
drop label;
run;

data srv_add&m;
merge srv_add&m._est srv_add&m._scale srv_add&m._fit1 srv_add&m._num1 srv_add&m._AUX srv_add&m._DIAG1;
by SNP_no;
run;

data srv_add&m (rename= ( 
                SurvEffect=SRV_ADD&m._VAR
                df=SRV_ADD&m._DF
                Estimate=SRV_ADD&m._BETA
                StdErr=SRV_ADD&m._SE
                ProbChiSq=SRV_ADD&m._P
                AIC=SRV_ADD&m._AIC
                Min_2LOG=SRV_ADD&m._MIN2LOG
                N=SRV_ADD&m._N 
                SYSERR1=SRV_ADD&m._SYSERR)) ;
set SRV_add&m;
keep SurvEffect Estimate StdErr df ProbChiSq SRV_ADD&m._SCALE SRV_ADD&m._SCALE SRV_ADD&m._SCALE_SE SNP_NO AIC Min_2LOG N SYSERR1 SRV_ADD&m._DIAG_BETA SRV_ADD&m._DIAG_se;

LABEL  SurvEffect               ="Variable in the Survival Model. SNP number is the modified SNP order in the final analysis geno file"
       df                       ="Degree of freedom of each SNP in the Survival Model"
       Estimate                 ="Beta Estimate of each SNP in the Survival Model"
       StdErr                   ="Standard Error of each SNP in the Survival Model"
       ProbChiSq                ="P-Value of each SNP in the Survival Model"
       SRV_ADD&m._SCALE         ="Scale parameter in the Survival Model"
       SRV_ADD&m._SCALE_SE      ="Standard Error of Scale parameter in the Survival Model"
       SNP_NO                   ="Modified SNP order in the final analysis work._geno file"
       AIC                      ="AIC (Akaike information criterion) from fit statistics in the full Survival Model"
       Min_2LOG                 ="-2 LOG Likelihood in the full Survival Model"
       N                        ="Number of observations used in the full Survival Model"
       SYSERR1                  ="System error message in the full Survival Model"
       SRV_ADD&m._DIAG_BETA     ="Count of Beta Estimates that are missing. Scale parameter is also included."
       SRV_ADD&m._DIAG_SE       ="Count of Standard Errors for Beta Estimates that are missing. Standard Error of scale parameter is also included";
run;

proc printto print = myLOG; run;
proc means data = SRV_add&m (drop=snp_no);
Title "Summary of the results from SRV_add&m for &nsnps SNPs";
Title2 "Lifereg model: esrdtime*esrd(0)= SNP &xvar / dist = &sdist";
Title3 "srvoptions:= &srvoptions" ;
proc printto; run;
       
       
%let merges=&merges SRV_add&m;
%mend SRV_add_finish;