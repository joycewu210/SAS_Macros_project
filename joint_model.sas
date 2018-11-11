/* Joint Model Macro */

%macro jm_weibull; 
      bounds gamma > 0;
      * linpsurv  = bs0 + bs1*SNP&i + bs2*u0 + bs3*u1 + bs4*Base_gfr;    /* Original syntax */
      alpha  = exp(-linpsurv);
      G_t    = exp(-(alpha*esrdtime)**gamma);
      g      = gamma*alpha*((alpha*esrdtime)**(gamma-1))*G_t;
%mend jm_weibull;

%macro jm_llogistic; 
      bounds sgma>0;
      alpha  = exp(-linpsurv/sgma);
      gamma = 1/sgma;
      tmp = 1+ alpha * (esrdtime ** gamma); 
      G_t = 1/tmp;
      g = alpha * gamma* esrdtime ** (gamma-1)/(tmp*tmp);
%mend jm_llogistic;  

%macro jm_lnormal; 
       bounds sigma>0;
       q1 = (log(esrdtime)-linpsurv)/sigma;
       g = exp(-0.5*(q1**2))/((esrdtime*(2*CONSTANT('PI'))**0.5)*sigma);
       G_t = 1- cdf('normal', q1);
%mend jm_lnormal; 

%macro nlmixed_start_init;
%lmm_add_init(xvar = &lmm_xvar, m = &auxm);
%SRV_add_init(xvar = u0 u1 &SRV_xvar, m = &auxm, data = one_sub_rand&auxm, sdist=&jmdist);
%mend  nlmixed_start_init;


%macro nlmixed_start_main;

%lmm_add_main(xvar = &lmm_xvar, m = &auxm, drop = F);
%SRV_add_main(xvar = u0 u1 &SRV_xvar, m = &auxm, data = one_sub_rand&auxm, drop = F, sdist=&jmdist);

%if %sysfunc(exist(lmm&auxm._est1)) and %sysfunc(exist(srv&auxm._est1)) %then %do;
data _null_;
set work.lmm&auxm._est1 end=last;
length linplong $500;
retain linplong;
if _n_=1 then linplong= "(bl0 + u0)+bl1*SNP&i+(bl2 + u1)*FUTIME";
if _n_> 3 then do;
linplong=cats(linplong,"+",Parameter,"*",Effect);
end;
if last then call symput("linplong",linplong);
run;

data srv&auxm._est1;
set srv&auxm._est1 end=last;
length linpsurv $500;
retain linpsurv;
jmdist=symget("jmdist");
if _n_=1 then linpsurv= "bs0";
if _n_> 1 and SurvEffect not in ("Scale","Weibull Shape") then 
linpsurv=cats(linpsurv,"+",Parameter,"*",SurvEffect);
if last then call symput("linpsurv",linpsurv);
if jmdist = "weibull" and  SurvEffect in ("Scale")then delete;
run;

/*----------------Create Initial values for joint model analysis ------------------------------------*/ 


data work.JM_init&m;                                   /* Initial estimates for JM model*/
set lmm&auxm._est1 lmm&auxm._cpe1 srv&auxm._est1;
     keep Parameter Estimate StdErr;
run;

%drop_data (data= one_sub_rand&auxm lmm&auxm._est1 lmm&auxm._cpe1 SRV&auxm._est1 lmm_est lmm_cpe lmm_fit lmm_num lmm_diag1 lmm_diag2 SRV_est SRV_scale SRV_fit SRV_num srv_diag1);
%end;
%mend nlmixed_start_main;

%macro  nlmixed_start_finish;
%lmm_add_finish(xvar = &lmm_xvar, m = &auxm);
%SRV_add_finish(xvar = SNP u0 u1 &SRV_xvar ,  m = &auxm, sdist=&jmdist);
%mend  nlmixed_start_finish;


%macro jm_model;
/* ---- Joint Model - NLMIXED-Weibull Distribution---- */
ods output  ParameterEstimates = JM_Est (drop=tValue Alpha Gradient Lower Upper);
ods output  FitStatistics = JM_Fit;
ods output  Dimensions = JM_Num;
ods output  ConvergenceStatus=JM_Conv (keep=status);
ods output  IterHistory=JM_Ite (drop=Flag1);

Proc nlmixed data=WORK.mul_sub &jmoptions;
 parms / data=work.JM_init&m; 

    /*-----All parameters not assigned starting values-------*/
    /*-----explicitly are assigned the default value (1)-----*/

    /*-----Compute log likelihood contribution of the  -----*/
    /*-----survival data part when the last observaion -----*/
    /*-----of a subject is reached                     -----*/

    /*-----NOTE: This parameterization yields estimates ----*/
    /*-----equivalent to those in LIFEREG with          ----*/
    /*-----Weibull Distribution   ----*/
    if (LAST_ID) then do;
    linpsurv  = &linpsurv;
    
    %let mname= jm_&jmdist;
    %&mname;

     
    llsurv = (ESRD=1)*log(g) + (ESRD=0)*log(G_t);
    end; else llsurv=0;

    /*---Cholesky parameterization of the random effects ---*/
    /*---variance matrix.                                ---*/
    /*---This ensures that the variance-covariance matrix---*/
    /*---of the random effects is non-negative definite  ---*/
    v11 = a11*a11;
    v12 = a11*a12;
    v22 = a12*a12 + a22*a22;
    SD1=sqrt(v11);
    SD2=sqrt(v22);
    Cor12=v12/(SD1*SD2);

    /*---  Compute the contribution of the longitudinal ---*/
    /*---  part. Every observation in the data set      ---*/
    /*---  makes a contribution. Notice that conditional --*/
    /*---  on the random effects we have independent    ---*/
    /*---  Gaussian contributions                       ---*/
    * linplong = (bl0 + u0)+bl1*SNP&i+(bl2 + u1)*FUTIME+bl3*SNP&i*FUTIME+bl4*Base_gfr; /* Original syntax */
    linplong = &linplong;
    resid = (gfr-linplong);
    if (abs(resid) > 1.3E100) or (S2 < 1e-12) then do;
       lllong = -1e20;
    end; else do;
       lllong = -0.5*(1.837876 + resid**2 / S2  + log(S2));
    end;

    /*---Any numeric variable in the data set can be used---*/
    /*---as the response in the MODEL statement. It has  ---*/
    /*---no bearing on the results                       ---*/
    model LAST_ID ~ general(lllong + llsurv);
    random u0 u1 ~ normal([0, 0],[v11,v12,v22]) subject=_numid;

    /*---  Compute median of the patient-specific      -----*/
    /*---  survival distributions                      -----*/
  

    estimate 'Var[U0]'    v11;
    estimate 'Cov[U0,U1]' v12;
    estimate 'Var[U1]'    v22;
    estimate 'SD1'        SD1;
    estimate 'SD2'        SD2;
    estimate 'Cor12'      Cor12;
run;

%let syserr1  = &syserr;
data _AUX; 
LENGTH SYSERR1 3 ;
LENGTH SNP_NO 8 ;
LENGTH _SNP $16;
SNP_NO=symget("i");
SYSERR1=symget("syserr1");
_SNP=compress("SNP",snp_no);
run;

%*if &i=3 %then %drop_data (data = JM_Est);
%if %sysfunc(exist(JM_Est)) and %sysfunc(exist(JM_Fit)) and %sysfunc(exist(JM_Num))  and %sysfunc(exist(JM_Conv))  and %sysfunc(exist(JM_Ite))%then %do;
%jm_add_main_append;

%end;
%mend jm_model;


%macro jm_add_main (lmm_xvar = , SRV_xvar = , m = 1, drop = T, jmoptions=%str(tech=trureg),jmdist=weibull, init_val=nlmixed_start);

  %let auxm = JM&m;
  %&init_val._main; /* nlmixed_start, nlxd_start, jmhybrid_start */
   /* conditionally execute jm_model depending on whether initial data exist or not  */
  %if  %sysfunc(exist(work.JM_init&m))  %then %jm_model;
  
%mend jm_add_main;


%macro jm_add_main_append;
ods  output Variables=JM_Ite_Var;
proc contents data=JM_Ite;run;

data JM_Ite_Var;
set JM_Ite_Var;
if num=6 then call symput("Ite_var",variable);
run;

data JM_Ite (rename=(&Ite_var=Index));
set JM_Ite;
run;

data JM_Ite;
set JM_Ite end=last;
if last;
run;

data JM_diag1;
 set JM_Est end=_last_;
 if Estimate = "." or Estimate = 0 then JM_ADD&m._DIAG1=1;
 else  JM_ADD&m._DIAG1=0;

 if StandardError="." or StandardError = 0 then JM_ADD&m._DIAG2=1;
else  JM_ADD&m._DIAG2=0;

JM_ADD&m._DIAG_BETA+JM_ADD&m._DIAG1;
JM_ADD&m._DIAG_SE+JM_ADD&m._DIAG2;
if _last_;
KEEP JM_ADD&m._DIAG_BETA JM_ADD&m._DIAG_SE;
run;


%drop_data (data = JM_init&m);

data jm_EST;
set jm_EST end=last;
retain 
JM_ADD&m._LMM_INT_BETA
JM_ADD&m._LMM_INT_SE
JM_ADD&m._LMM_SNP_BETA
JM_ADD&m._LMM_SNP_SE
JM_ADD&m._LMM_SNP_P
JM_ADD&m._LMM_TIME_BETA
JM_ADD&m._LMM_TIME_SE
JM_ADD&m._LMM_TIME_P
JM_ADD&m._LMM_SNP_TIME_BETA
       JM_ADD&m._LMM_SNP_TIME_SE 
       JM_ADD&m._LMM_SNP_TIME_P 
       JM_ADD&m._SRV_SNP_BETA 
       JM_ADD&m._SRV_SNP_SE 
       JM_ADD&m._SRV_SNP_P 
       JM_ADD&m._GAMMA_SIGMA 
       JM_ADD&m._DF;     
       
 if Parameter = "bl0" then do;
              JM_ADD&m._LMM_INT_BETA=Estimate;
              JM_ADD&m._LMM_INT_SE=StandardError;
            
end;

 if Parameter = "bl1" then do;
              JM_ADD&m._LMM_SNP_BETA=Estimate;
              JM_ADD&m._LMM_SNP_SE=StandardError;
              JM_ADD&m._LMM_SNP_P=Probt;
            
end;

 if Parameter = "bl2" then do;
              JM_ADD&m._LMM_TIME_BETA=Estimate;
              JM_ADD&m._LMM_TIME_SE=StandardError;
              JM_ADD&m._LMM_TIME_P=Probt;
            
end;
       
if Parameter = "bl3" then do;
       JM_ADD&m._LMM_SNP_TIME_BETA=Estimate;
       JM_ADD&m._LMM_SNP_TIME_SE=StandardError;
       JM_ADD&m._LMM_SNP_TIME_P=Probt;
end;

if Parameter = "bs1" then do;
       JM_ADD&m._SRV_SNP_BETA=Estimate;
       JM_ADD&m._SRV_SNP_SE=StandardError;
       JM_ADD&m._SRV_SNP_P=Probt;
       JM_ADD&m._DF=DF;
end;

if last then JM_ADD&m._GAMMA_SIGMA =Estimate;
if last;
drop Parameter Estimate StandardError Probt DF;
run;

%append_snp(JM_ADD&m._est, JM_est, &drop);
%append_snp(JM_ADD&m._num, JM_num, &drop);
%append_snp(JM_ADD&m._fit, JM_fit, &drop);
%append_snp(JM_ADD&m._Conv, JM_Conv, &drop);
%append_snp(JM_ADD&m._AUX, _AUX,   &drop);
%append_snp(JM_ADD&m._diag1, JM_diag1, &drop);
%append_snp(JM_ADD&m._Ite, JM_Ite, &drop);
%mend jm_add_main_append;

%macro jm_add_init(lmm_xvar = , SRV_xvar = , m = 1, drop = T, jmoptions=%str(tech=trureg),jmdist=weibull, init_val=nlmixed_start);
%let auxm=JM&m;
%&init_val._init;


data JM_ADD&m._est;

 length JM_ADD&m._LMM_INT_BETA 8;
 length JM_ADD&m._LMM_INT_SE 8;
 
  length JM_ADD&m._LMM_SNP_BETA 8;
  length JM_ADD&m._LMM_SNP_SE 8;
 length JM_ADD&m._LMM_SNP_P 8;
 
  length JM_ADD&m._LMM_TIME_BETA 8;
  length JM_ADD&m._LMM_TIME_SE 8;
 length JM_ADD&m._LMM_TIME_P 8;
 
 length JM_ADD&m._LMM_SNP_TIME_BETA 8;
 length JM_ADD&m._LMM_SNP_TIME_SE 8;
 length JM_ADD&m._LMM_SNP_TIME_P 8;
 length JM_ADD&m._SRV_SNP_BETA 8;
 length JM_ADD&m._SRV_SNP_SE 8;
 length JM_ADD&m._SRV_SNP_P 8;
 length JM_ADD&m._GAMMA_SIGMA 8;
 length JM_ADD&m._DF 8;
 length _SNP $16;
 length SNP_NO 8;
 if _n_<1;
run;

data JM_ADD&m._fit;
 length Descr $40;
 length Value 8;
 length _SNP $16;
 length SNP_NO 8;
 if _n_<1;
run;

data JM_ADD&m._Num;
 length Descr $40;
 length Value 8;
 length _SNP $16;
 length SNP_NO 8;
 if _n_<1;
run;

data JM_ADD&m._Conv;
 length Status 8;
 length _SNP $16;
 length SNP_NO 8;
 if _n_<1;
run;

data JM_ADD&m._AUX;
 length SYSERR1 3;
 length _SNP $16;
 length SNP_NO 8;
 if _n_<1;
run;

data JM_ADD&m._diag1;
length JM_ADD&m._DIAG_BETA 8;
length JM_ADD&m._DIAG_SE 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

data JM_ADD&m._Ite;
length Iter 8;
length Calls 8;
length NegLogLike 8;
length Diff 8;
length MaxGrad 8;
length Index 8;
length _SNP $16;
length SNP_NO 8;
if _n_<1;
run;

%mend jm_add_init;


%macro jm_add_finish(lmm_xvar = , SRV_xvar = , m = 1, drop = T, jmoptions=%str(tech=trureg),jmdist=weibull, init_val=nlmixed_start);
%let auxm = JM&m;
%&init_val._finish;

data JM_ADD&m._Num1(rename=(Value=N));
 set JM_ADD&m._Num; 
 where Descr = "Observations Used";
 drop Descr ;
run;


data JM_ADD&m._fit1;
 set JM_ADD&m._fit;
 retain min_2log;
 if Descr="AIC (smaller is better)" then AIC=Value;
 if  Descr="-2 Log Likelihood" then min_2log=Value;
 if  Descr="AIC (smaller is better)" ;
 drop Descr Value;
run;


data JM_ADD&m;
 merge JM_ADD&m._est JM_ADD&m._fit1 JM_ADD&m._Num1 JM_ADD&m._Ite JM_ADD&m._Conv JM_ADD&m._AUX  JM_ADD&m._diag1 ;
 by SNP_no;
run;

data JM_ADD&m(rename = (    
                         AIC=JM_ADD&m._AIC 
                         Min_2log=JM_ADD&m._MIN2LOG 
                         N=JM_ADD&m._N 
                         STATUS=JM_ADD&m._CONV 
                         SYSERR1=JM_ADD&m._SYSERR
                         Iter=JM_ADD&m._Iter
                         Calls=JM_ADD&m._Calls
                         NegLogLike=JM_ADD&m._NegLogLike
                         Diff=JM_ADD&m._Diff
                         MaxGrad=JM_ADD&m._MaxGrad
                         Index=JM_ADD&m._Index
                        ));

set JM_ADD&m;
drop _snp;


LABEL  
       JM_ADD&m._LMM_INT_BETA  ="Beta Estimate of intercept of LMM part in the Joint Model"
       JM_ADD&m._LMM_INT_SE    ="Standard Error of intercept of LMM part in the Joint Model"
       JM_ADD&m._LMM_INT_P     ="P-Value of intercept of LMM part in the Joint Model"
       
       JM_ADD&m._LMM_SNP_BETA  ="Beta Estimate of SNP of LMM part in the Joint Model"
       JM_ADD&m._LMM_SNP_SE    ="Standard Error of SNP of LMM part in the Joint Model"
       JM_ADD&m._LMM_SNP_P     ="P-Value of SNP of LMM part in the Joint Model"
       
       JM_ADD&m._LMM_TIME_BETA  ="Beta Estimate of FUTIME of LMM part in the Joint Model"
       JM_ADD&m._LMM_TIME_SE    ="Standard Error of FUTIME of LMM part in the Joint Model"
       JM_ADD&m._LMM_TIME_P     ="P-Value of TIME of FUTIME part in the Joint Model"
       
       JM_ADD&m._LMM_SNP_TIME_BETA  ="Beta Estimate of each interaction term (SNP*FUTIME) of LMM part in the Joint Model"
       JM_ADD&m._LMM_SNP_TIME_SE    ="Standard Error of each interaction term (SNP*FUTIME)of LMM part in the Joint Model"
       JM_ADD&m._LMM_SNP_TIME_P     ="P-Value of each interaction term (SNP*FUTIME)of LMM part in the Joint Model"
       JM_ADD&m._SRV_SNP_BETA       ="Beta Estimate of each SNP of Survival Model part in the Joint Model"
       JM_ADD&m._SRV_SNP_SE         ="Standard Error of each SNP of Survival Model part in the Joint Model"
       JM_ADD&m._SRV_SNP_P          ="P-Value of each SNP of Survival Model part in the Joint Model"
       JM_ADD&m._GAMMA_SIGMA        ="Estimate of gamma or sigma of Survival Model part in the Joint Model"
       JM_ADD&m._DF                 ="Degree of freedom  in the Joint Model"
       SNP_NO                       ="Modified SNP order in the final analysis work._geno file"
       AIC                          ="AIC (Akaike information criterion) from fit statistics in the full Joint Model"
       Min_2log                     ="-2 Log Likelihood in the full Joint Model"
       N                            ="Number of observations used in the full Joint Model"
       STATUS                       ="GCONV convergence status in the Joint Model. When status = 1, ERROR = Optimization cannot be completed"
       SYSERR1                      ="System error message in the full Joint Model"
       Iter                         ="The iteration number from Iteration History"
       Calls                        ="The number of function calls from Iteration History"
       NegLogLike                   ="The value of the objective function from Iteration History"
       Diff                         ="The difference between adjacent function values from Iteration History"
       MaxGrad                      ="The maximum of the absolute (projected) gradient components from Iteration History"  
       Index                        ="When tech=trureg, Index=radius. When tech=nrridg, Index=Rho. When tech=quanew or newrap, Index=Slope"
       JM_ADD&m._DIAG_BETA          ="Count of Beta Estimates that are missing. Scale parameter is also included"
       JM_ADD&m._DIAG_SE            ="Count of Standard Errors for Beta Estimates that are missing. Standard Error of scale parameter is also included";
run;

/* When tech=trureg, Index=radius: The radius of the trust region (TRUREG only) from Iteration History.
   When tech=nrridg, Index=Rho: The ratio between the achieved and predicted values of Diff (NRRIDG only)from Iteration History.
   When tech=quanew, Index=Slope: The slope of the search direction at the current parameter iterate (QUANEW only) from Iteration History */
   
proc printto print = mylog; run;
proc means data = JM_ADD&m (drop=snp_no);
Title "Summary of the results from JM_ADD&m for &nsnps SNPs. Optimization technique can be chosen by using tech= in the 'JM_ADD' macro";
Title2 "lmm_xvar := &lmm_xvar (+ SNP + FUTIME + SNP*FUTIME)";
Title3 "srv_xvar := &srv_xvar (+ SNP + u0 +u1) /&jmdist";
Title4 "jmdist := &jmdist, jmoptions : =  &jmoptions";

proc printto; run;

%let merges=&merges jm_add&m;
%mend jm_add_finish;