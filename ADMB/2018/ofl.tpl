FUNCTION get_reference_points
// cout<<"start reference points"<<endl;
  dvariable mr;
  ref_Ft = 0.0;
  ref_Ff = 0.0;
  ref_Ftc = 0.0;
  ref_sel_fit.initialize();
  ref_sel_ret.initialize();
  ref_moltp.initialize();
  ref_sel_trawl.initialize();

  mr = 0;
  for(i=ryear(1); i<=ryear(2); i++)
  {
     mr += mfexp(mean_log_rec+rec_dev(1,i))+mfexp(mean_log_rec+rec_dev(2,i));
  }
  mr = mr/((ryear(2)-ryear(1)+1.0)*2.0);           //mean recruitment for B35 estimation
  for (i = endyr-9; i<=endyr; i++)
  {
    ref_Ft += fmortt(i);
    ref_Ff += fmortf(i);
    ref_Ftc += Ftcm(i);
    for (k=1; k<=nlenm; k++)
    {
       ref_sel_trawl(k) += sel_trawl_m(i,k);
       ref_sel_fix(k) += sel_fix_m(i,k);
       ref_moltp(k) += moltp(i-1,k);
    }
  }
  for (i = endyr-2; i<=endyr-1; i++)
  {
    for (k=1; k<=nlenm; k++)
    {
      ref_sel_fit(k) += sel_fit(i,k);
      ref_sel_ret(k) += sel(i,k)-sel_fit(i,k);
    }
  }
  ref_Ft = ref_Ft/10.0;
  ref_Ff = ref_Ff/10.0;
  ref_Ftc = ref_Ftc/10.0;
  ref_sel_trawl = ref_sel_trawl/10.0;
  ref_sel_fix = ref_sel_fix/10.0;
  ref_moltp = ref_moltp/10.0;
  ref_sel_fit = ref_sel_fit/2.0;
  ref_sel_ret = ref_sel_ret/2.0;

  for (j = 1; j<= 101; j++)
 {
   ref_F = 0.01*j-0.01;
   ref_Fret=ref_sel_fit*ref_F;
   ref_S = mfexp(-1.0*(ref_Fret+ref_sel_ret*ref_F));
   tem_S = ref_S;
   ref_Ftot = ref_Ftc*sel_tcm+ref_sel_trawl*ref_Ft + ref_sel_fix*ref_Ff;
   if (j==1) ref_Ftot = 0.0000000000001;
   tem_Stc = mfexp(-1.0*ref_Ftot);
   if (j==1)
   {
     tem_S = 1.0;
     tem_Stc = 1.0;
   }
   ref_catch.initialize();
   ref_catch_ret.initialize();
   ref_catch_m.initialize();
   ref_na_m.initialize();
   ref_na.initialize();
// cout<<ref_catch<<endl;
//initial year
   ref_na_m(1,1) = na_m(1,endyr);
   ref_na_m(2,1) = na_m(2,endyr);
   ref_na(1)     = ref_na_m(1,1)  + ref_na_m(2,1);

// Now do Recruitments..........
   for (i=2;i<=100;i++)
   {
      ref_na_m(1,i) += 1000000.0 *rec_len(2);
   }
//numbers at length
   for (i=1;i<100;i++)
   {
     dvar_vector ref_tmp = elem_prod(ref_moltp*mfexp(-(1.0-tc_cm(endyr))*M(2,endyr)),elem_prod(mfexp(-(tc_cm(endyr)-cm(endyr))*M(2,endyr))*elem_prod(tem_S,mfexp(-cm(endyr)*M(2,endyr))*ref_na(i)),tem_Stc));
     ref_na_m(1,i+1) +=  ref_tmp * len_len(4);
     ref_na_m(2,i+1) = elem_prod((1.0-ref_moltp)*mfexp(-(1-tc_cm(endyr))*M(2,endyr)),elem_prod(mfexp(-(tc_cm(endyr)-cm(endyr))*M(2,endyr))*elem_prod(tem_S,mfexp(-cm(endyr)*M(2,endyr))*ref_na(i)),tem_Stc));
     ref_na(i+1)     = ref_na_m(1,i+1)  + ref_na_m(2,i+1);
     ref_na_fishtime(i) = mfexp(-0.30*M(2,endyr))*ref_na(i);
     ref_na_fishtime_tc(i) = elem_prod(ref_na_fishtime(i),tem_S)*mfexp(-(tc_cm(endyr)-cm(endyr))*M(2,endyr));
     ref_mbio215(i) = (elem_prod(mfexp(-(0.694-tc_cm(endyr))*M(2,endyr))*mat,elem_prod(mfexp(-(tc_cm(endyr)-cm(endyr))*M(2,endyr))*elem_prod(tem_S,mfexp(-cm(endyr)*M(2,endyr))*ref_na(i)),tem_Stc)))*wt(2);
    for (k = 1; k<= nlenm; k++)
    {
       ref_catch_lmale(i,k) = ref_na_fishtime(i,k)*(1.0-ref_S(k));
       ref_catch(i) += ref_catch_lmale(i,k)*wt(2,k);
       ref_catch_male_ret(i,k) = ref_na_fishtime(i,k)*(1.0-mfexp(-1.0*ref_Fret(k)));
       ref_catch_ret(i) += ref_catch_male_ret(i,k)*wt(2,k);
    }
    ref_catch_m(i) = ref_catch(i) - ref_catch_ret(i);
     ref_catch_trawl(i)= ref_na_fishtime_tc(i)*elem_prod(elem_prod(1.0-exp(-1.0*ref_Ftot),elem_div(ref_sel_trawl*ref_Ft,ref_Ftot)),wt(2));
     ref_catch_fix(i)  = ref_na_fishtime_tc(i)*elem_prod(elem_prod(1.0-exp(-1.0*ref_Ftot),elem_div(ref_sel_fix*ref_Ff,ref_Ftot)),wt(2));
   }

   ref_na_fishtime(99) = mfexp(-0.30*M(2,endyr))*ref_na(99);
   ref_mbio(j) = ref_mbio215(99)/1000.0;                    //kg/R
   ref_totc(j) = ref_catch(99)/1000.0;
   ref_retc(j) = ref_catch_ret(99)/1000.0;
 }
    i35 = 0;
    i40 = 0;
  for (j = 1; j<= 101; j++)
 {
    if (i35 < 1.0)
    {
      if (ref_mbio(j) <= 0.35*ref_mbio(1))
      {
         f35 = 0.01*j-0.01;
         eb35 = ref_mbio(j)*mr/1000.0;
         b25 = 0.25*eb35;
         i35 = 2.0;
      }
    }
    if (i40 < 1.0)
    {
      if (ref_mbio(j) <= 0.40*ref_mbio(1))
      {
         f40 = 0.01*j-0.01;
         i40 = 2.0;
      }
    }
 }
  ref_f = 0.0;
  ref_sel_f = 0.0;
  for (i = endyr-2; i<=endyr-1; i++)
  {
    ref_f += fmortdf(i);
    for (k=1; k<=nlenm; k++)
    {
      ref_sel_f(k) += sel_discf(i,k);
    }
  }
  ref_sel_f = ref_sel_f/2.0;
  ref_f = ref_f/2.0;
  ofl_f = f35;
  ref_Ftot = ref_Ftc*sel_tcm+ref_sel_trawl*ref_Ft + ref_sel_fix*ref_Ff;
  ref_S = mfexp(-1.0*(ref_sel_fit*ofl_f+ref_sel_ret*ofl_f));
  tem_Stc = mfexp(-1.0*ref_Ftot);
  if (last_mmb < eb35)
  {
     for (k = 1; k<10; k++)
     {
        ref_Ftot = ref_Ftc*sel_tcm+ref_sel_trawl*ref_Ft + ref_sel_fix*ref_Ff;
        ref_S = mfexp(-1.0*(ref_sel_fit*ofl_f+ref_sel_ret*ofl_f));
        last_mmb = (elem_prod(mfexp(-(0.694-tc_cm(endyr))*M(2,endyr))*mat,elem_prod(mfexp(-(tc_cm(endyr)-cm(endyr))*M(2,endyr))*elem_prod(ref_S,mfexp(-cm(endyr)*M(2,endyr))*na(endyr)),tem_Stc)))*wt(2);
        if (last_mmb < b25)
        {
           ofl_f = 0.0;
        }
        else
        {
            ofl_f = f35*(last_mmb/eb35-0.1)/0.9;
        }
     }
  }
   ref_Ftot = ref_Ftc*sel_tcm+ref_sel_trawl*ref_Ft+ref_sel_fix*ref_Ff;
   ref_Fret=ref_sel_fit*ofl_f;
   ref_S = mfexp(-1.0*(ref_sel_fit*ofl_f+ref_sel_ret*ofl_f));
   ref_f = ref_f * ofl_f/f35;
   ref_Ftotf = ref_Ftc*sel_tcf+ ref_sel_trawl*ref_Ft+ref_sel_fix*ref_Ff;
   ref_catch(1) = 0.0;
   ref_catch_ret(1) = 0.0;
   for (k = 1; k<= nlenm; k++)
   {
      ref_catch_lmale(1,k) = na_fishtime(2,endyr,k)*(1.0-ref_S(k));
      ref_catch(1) += ref_catch_lmale(1,k)*wt(2,k);
      ref_catch_male_ret(1,k) =  na_fishtime(2,endyr,k)*(1.0-mfexp(-1.0*ref_Fret(k)));
      ref_catch_ret(1) += ref_catch_male_ret(1,k)*wt(2,k);
   }
   ref_catch_trawl(1)= na_fishtime_tc(2,endyr)*elem_prod(elem_prod(1.0-exp(-1.0*ref_Ftot),elem_div(ref_sel_trawl*ref_Ft,ref_Ftot)),wt(2));
   ref_catch_fix(1)  = na_fishtime_tc(2,endyr)*elem_prod(elem_prod(1.0-exp(-1.0*ref_Ftot),elem_div(ref_sel_fix*ref_Ff,ref_Ftot)),wt(2));
   ref_catch_m(1) = ref_catch(1) - ref_catch_ret(1);
   ref_catch_trawl(1) += na_fishtime_tc(1,endyr)*elem_prod(elem_prod(1.0-exp(-1.0*ref_Ftotf),elem_div(ref_sel_trawl*ref_Ft,ref_Ftotf)),wt(1));
   ref_catch_fix(1) += na_fishtime_tc(1,endyr)*elem_prod(elem_prod(1.0-exp(-1.0*ref_Ftotf),elem_div(ref_sel_fix*ref_Ff,ref_Ftotf)),wt(1));
   ref_catch_fix(2) = na_fishtime_tc(1,endyr)*elem_prod(elem_prod(1.0-exp(-1.0*ref_Ftotf),elem_div(sel_tcf*ref_Ftc,ref_Ftotf)),wt(1));
   ref_catch_fix(2) += na_fishtime_tc(2,endyr)*elem_prod(elem_prod(1.0-exp(-1.0*ref_Ftot),elem_div(sel_tcm*ref_Ftc,ref_Ftot)),wt(2));
   ref_catch_f = na_fishtime(1,endyr)*elem_prod((1.0-exp(-1.0*ref_sel_f*ref_f)),wt(1));
   ref_catch(1) = ref_catch_ret(1)+ref_catch_m(1)+ref_catch_f + ref_catch_trawl(1) + ref_catch_fix(1)+ ref_catch_fix(2);
   ofl_catch = ref_catch(1);

  ofstream report1("refp75172b.out");

  report1 <<eb35<<"  "<<last_mmb<<"  "<<ref_catch(1)<<endl;
  report1 <<"Total MMB (>119mm) 2/15 as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_mbio<<endl;
  report1 <<"Total MMB (>129mm) 2/15 as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_mbio1<<endl;
  report1 <<"Total MMB (>109mm) 2/15 as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_mbio2<<endl;
  report1 <<"Total catch as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_totc<<endl;
  report1 <<"Retained catch as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_retc<<endl;
  report1 <<"F35: "<<f35<<"  B35: "<<eb35<<" Mean R: "<<mr<<endl;
  report1 <<"F40: "<<f40<<endl;
  report1 <<"ref_Ft = (mean of 10 years)  "<<endl;
  report1 <<ref_Ft<<endl;
  report1 <<"ref_Ff = (mean of 10 years)  "<<endl;
  report1 <<ref_Ff<<endl;
  report1 <<"ref_sel_trawl = (mean of 10 years)  "<<endl;
  report1 <<ref_sel_trawl<<endl;
  report1 <<"ref_sel_fix = (mean of 10 years)  "<<endl;
  report1 <<ref_sel_fix<<endl;
  report1 <<"ref_sel_retained = (mean of 2 years)  "<<endl;
  report1 <<ref_sel_fit<<endl;
  report1 <<"ref_sel_discarded = (mean of 2 years)  "<<endl;
  report1 <<ref_sel_ret<<endl;
// cout<<"end of reference points"<<endl;
  report1 <<"ref_sel_disc_females = (mean of 2 years)  "<<endl;
  report1 <<ref_sel_f<<endl;
   report1 <<"OFL =  "<<ref_catch(1)<<endl;
   report1 <<"Retained catch =  "<<ref_catch_ret(1)<<endl;
   report1 <<"Male pot bycatch =  "<<ref_catch_m(1)<<endl;
   report1 <<"Female pot bycatch =  "<<ref_catch_f<<endl;
   report1 <<"Trawl bycatch =  "<<ref_catch_trawl(1)<<endl;
   report1 <<"Fixed gear bycatch =  "<<ref_catch_fix(1)<<endl;
   report1 <<"OFL F =  "<<ofl_f<<endl;
   report1 <<"MMB at terminal year =  "<<last_mmb<<endl;
