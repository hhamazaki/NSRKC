//==========================================================================
// 11.0  Tire 3 OFL  Calculation  
//==========================================================================
FUNCTION get_reference_points
  cout<<"start reference points"<<endl;
  int i,j,k,kk;
  dvar_matrix ref_catch_s(1,100,1,na);
  dvar_matrix ref_catch_sd(1,100,1,na);
  dvar_matrix ref_catch_w_c(1,100,1,na);
  dvar_matrix ref_catch_w_s(1,100,1,na);
  dvar_matrix ref_catch_w_cd(1,100,1,na);
  dvar_matrix ref_catch_w_sd(1,100,1,na); 
  dvar_matrix ref_catch(1,100,1,na);
  dvar_matrix ref_ps(1,100,1,na);
  dvar_matrix ref_nps(1,100,1,na);
  dvar_matrix ref_ops(1,100,1,na);
  dvar_matrix ref_npw(1,100,1,na);
  dvar_matrix ref_opw(1,100,1,na); 
  dvar_vector ref_catch_b(1,100);
  dvar_vector ref_catch_n(1,100);
  dvariable pp, wpm, wpim, tt1, ttt0, ttt1, sb, tw;
  dvar_vector nscaf(1,na), tsc(1,na),tscd(1,na), tt0(1,na);    
  dvar_vector sselw(1,na),dselw(1,na),lselw(1,na),dlselw(1,na);
  dvar_vector sselc(1,na),dselc(1,na); 
  dvar_vector nwc(1,na),owc(1,na),nwcd(1,na),owcd(1,na);
  dvar_vector nws(1,na),ows(1,na),nwsd(1,na),owsd(1,na);
  dvar_vector nsc(1,na),osc(1,na),nscd(1,na),oscd(1,na);  
  dvar_vector par1(1,na),par2(1,na),fx(1,na), ref_Hs(1,na),ref_Hw(1,na);   
  dvariable TotalCatch, CatNum, cofl,fl;    
  dvar_matrix ref_ws(1,na,1,4);  
  dvariable ref_F, ref_FSub, ref_FSubd, ref_rec;
  dvariable ref_Hwc;
  dvariable ref_Hts;
  dvar_vector ref_Htw(1,101);
  dvar_vector ref_mbio(1,101);
  dvar_vector ref_mmb(1,100);
  dvar_vector ref_fcatch_s(1,101);
  dvar_vector ref_fcatch_sd(1,101);
  dvar_vector ref_fcatch_w_c(1,101);
  dvar_vector ref_fcatch_w_s(1,101);
  dvar_vector ref_fcatch_w_cd(1,101);
  dvar_vector ref_fcatch_w_sd(1,101);   
  dvar_vector ref_legaln(1,101);
  dvar_vector ref_Hrt(1,100);
  dvar_vector ref_Hr(1,101);
  dvar_vector ref_legal(1,100);
  dvar_vector ref_legal_b(1,100);
  dvariable ref_ys;
  dvariable f35,f40,h35,h40,i35,i40;

  
  ref_Htw.initialize();
  ref_FSub.initialize();
  ref_FSubd.initialize();  
  ref_ys = ys(lyear);
  ref_Hwc = 0.0;
//--------------------------------------------------------------------------
// Determine average subsistence harvest rate (10 years average)  
//--------------------------------------------------------------------------
  for (i = lyear-10; i<= lyear-1; i++)
   {
    ref_FSub += FSub(i);
	ref_FSubd += FSubd(i);
   }
  ref_FSub = ref_FSub/10.0; 
  ref_FSubd = ref_FSubd/10.0;   
   ref_ys = ys(lyear);
// Calculate winter subsistence harvest selectivity
// Assume subsistence harvest take all mature crab caught
   sselw = elem_prod(matc,selw); 
// Assume subsistence harvest discards all immature crab caught   
   dselw = elem_prod(imatc,selw); 
// Calculate legal crab caught by winter gear    
   lselw = elem_prod(lg,selw);   
// Calculate sublegal crab caught by winter gear (discards)  
   dlselw = elem_prod(slg,selw);      
// Calculate legal crab caught by Summer Commercial
   sselc = elem_prod(sel(lyear),lg);   
// Calculate sublegal crab caught by Summer commercial        
   dselc = elem_prod(sel(lyear),slg);  
// Calculate average 5 years of recruit; 
   ref_rec = rec(lyear);


//--------------------------------------------------------------------------
// Simulation for 100 years   
//--------------------------------------------------------------------------
  ref_Hts = 0.0;
  for (kk = 1; kk <= 101; kk++)
   {
// Set Fishery mortality
  ref_F = 0.01*kk-0.01;

//--------------------------------------------------------------------------
// 8.1 Calculate Winter and Summer commercial Harvest rate:   
//--------------------------------------------------------------------------	    	       
  par1 = 1-exp(-ref_F-0.42*Mn);
  par2 = 1-pwh*(1-exp(-0.42*Mn));
//Winter commercial harvest rate   
  ref_Hw = elem_div(pwh*(1-exp(-ref_F))*exp(-0.42*Mn),par2);
// Summer commercial harvest rate  
  ref_Hs = elem_div((1-pwh)*(1-exp(-ref_F))*exp(-0.42*Mn),par2); 
// Initialize temp matrix
  ref_nps.initialize();
  ref_ops.initialize();
  ref_npw.initialize();
  ref_opw.initialize();
// Set First year
  ref_npw(1) = npw(lyear);
  ref_opw(1) = opw(lyear);

//==========================================================================
// 8.7  Start 100 year Simulation 
//==========================================================================
    for (i=1;i<100;i++)
     { 
//--------------------------------------------------------------------------
// 8.1 Winter commercial catch and discards size composition 
//--------------------------------------------------------------------------	    	       
// Calculate number of newshell and oldshell in commercial catch 	  
	nwc = elem_prod(elem_prod(ref_npw(i),lselw),ref_Hw);
	owc = elem_prod(elem_prod(ref_opw(i),lselw),ref_Hw);
// Calculate number of newshell and oldshell discards in commercial  	  
	nwcd = elem_prod(elem_prod(ref_npw(i),dlselw),ref_Hw);
	owcd = elem_prod(elem_prod(ref_opw(i),dlselw),ref_Hw);
    ref_catch_w_c(i) = nwc + owc;
    ref_catch_w_cd(i) = nwcd + owcd;	
//--------------------------------------------------------------------------
// 8.2 Winter subsistence catch and discards size composition 
//     Assume that Subsistence fishers took all length >= 3 crab  
//	   and discarded length 1 & 2 crab
//--------------------------------------------------------------------------
// Calculate number of newshell and oldshell in subsistence catch 
		nws = elem_prod(ref_npw(i),sselw)*ref_FSub;
		ows = elem_prod(ref_opw(i),sselw)*ref_FSub;
// Calculate number of newshell and oldshell in subsistence discards  		
		nwsd = elem_prod(ref_npw(i),dselw)*ref_FSubd;
        owsd = elem_prod(ref_opw(i),dselw)*ref_FSubd;
    ref_catch_w_s(i) = nws + ows;
    ref_catch_w_sd(i) = nwsd + owsd;
//--------------------------------------------------------------------------
// 8.3 Calculate abundance for summer July 1st:
//     Summer population is survivors of winter fisheries 
//     and natural mortality
//--------------------------------------------------------------------------
  ref_nps(i) = elem_prod(ref_npw(i)-nwc-nws-nwcd*hm(2)-nwsd*hm(2),exp(-0.417*Mn));  
  ref_ops(i) = elem_prod(ref_opw(i)-owc-ows-owcd*hm(2)-owsd*hm(2),exp(-0.417*Mn));  
		 
//--------------------------------------------------------------------------
// 8.5  Calculate Summer Commercial catch  
//--------------------------------------------------------------------------    
// Calculate Summer commercial legal and sublegal catch by length;  
	nsc = elem_prod(elem_prod(ref_nps(i),sselc),ref_Hs);
	osc = elem_prod(elem_prod(ref_ops(i),sselc),ref_Hs); 
// Calculate number of newshell and oldshell discards in commercial  	  
	nscd = elem_prod(elem_prod(ref_nps(i),dselc),ref_Hs);
	oscd = elem_prod(elem_prod(ref_ops(i),dselc),ref_Hs);
// Total summer commercial catch 	  
	ref_catch_s(i) = nsc+osc;  	
	ref_catch_sd(i) = nscd+oscd;

//--------------------------------------------------------------------------
// 8.6  Calculate number of Summer Crab after summer fishery (nscaf)
//      nscaf = July 1st abundance*mortality till fishery 
//               - com catch and discards. 
//--------------------------------------------------------------------------
  nscaf = elem_prod(ref_nps(i)+ref_ops(i),exp(-ref_ys*Mn))-ref_catch_s(i)-ref_catch_sd(i)*hm(1);	  
  	
// Stop gap measure: abundance does not go below zero
    for (k=1;k<=na;k++)
     {	  
      if (nscaf(k) < 0.0) nscaf(k)= 0.001;  
	  }
	ref_Hrt(i) = sum(nscaf)/sum(ref_nps(i)+ref_ops(i));  
//--------------------------------------------------------------------------
// 8.7  Calculate Crab abundance on Feb 1st 
//      Newshell: molted*mortality + recruit
//      Oldshell: unmolted*mortality
//--------------------------------------------------------------------------    
// Calculate New Shell population abundance by length class:
    for (j=1;j<=na;j++)
     {
      pp = 0.0; 
	  //Each crab molts right after fishery ended  
      for (k=1;k<=j;k++) pp += tgr(k,j)*nscaf(k)*mp1(k); 	  
	  //New shell crab molted*Mortality + recruits
      ref_npw(i+1,j) = pp*exp(-(0.583-ref_ys)*Mn(j)) + rpp(j)*ref_rec;    
	  //Old shell crab are unmolted 
      ref_opw(i+1,j) = nscaf(j)*(1.0-mp1(j))*exp(-(0.583-ref_ys)*Mn(j));    
     }
//==========================================================================
// 8.7  Calculate February 1st legal Crab and mature male biomass  for assessment
//==========================================================================

  ref_mmb(i) = sum(elem_prod(elem_prod(ref_npw(i)+ref_opw(i),matc),wm));  // Mature male biomass
  ref_legal(i) = sum(elem_prod(elem_prod(ref_npw(i)+ref_opw(i),sselc),wm));  // number of legal crab  
   }	  

//==========================================================================
// Extract Harvest and biomass
//==========================================================================
    ref_mbio(kk) = ref_mmb(99);                    //lbs/R
//    ref_totc(j) = ref_catch_b(100);
    ref_legaln(kk) = ref_legal(99);
	ref_fcatch_s(kk) = sum(ref_catch_s(99));
    ref_fcatch_sd(kk) = sum(ref_catch_sd(99));
    ref_fcatch_w_c(kk) = sum(ref_catch_w_c(99));
    ref_fcatch_w_s(kk) = sum(ref_catch_w_s(99));
    ref_fcatch_w_cd(kk) = sum(ref_catch_w_cd(99));
    ref_fcatch_w_sd(kk) = sum(ref_catch_w_sd(99));  
    ref_Hr(kk) = ref_Hrt(99);	 
//    cout << "ncrab" << ref_opw(99)+ref_npw(99)<<endl; 
  }

//==========================================================================
// 9.0  Find F35%
//==========================================================================
   
   i35 = 0; i40 = 0;
  
  for (j = 1; j<= 101; j++)
   {
    if (i35 < 1.0)
// Find F35 F35 = 35% of reference mmb		
     if (ref_mbio(j) <= 0.35*ref_mbio(1))
      {
       f35 = 0.01*j-0.01;
//     h35 = ref_Hr(j);
       i35 = 2.0;
      }
    if (i40 < 1.0)
     if (ref_mbio(j) <= 0.40*ref_mbio(1))
      {
       f40 = 0.01*j-0.01;
//       h40 = ref_Hr(j);
       i40 = 2.0;
      }
   }

  ofstream report1("refp.out");
  report1 <<"Harvest rate in term of legal males as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_Hr<<endl;
  report1 <<"Total MMB (>93mm) Feb 01 as F = 0.00, 0.01, ... 1.0 (lbs/R)"<<endl;
  report1 << ref_mbio<<endl;
  report1 <<"Total legal as F = 0.00, 0.01, ... 1.0 (lbs/R)"<<endl;
  report1 << ref_legaln<<endl;
  report1 <<"Total Summer Com catch as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_fcatch_s<<endl;
  report1 <<"Total Summer dis as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_fcatch_sd<<endl;
  report1 <<"Total Winter Com catch as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_fcatch_w_c<<endl;
  report1 <<"Total Winter Com dis as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_fcatch_w_cd<<endl;
  report1 <<"Total Winter Sub catch as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_fcatch_w_s<<endl;
  report1 <<"Total Winter Sub dis as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_fcatch_w_sd<<endl;  
  report1 <<"Total legals as F = 0.00, 0.01, ... 1.0 (crabs/R)"<<endl;
  report1 << ref_legaln<<endl;
  report1 <<"F35: "<<f35<<"  H35%: "<<h35<<endl;
  report1 <<"F40: "<<f40<<"  H40%: "<<h40<<endl;
//  report1 <<"ref_Hwc = (mean of 10 years):  "<<ref_Hwc<<endl;
  cout<<"end of reference points"<<endl;

//==========================================================================

