// OFL Tier 3 Calculation See Andre's method back in 2011. 

// Extract Winter Sub, Com exploitation rate
    FSum(i) = tc(i)/tb(i);
    FSub(i) = tsc(i)/pp;
    FWin(i) = twc(i)/tt1;
	
  if (mceval_phase())
   {
    mcmc_count += 1;
    if (mcmc_count > BurnIn & mcmc_count <= LastSim)
     {
      FratSub = (FSub(ny-5)+FSub(ny-4)+FSub(ny-3)+FSub(ny-2)+FSub(ny-1))/(FSum(ny-5)+FSum(ny-4)+FSum(ny-3)+FSum(ny-2)+FSum(ny-1));
      FratWin = (FWin(ny-5)+FWin(ny-4)+FWin(ny-3)+FWin(ny-2)+FWin(ny-1))/(FSum(ny-5)+FSum(ny-4)+FSum(ny-3)+FSum(ny-2)+FSum(ny-1));
      if (ProjType == 1)
       {
        YrPass = ny;
        Find_OFL();
       }
//      else
//--------------------------------------------------------------------------
//	  {
//        SR_rel = 1;
//        Fmult = 0;
//        PrintDiag = 0;
//        YrPass = ny;
//        ProjConstF();
//        SBPR = mmb(ny+Nproj-1)/rec(ny+Nproj-1);

        // Specify steepness as given above
//        Get_Steepness();

        // Real stock-recruitment relationship
//        SR_rel = SR_FORM;

        // Now do all those projections
//        Nproj = NprojSim;
//        if (R0 > 0)
//        {
//          mcmc2 << Steep << " " << MMB0 << " " << BmsyProx << endl;
//          mcmc2 << " " << -1 << " " << MMB0 << " 1 0 0 0 " << R0 << endl;
//          for (i=1;i<=ny-1;i++)
//           mcmc2 << " " << FirstYr+i-1 << " " << mmb(i) << " " << mmb(i)/MMB0 << " -1 -1 " << tc(i)+twc(i)+tsc(i) << " " << rec(i) << endl;
//          ProjAll();
//         }
//       }
//--------------------------------------------------------------------------	   
     }
   }


// =========================================================================
//  Ignore this section 
// =========================================================================
FUNCTION Get_Steepness
  dvariable SBPRF0, SBPRFmsy, Fmsy, Mbar, RbarFmsy, nn, Term1, Term2;
  dvariable BioDep;
  dvariable Cat1, Cat2;
  dvariable MinSteep,MaxSteep,DerivMin,DerivMax,Deriv;
  int II,iy;

  // Set FMSY
  Fmsy = M;

  // Set projection length
  Nproj = NprojFmsy;

  SR_rel = 1;
  Fmult = 0;
  PrintDiag = 0;
  YrPass = ny;
  ProjConstF();
  SBPRF0 = mmb(ny+Nproj-1)/rec(ny+Nproj-1);
  MMB0 = mmb(ny+Nproj-1);
  R0 = rec(ny+Nproj-1);
  Fmult = Fmsy;
  ProjConstF();
  SBPRFmsy = mmb(ny+Nproj-1)/rec(ny+Nproj-1);

  // Find mean biomass and recruitment at BMSY
  Mbar = 0; nn= 0;
  for (iy=Bmsy_yr1;iy<=Bmsy_yr2;iy++)
   { Mbar += mmb(iy); nn+= 1; }
  BmsyProx = Mbar/nn;
  RbarFmsy = BmsyProx/SBPRFmsy;

  MMB0 = BmsyProx / 0.35;
  R0 = MMB0 / SBPRF0;

  // Projections from here are based on a stock-recruitment relationship
  SR_rel = 2;

  DerivMin = -1.0e20; DerivMax = 1.0e20;
  MinSteep = 0.21; MaxSteep = 0.99;
  for (II=21;II<=99;II++)
   {
    Steep = float(II)*0.01;;
    Fmult = Fmsy + 0.001;
    ProjConstF();
    Cat1 = PassCatch;
    Fmult = Fmsy - 0.001;
    ProjConstF();
    Cat2 = PassCatch;
    Deriv = (Cat1-Cat2)/0.002;
    if (Deriv < 0 & Deriv > DerivMin)
     { MinSteep = Steep; DerivMin = Deriv; }
    if (Deriv > 0 & Deriv < DerivMax)
     { MaxSteep = Steep; DerivMax = Deriv; }
   }

  // Solve for FMSY (fine search)
  for (II=1;II<=40;II++)
   {
    Steep = value((MinSteep+MaxSteep)/2.0);
    Fmult = Fmsy + 0.001;
    ProjConstF();
    Cat1 = PassCatch;
    Fmult = Fmsy - 0.001;
    ProjConstF();
    Cat2 = PassCatch;
    Deriv = (Cat1-Cat2)/0.001;
    if (Deriv < 0) MinSteep = Steep; else MaxSteep = Steep;
   }
  Fmsy = Fmult;
  ProjConstF();
  BioDep = mmb(ny+Nproj-1)/MMB0;
  cout << "Final " << Steep << " " << R0 << " " << MMB0 << " " << BioDep << " " << Deriv << endl;

  // Now compyte B0 given R at BMSY;
  Term1 = (1-Steep) + (5*Steep-1)*BioDep;
  Term2 = 4*Steep*BioDep;
  R0 = RbarFmsy*Term1/Term2;
  MMB0 = R0*SBPRF0;

  for(II=0;II<=141;II++)
   {
    Fmult = (float(II)/21.0)*Fmsy;
//    ProjConstF();
//    cout << "T " << Fmult/Fmsy << " " << PassCatch << " " << mmb(ny+Nproj-1)/MMB0 << " "  << mmb(ny+Nproj-1)/BmsyProx << endl;
   }
// =========================================================================


FUNCTION Find_OFL
  dvariable Bmsy,Fmsy,par1,par2, 
  int ii,iy;

// Define Bmsy:  Bmsy is average of mmb
  Bmsy = 0;
  for (iy=Bmsy_yr1;iy<=Bmsy_yr2;iy++) Bmsy += mmb(iy);
  Bmsy /= float(Bmsy_yr2-Bmsy_yr1+1);

// Specify recruitment: Recruitment is average of past 5 years of recruitment
  RecPass = 0;
  for (iy=YrPass-5;iy<=YrPass-1;iy++) RecPass += rec(iy);
  RecPass /= 5;
// Fmsy is natural mortality
  Fmsy = M;
// Fmult is natural mortality
  par1 = 1-exp(-Fmsy-0.42*Mn);
  par2 = 1-pwh*(1-exp(-0.42*Mn));
  fx = elem_div(1-pwh*par1,par2); 
  Fw = 1-fx;
  Fs = elem_prod(fx-exp(-Fmsy),exp(-0.42*Mn)); 
  
// Start Projection   
  Proj1Yr();
//   
  if (mmb(YrPass+1) < Bmsy)
   {
    Fmult = 0;
    Proj1Yr();
    if (mmb(YrPass+1) > beta*Bmsy)
     {
      Fmult = M/2;
      Proj1Yr();
      for (ii=1;ii<=10;ii++)
       {
        Fmult = Fmsy*(mmb(YrPass+1)/Bmsy-alpha)/(1-alpha);
        Proj1Yr();
       }
     }
   }
 mcmc1 << mcmc_count << " " << Bmsy << " " << mmb(YrPass+1) << " " << Fmult << " " << PassCatch << endl;

//--------------------------------------------------------------------------
// Function Proj1Yr
// Project 1 year 
//--------------------------------------------------------------------------
FUNCTION Proj1Yr
  int Yr,Len,Len1;                               // Counter
  dvariable Total;                               // Totals
  dvar_vector ensc(1,na),eosc(1,na);             // Proportions for SUMMER fishery
  dvar_vector enws(1,na),eows(1,na);             // Proportions for SUBSISTENCE fishery
  dvar_vector enwc(1,na),eowc(1,na);             // Proportions for WINTER fishery
  dvar_vector tt0(1,na);                         // Abundance after summer fishery
  dvariable TotalCatch, CatNum;                  // Total catch (in mass)
  dvariable pp,tt1;
  dvariable yss;                                 // Time to summer fishery

 // Selection is for last yr as is time to fishery (as is the timing of the fishery)
  yss = ys(ny);

  // Catch in weight
  TotalCatch = 0;

// proportions for Winter  fishery (and the catch)
  Total = 0.0;
  
//==========================================================================
// 8.1 Winter commercial catch and discards size composition 
//==========================================================================
	    	       
// Calculate Total Winter Legal Crab Catchable to commercial fisheries 
    tw = sum(elem_prod(npw(YrPass)+opw(YrPass),lselw));	
// Calculate Total Winter discards to commercial fisheries 
    twd = sum(elem_prod(npw(YrPass)+opw(YrPass),dlselw));	 
// Calculate proportion of newshell and oldshell in commercial catch   	  
    enwc(YrPass) = elem_prod(npw(YrPass),lselw)/tw;
    eowc(YrPass) = elem_prod(opw(YrPass),lselw)/tw;
// Calculate number of newshell and oldshell in commercial catch 	  
    CatWin = sum(elem_prod(elem_prod(npw(YrPass),lselw),Fw)
	       + elem_prod(elem_prod(opw(YrPass),lselw),Fw);
	nwc = CatWin*enwc;
	owc = CatWin*eowc; 
// Calculate proportion of newshell and oldshell in commercial catch   	  
    enwc_d(YrPass) = elem_prod(npw(YrPass),dlselw)/twd;
    eowc_d(YrPass) = elem_prod(opw(YrPass),dlselw)/twd;
// Calculate number of newshell and oldshell discards in commercial  	  
    CatWin_d = sum(elem_prod(elem_prod(npw(YrPass),dlselw),Fw)
	         + elem_prod(elem_prod(opw(YrPass),dlselw),Fw));
	nwd = CatWin_d*enwc_d;
	owd = CatWin_d*eowc_d;

//==========================================================================
// 8.2 Winter subsistence catch and discards size composition 
//     Assume that Subsistence fishers took all length >= 3 crab  
//	   and discarded length 1 & 2 crab
//==========================================================================
// Calculate total mature population catchable to winter fisheries 
    wpm = sum(elem_prod(npw(YrPass)+opw(YrPass),sselw));
// Calculate total immature population cathable to winter sfisheries 
    wpim = sum(elem_prod(npw(YrPass)+opw(YrPass),dselw));
// Calculate proportion of newshell and oldshell in winter subsistence  
	enwp(YrPass) = elem_prod(npw(YrPass),sselw)/wpm;
    eowp(YrPass) = elem_prod(opw(YrPass),sselw)/wpm;
// Calculate proportion of newshell and oldshell in winter subsistence  
	enwimp(YrPass) = elem_prod(npw(YrPass),dselw)/wpim;
    eowimp(YrPass) = elem_prod(opw(YrPass),dselw)/wpim;
// Calculate number of newshell and oldshell in winter subsistence 	  
	CatSub = sum(elem_prod(elem_prod(npw(YrPass),sselw),Fw)*FratSub
			 + elem_prod(elem_prod(opw(YrPass),sselw),Fw)*FratSub);
// Calculate number of newshell and oldshell in winter subsistence 	  
	CatSub_d = sum(elem_prod(elem_prod(npw(YrPass),dselw),Fw)*FratSub
			 + elem_prod(elem_prod(opw(YrPass),dselw),Fw)*FratSub);	 
// Calculate proportion of new and old shell discards in winter subsistence  
		nws = CatSub*enwp;
		ows = CatSub*eowp;
		nwsd = CaSub_d*enwimp;
        owsd = CaSub_d*eowimp;

//==========================================================================
// 8.4 Calculate abundance for summer July 1st:
//     Summer population is survivors of winter fisheries 
//     and natural mortality
//==========================================================================
  nps(YrPass) = elem_prod(npw(YrPass)-nwc-nws-nwcd*hm(2)-nwsd*hm(2),exp(-0.417*Mn));  
  ops(YrPass) = elem_prod(opw(YrPaSS)-owc-ows-owcd*hm(2)-owsd*hm(2),exp(-0.417*Mn));  

//==========================================================================
// 8.4  Calculate Commercial catch  
//==========================================================================
// Calculate summer crab legal and sublegal crab selectivity
   sselc = elem_prod(sel(ny),lg);   
   dselc = elem_prod(sel(ny),slg);     
// Total number of crab catchable to commercial fishery    
   tb(YrPass) = sum(elem_prod(nps(YrPass)+ops(YrPass),sselc));	
   tbd(YrPass) = sum(elem_prod(nps(YrPass)+ops(YrPass),dselc));	

// Calculate proportion of newshell and oldshell in commercial		
	enc(YrPass) = elem_prod(nps(YrPass),sselc)/tb(YrPass);
	eoc(Yrpass) = elem_prod(ops(YrPass),sselc)/tb(YrPass);

// Calculate Summer commercial legal and sublegal catch by length;  
    CatSum = sum(elem_prod(elem_prod(nps(YrPass),sselc),Fs)
	       + elem_prod(elem_prod(ops(YrPass),sselec),Fs);
	nsc = CatSum*enc;
	osc = CatSum*eoc; 
// Calculate proportion of newshell and oldshell in commercial catch   	  
    enc_d(YrPass) = elem_prod(nps(YrPass),dselc)/twd;
    eoc_d(YrPass) = elem_prod(ops(YrPass),dselc)/twd;
// Calculate number of newshell and oldshell discards in commercial  	  
    CatSum_d = sum(elem_prod(elem_prod(npw(YrPass),dlselw),Fw)
	         + elem_prod(elem_prod(opw(YrPass),dlselw),Fw));
	nscd = CatSum_d*enc_d;
	oscd = CatSum_d*eoc_d;

    tsc = nsc+osc;
    tscd = nscd+oscd;

//==========================================================================
// 8.5  Calculate number of Summer Crab after summer fishery (nscaf)
//      nscaf = July 1st abundance*mortality till fishery 
//               - com catch and discards. 
//==========================================================================
  		
    nscaf = elem_prod(nps(YrPass)+ops(YrPass),exp(-ys(i)*Mn))-tsc-hm(1)*tscd;	  
	
// Stop gap measure: abundance does not go below zero
    for (k=1;k<=na;k++)
     {	  
      if (nscaf(k) < 0.0) nscaf(k)= 0.001;  
	  }
 
 
//==========================================================================
// 8.6  Calculate Crab abundance on Feb 1st 
//      Newshell: molted*mortality + recruit
//      Oldshell: numolted*mortality
//==========================================================================
     
// Calculate New Shell popululation abundance by length class:
    for (j=1;j<=na;j++)
     {
      pp = 0.0; 
	  //Each crab molts right after fishery ended  
      for (k=1;k<=j;k++) pp += tgr(k,j)*nscaf(k)*molp(i,k); 	  
	  //New shell crab molted + recruits
      npw(YrPass+1,j) = pp*exp(-(0.583-ys(i))*Mn(j)) + rpp(j)*rec;    
	  //Old shell crab are unmolted 
      opw(YrPass+1,j) = nscaf(j)*(1.0-molp(YrPass,j))*exp(-(0.583-yss)*Mn(j));    
     }
	
//==========================================================================
// 8.7  Calculate February 1st legal Crab and mature male biomass  for assessment
//==========================================================================

    mmb(YrPass+1) = sum(elem_prod(elem_prod(npw(YrPass+1)+opw(YrPass+1),matc),wm));  // Mature male biomass
    legaln(YrPass+1) = sum(elem_prod(npw(YrPass+1)+opw(YrPass+1),lg));  // number of legal crab
    legalb(YrPass+1) = sum(elem_prod(npw(YrPass+1)+opw(YrPass+1),wmlg)); // legal biomass   

  PassCatch = TotalCatch;


	
// ===========================================================================

FUNCTION Find_OFL_EST
  dvariable BmsyEst,Fmsy;
  int ii,iy;

  // Define Bmsy
  BmsyEst = 0;
  for (iy=Bmsy_yr1;iy<=Bmsy_yr2;iy++) BmsyEst += mmbEst(iy);
  BmsyEst /= float(Bmsy_yr2-Bmsy_yr1+1);

  // Specify recruitment
  RecPass = 0;
  for (iy=YrPass-5;iy<=YrPass-1;iy++) RecPass += recEst(iy);
  RecPass /= 5;

  Fmsy = M;

  Fmult = Fmsy;
  Proj1YrEst();
  if (mmbEst(YrPass+1) < BmsyEst)
   {
    Fmult = 0;
    Proj1YrEst();
    if (mmbEst(YrPass+1) > beta*BmsyEst)
     {
      Fmult = M/2;
      Proj1YrEst();
      for (ii=1;ii<=10;ii++)
       {
        Fmult = Fmsy*(mmbEst(YrPass+1)/BmsyEst-alpha)/(1-alpha);
        Proj1YrEst();
       }
     }
   }
  Proj1YrEst();
  if (PrintAll==1) cout << mcmc_count << " " << YrPass << " " << BmsyEst << " " << mmbEst(YrPass+1) << " " << Fmult << " ";

// ===========================================================================

FUNCTION Proj1YrEst
  int Yr,Len,Len1;                                // Counter
  dvariable Total;                                // Totals
  dvar_vector encf(1,na),eocf(1,na);              // Proportions for SUMMER fishery
  dvar_vector en0f(1,na),eo0f(1,na);              // Proportions for SUBSISTENCE fishery
  dvar_vector en1f(1,na),eo1f(1,na);              // Proportions for WINTER fishery
  dvar_vector tt0(1,na);                          // Abundance after summer fishery
  dvariable TotalCatch, CatNum;                   // Total catch (in mass)
  dvariable pp,tt1;
  dvariable yss;                                  // Time to summer fishery

  // Selection is for last yr as is time to fishery (as is the timing of the fishery)
  yss = ys(ny);

  // Catch in weight
  TotalCatch = 0;

  // proportions for summer fishery (and the catch)
  Total = 0.0;
  for (Len=1;Len<=na;Len++)
   Total += (npsEst(Len,YrPass)+opsEst(Len,YrPass))*sel(Len,ny)*lg(Len);
  if (Total < 0.01) Total = 0.001;
  CatSum = 0;
  for (Len=1;Len<=na;Len++)
   {
    encf(Len) = npsEst(Len,YrPass)*sel(Len,ny)*lg(Len)/Total;
    CatNum = npsEst(Len,YrPass)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
    CatSum += CatNum;
    TotalCatch += CatNum*wm(Len);
    eocf(Len) = opsEst(Len,YrPass)*sel(Len,ny)*lg(Len)/Total;
    CatNum = opsEst(Len,YrPass)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
    CatSum += CatNum;
    TotalCatch += CatNum*wm(Len);
   }

  // Survivors from the summer fishery
  for (Len=1;Len<=na;Len++)
   {
    tt0(Len) = ((npsEst(Len,YrPass)+opsEst(Len,YrPass))*exp(-yss*Mn(Len))-CatSum*(encf(Len)+eocf(Len)))*exp(-(0.583-yss)*Mn(Len));
    if (tt0(Len) < 0.0) tt0(Len)= 0.0;
   }

  // Upodate dynamics (after summer fishery)
  for (Len=1;Len<=na;Len++)
   {
    pp = 0;
    for (Len1=1;Len1<=Len;Len1++) pp += gr(Len1,Len)*tt0(Len1)*mp0(Len1);
    npwEst(Len,YrPass) = pp;
    if (Len==1) npwEst(1,YrPass) += r1*RecPass;
    if (Len==2) npwEst(1,YrPass) += (1.0-r1)*RecPass;
    opwEst(Len,YrPass) = tt0(Len)*(1.0-mp0(Len));
    if (opwEst(Len,YrPass)<0.0) opwEst(Len,YrPass) = 0.0001;
   }

  //propotions for winter fishery and subsistence fishery
  pp = 0.0; tt1 = 0.0;
  for (Len=3;Len<=na;Len++)
   {
    pp += (npwEst(Len,YrPass)+opwEst(Len,YrPass)*sow)*selw(Len);
    tt1 += (npwEst(Len,YrPass)+opwEst(Len,YrPass)*sow)*selw(Len)*lg(Len);
   }
  if (pp <= 0.0) pp = 0.000001;
  if (tt1 <= 0.0) tt1 = 0.000001;

  CatSub = 0; CatWin = 0;
  en0f.initialize();
  eo1f.initialize();
  for (Len=3;Len<=na;Len++)
   {
    en0f(Len) = npwEst(Len,YrPass)*selw(Len)/pp;
    CatNum = npwEst(Len,YrPass)*selw(Len)*Fmult*FratSub;
    CatSub += CatNum;
    TotalCatch += CatNum*wm(Len);
    eo0f(Len) = opwEst(Len,YrPass)*sow*selw(Len)/pp;
    CatNum = opwEst(Len,YrPass)*sow*selw(Len)*Fmult*FratSub;
    CatSub += CatNum;
    TotalCatch += CatNum*wm(Len);
    en1f(Len) = npwEst(Len,YrPass)*selw(Len)*lg(Len)/tt1;
    CatNum = npwEst(Len,YrPass)*selw(Len)*lg(Len)*Fmult*FratWin;
    CatWin += CatNum;
    TotalCatch += CatNum*wm(Len);
    eo1f(Len) = opwEst(Len,YrPass)*sow*selw(Len)*lg(Len)/tt1;
    CatNum = opwEst(Len,YrPass)*sow*selw(Len)*lg(Len)*Fmult*FratWin;
    CatWin += CatNum;
    TotalCatch += CatNum*wm(Len);
   }

  //update to abundance in next summer (new shell)
  for (Len=1;Len<=na;Len++)
   npsEst(Len,YrPass+1) = (npwEst(Len,YrPass)-en0f(Len)*CatSub-en1f(Len)*CatWin)*exp(-0.417*Mn(Len));
  //update to abundance in next summer (old shell)
  for (Len=1;Len<=na;Len++)
   opsEst(Len,YrPass+1) = (opwEst(Len,YrPass)-eo0f(Len)*CatSub-eo1f(Len)*CatWin)*exp(-0.417*Mo(Len));
  //for (Len=1;Len<=na;Len++) cout << npsEst(Len,YrPass+1) << " ";

  // Compute legal and MMB
  mmbEst(YrPass+1) = 0.0;
  for (Len=3; Len<=na; Len++) mmbEst(YrPass+1) += (npsEst(Len,YrPass+1)+opsEst(Len,YrPass+1))*wm(Len);

  PassCatch = TotalCatch;

// ==========================================================================

FUNCTION ProjAll
  int FutYr,Len,Len1,iy,II;                       // Counters
  dvariable yss;                                  // Time to summer fishery
  dvariable Total;                                // Totals
  dvar_vector encf(1,na),eocf(1,na);              // Proportions for SUMMER fishery
  dvar_vector en0f(1,na),eo0f(1,na);              // Proportions for SUBSISTENCE fishery
  dvar_vector en1f(1,na),eo1f(1,na);              // Proportions for WINTER fishery
  dvar_vector tt0(1,na);                          // Abundance after summer fishery
  dvariable pp,tt1,RecGen;                        // Various variable
  dvariable Term1,Term2;                          // Need for the SR calculations
  dvariable TrueOFL,Catch,FutCatch;               // Catches
  dvariable PredCatchW,CatNum,Fmin,Fmax,FutMort;  // Solve for F

  // Selection is for last yr as is time to fishery (as is the timing of the fishery)
  yss = ys(ny);

  for (FutYr=ny;FutYr<=ny+Nproj-1;FutYr++)
   {

    // Generate this year's recruitment (step 4.c.vi)
    if (SR_rel == 2)
     {
      Term1 = 4*R0*Steep*mmb(FutYr-Lag)/MMB0;
      Term2 = (1-Steep) + (5*Steep-1)*mmb(FutYr-Lag)/MMB0;
      rec(FutYr) = Term1/Term2*EstRec(mcmc_count,FutYr);
     }
    else
     {
      RecGen = rec(RecInd(mcmc_count,FutYr));
      rec(FutYr) = RecGen;
     }

    // Store stuff
    npsSto = nps; opsSto = ops;
    npwSto = npw; opwSto = opw;
    mmbSto = mmb; recSto = rec;

    // Compute the true OFL (step 4.c.i)
    YrPass = FutYr;
    Find_OFL();
    TrueOFL = PassCatch;

    // Store stuff
    nps = npsSto; ops = opsSto;
    npw = npwSto; opw = opwSto;
    mmb = mmbSto; rec = recSto;

    // Estimated values (step 4.c.ii.2)
    for (iy=1;iy<=FutYr;iy++)
     for (Len=1;Len<=na;Len++)
      { npsEst(Len,iy) = EstInd(mcmc_count,FutYr)*nps(Len,iy);
        opsEst(Len,iy) = EstInd(mcmc_count,FutYr)*ops(Len,iy); }
    mmbEst = EstInd(mcmc_count,FutYr)*mmb;
    recEst = EstInd(mcmc_count,FutYr)*rec;

    // Compute the OFL (step 4.c.iii) and then adjust by the buffer (step 4.c.iv)
    YrPass = FutYr;
    Find_OFL_EST();
    Catch = PassCatch;
    FutCatch = PassCatch*Buffer;

    Fmin = 0; Fmax = 40;
    for (II=1;II<=40;II++)
     {
      FutMort = (Fmin+Fmax)/2.0;

      // Reset predicted catch
      PredCatchW = 0;

      // proportions for summer fishery (and the catch)
      Total = 0.0;
      for (Len=1;Len<=na;Len++)
       Total += (nps(Len,FutYr)+ops(Len,FutYr))*sel(Len,ny)*lg(Len);
      if (Total < 0.01) Total = 0.001;
      FutCatSum = 0;
      for (Len=1;Len<=na;Len++)
       {
        encf(Len) = nps(Len,FutYr)*sel(Len,ny)*lg(Len)/Total;
        CatNum = nps(Len,FutYr)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*FutMort;
        FutCatSum += CatNum;
        PredCatchW += CatNum*wm(Len);
        eocf(Len) = ops(Len,FutYr)*sel(Len,ny)*lg(Len)/Total;
        CatNum = ops(Len,FutYr)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*FutMort;
        FutCatSum += CatNum;
        PredCatchW += CatNum*wm(Len);
       }

      // Survivors from the summer fishery
      for (Len=1;Len<=na;Len++)
       {
        tt0(Len) = ((nps(Len,FutYr)+ops(Len,FutYr))*exp(-yss*Mn(Len))-FutCatSum*(encf(Len)+eocf(Len)))*exp(-(0.583-yss)*Mn(Len));
        if (tt0(Len) < 0.0) tt0(Len)= 0.0;
       }

      // Update dynamics (after summer fishery)
      for (Len=1;Len<=na;Len++)
       {
        pp = 0;
        for (Len1=1;Len1<=Len;Len1++) pp += gr(Len1,Len)*tt0(Len1)*mp0(Len1);
        npw(Len,FutYr) = pp;
        if (Len==1) npw(1,FutYr) += r1*rec(FutYr);
        if (Len==2) npw(1,FutYr) += (1.0-r1)*rec(FutYr);
        opw(Len,FutYr) = tt0(Len)*(1.0-mp0(Len));
        if (opw(Len,FutYr)<0.0) opw(Len,FutYr) = 0.0001;
       }

      //propotions for winter fishery and subsistence fishery
      pp = 0.0; tt1 = 0.0;
      for (Len=3;Len<=na;Len++)
       {
        pp += (npw(Len,FutYr)+opw(Len,FutYr)*sow)*selw(Len);
        tt1 += (npw(Len,FutYr)+opw(Len,FutYr)*sow)*selw(Len)*lg(Len);
       }
      if (pp <= 0.0) pp = 0.000001;
      if (tt1 <= 0.0) tt1 = 0.000001;

      FutCatSub = 0; FutCatWin = 0;
      en0f.initialize();
      eo1f.initialize();
      for (Len=3;Len<=na;Len++)
       {
        en0f(Len) = npw(Len,FutYr)*selw(Len)/pp;
        CatNum = npw(Len,FutYr)*selw(Len)*FutMort*FratSub;
        FutCatSub += CatNum;
        PredCatchW += CatNum*wm(Len);
        eo0f(Len) = opw(Len,FutYr)*sow*selw(Len)/pp;
        CatNum = opw(Len,FutYr)*sow*selw(Len)*FutMort*FratSub;
        FutCatSub +=  CatNum;
        PredCatchW += CatNum*wm(Len);
        en1f(Len) = npw(Len,FutYr)*selw(Len)*lg(Len)/tt1;
        CatNum = npw(Len,FutYr)*selw(Len)*lg(Len)*FutMort*FratWin;
        FutCatWin += CatNum;
        PredCatchW += CatNum*wm(Len);
        eo1f(Len) = opw(Len,FutYr)*sow*selw(Len)*lg(Len)/tt1;
        CatNum = opw(Len,FutYr)*sow*selw(Len)*lg(Len)*FutMort*FratWin;
        FutCatWin += CatNum;
        PredCatchW += CatNum*wm(Len);
       }

      if (PredCatchW > FutCatch)
       Fmax = FutMort;
      else
       Fmin = FutMort;

     }

    //update to abundance in next summer (new shell)
    for (Len=1;Len<=na;Len++)
     nps(Len,FutYr+1) = (npw(Len,FutYr)-en0f(Len)*FutCatSub-en1f(Len)*FutCatWin)*exp(-0.417*Mn(Len));
    //update to abundance in next summer (old shell)
    for (Len=1;Len<=na;Len++)
     ops(Len,FutYr+1) = (opw(Len,FutYr)-eo0f(Len)*FutCatSub-eo1f(Len)*FutCatWin)*exp(-0.417*Mo(Len));

    // Compute legal and MMB
    mmb(FutYr+1) = 0.0;
    for (Len=3; Len<=na; Len++) mmb(FutYr+1) += (nps(Len,FutYr+1)+ops(Len,FutYr+1))*wm(Len);
    if (PrintAll == 1) cout << mmb(FutYr+1) << " " << rec(FutYr) << " " << " " << FutCatSum+FutCatWin+FutCatSub << endl;

    mcmc2 << " " << FirstYr+FutYr-1 << " " << mmb(FutYr) << " " << mmb(FutYr)/MMB0 << " " << TrueOFL << " " << Catch << " " << FutCatch << " " << rec(FutYr) << " " << FutMort << endl;

   }

// ==========================================================================
//  Ignore this section 
// ==========================================================================
FUNCTION ProjConstF
  int FutYr,Len,Len1,iy;                          // Counters
  dvariable yss;                                  // Time to summer fishery
  dvariable Total;                                // Totals
  dvar_vector encf(1,na),eocf(1,na);              // Proportions for SUMMER fishery
  dvar_vector en0f(1,na),eo0f(1,na);              // Proportions for SUBSISTENCE fishery
  dvar_vector en1f(1,na),eo1f(1,na);              // Proportions for WINTER fishery
  dvar_vector tt0(1,na);                          // Abundance after summer fishery
  dvariable pp,tt1,RecGen;                        // Various variable
  dvariable TotalCatch,NumYears,CatNum;           // Used to compute average catch
  dvariable FutRec;                               // Future recruitment
  dvariable Term1,Term2;

  // Selection is for last yr as is time to fishery (as is the timing of the fishery)
  yss = ys(ny);

  // Specify recruitment
  FutRec = 0;
  for (iy=ny-5;iy<=ny-1;iy++) FutRec += rec(iy);
  FutRec /= 5;

  TotalCatch = 0; NumYears = 0;
  for (FutYr=ny;FutYr<=ny+Nproj-1;FutYr++)
   {
    if (SR_rel == 2)
     {
      Term1 = 4*R0*Steep*mmb(FutYr-Lag)/MMB0;
      Term2 = (1-Steep) + (5*Steep-1)*mmb(FutYr-Lag)/MMB0;
      rec(FutYr) = Term1/Term2;
     }
    else
     rec(FutYr) = FutRec;

    // Catch in weight
    TotalCatch = 0;

    // proportions for summer fishery (and the catch)
    Total = 0.0;
    for (Len=1;Len<=na;Len++)
     Total += (nps(Len,FutYr)+ops(Len,FutYr))*sel(Len,ny)*lg(Len);
    if (Total < 0.01) Total = 0.001;
    CatSum = 0;
    for (Len=1;Len<=na;Len++)
     {
      encf(Len) = nps(Len,FutYr)*sel(Len,ny)*lg(Len)/Total;
      CatNum = nps(Len,FutYr)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
      CatSum += CatNum;
      TotalCatch += CatNum*wm(Len);
      eocf(Len) = ops(Len,FutYr)*sel(Len,ny)*lg(Len)/Total;
      CatNum = ops(Len,FutYr)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
      CatSum += CatNum;
      TotalCatch += CatNum*wm(Len);
     }

    // Survivors from the summer fishery
    for (Len=1;Len<=na;Len++)
     {
      tt0(Len) = ((nps(Len,FutYr)+ops(Len,FutYr))*exp(-yss*Mn(Len))-CatSum*(encf(Len)+eocf(Len)))*exp(-(0.583-yss)*Mn(Len));
      if (tt0(Len) < 0.0) tt0(Len)= 0.0;
     }

    // Update dynamics (after summer fishery)
    for (Len=1;Len<=na;Len++)
     {
      pp = 0;
      for (Len1=1;Len1<=Len;Len1++) pp += gr(Len1,Len)*tt0(Len1)*mp0(Len1);
      npw(Len,FutYr) = pp;
      if (Len==1) npw(1,FutYr) += r1*rec(FutYr);
      if (Len==2) npw(1,FutYr) += (1.0-r1)*rec(FutYr);
      opw(Len,FutYr) = tt0(Len)*(1.0-mp0(Len));
      if (opw(Len,FutYr)<0.0) opw(Len,FutYr) = 0.0001;
     }

    //propotions for winter fishery and subsistence fishery
    pp = 0.0; tt1 = 0.0;
    for (Len=3;Len<=na;Len++)
     {
      pp += (npw(Len,FutYr)+opw(Len,FutYr)*sow)*selw(Len);
      tt1 += (npw(Len,FutYr)+opw(Len,FutYr)*sow)*selw(Len)*lg(Len);
     }
    if (pp <= 0.0) pp = 0.000001;
    if (tt1 <= 0.0) tt1 = 0.000001;

    CatSub = 0; CatWin = 0;
    en0f.initialize();
    eo1f.initialize();
    for (Len=3;Len<=na;Len++)
     {
      en0f(Len) = npw(Len,FutYr)*selw(Len)/pp;
      CatNum = npw(Len,FutYr)*selw(Len)*Fmult*FratSub;
      CatSub += CatNum;
      TotalCatch += CatNum*wm(Len);
      eo0f(Len) = opw(Len,FutYr)*sow*selw(Len)/pp;
      CatNum = opw(Len,FutYr)*sow*selw(Len)*Fmult*FratSub;
      CatSub +=  CatNum;
      TotalCatch += CatNum*wm(Len);
      en1f(Len) = npw(Len,FutYr)*selw(Len)*lg(Len)/tt1;
      CatNum = npw(Len,FutYr)*selw(Len)*lg(Len)*Fmult*FratWin;
      CatWin += CatNum;
      TotalCatch += CatNum*wm(Len);
      eo1f(Len) = opw(Len,FutYr)*sow*selw(Len)*lg(Len)/tt1;
      CatNum = opw(Len,FutYr)*sow*selw(Len)*lg(Len)*Fmult*FratWin;
      CatWin += CatNum;
      TotalCatch += CatNum*wm(Len);
     }

    //update to abundance in next summer (new shell)
    for (Len=1;Len<=na;Len++)
     nps(Len,FutYr+1) = (npw(Len,FutYr)-en0f(Len)*CatSub-en1f(Len)*CatWin)*exp(-0.417*Mn(Len));
    //update to abundance in next summer (old shell)
    for (Len=1;Len<=na;Len++)
     ops(Len,FutYr+1) = (opw(Len,FutYr)-eo0f(Len)*CatSub-eo1f(Len)*CatWin)*exp(-0.417*Mo(Len));

    // Compute legal and MMB
    mmb(FutYr+1) = 0.0;
    for (Len=3; Len<=na; Len++) mmb(FutYr+1) += (nps(Len,FutYr+1)+ops(Len,FutYr+1))*wm(Len);
    if (PrintDiag == 1) cout << FutYr << " " << Fmult << " " << mmb(FutYr+1) << " " << rec(FutYr) << " " << " " << CatSum+CatWin+CatSub << " " << mmb(FutYr+1)/rec(FutYr) << endl;
    PassCatch = TotalCatch;
   }

// ==========================================================================

FUNCTION get_reference_points
  cout<<"start reference points"<<endl;
  int i,j,k,kk;
  dvar_matrix ref_catch_s(1,na,1,100);
  dvar_matrix ref_catch_sd(1,na,1,100);
  dvar_matrix ref_catch_w_c(1,na,1,100);
  dvar_matrix ref_catch_w_s(1,na,1,100);
  dvar_matrix ref_catch_w_cd(1,na,1,100);
  dvar_matrix ref_catch_w_sd(1,na,1,100); 
  dvar_matrix ref_catch(1,na,1,100);
  dvar_matrix ref_ps(1,na,1,100);
  dvar_matrix ref_nps(1,na,1,100);
  dvar_matrix ref_ops(1,na,1,100);
  dvar_matrix ref_npw(1,na,1,100);
  dvar_matrix ref_opw(1,na,1,100);
  dvar_vector ref_catch_b(1,100);
  dvar_vector ref_catch_n(1,100);
  dvariable pp, wpm, wpim, tt1, ttt0, ttt1, sb, tw;
  dvar_vector nscaf(1,na), tsc(1,na),tscd(1,na), tt0(1,na);    
  dvar_vector sselw(1,na),dselw(1,na),lselw(1,na),dlselw(1,na);
  dvar_vector sselc(1,na),dselc(1,na); 
  dvar_vector nwc(1,na),owc(1,na),nwcd(1,na),owcd(1,na);
  dvar_vector nws(1,na),ows(1,na),nwsd(1,na),owsd(1,na); 
  dvar_vector par1(1,na),par2(1,na),fx(1,na), Fs(1,na),Fw(1,na);   
  dvariable TotalCatch, CatNum, cofl,fl;    
  dvar_matrix ref_ws(1,na,1,4);  
  dvariable ref_F;
  dvariable ref_Hwc;
  dvariable ref_Hts;
  dvar_vector ref_Htw(1,ny);
  dvar_vector ref_mbio(1,101);
  dvar_vector ref_mmb(1,100);
  dvar_vector ref_totc(1,101);
  dvar_vector ref_legaln(1,101);
  dvar_vector ref_Hr(1,101);
  dvar_vector ref_legal(1,100);
  dvar_vector ref_legal_b(1,100);
  dvariable ref_ys;
  dvariable f35,f40,h35,h40,i35,i40;
  dvariable tt0,pp;
  
  ref_Htw.initialize();
  
  ref_ys = ys(ny-1);
  ref_Hwc = 0.0;
//--------------------------------------------------------------------------
// Determine average subsistence harvest rate (10 years average)  
//--------------------------------------------------------------------------
  for (i = ny-10; i<= ny-1; i++)
   {
    ref_FSub += FSub(i);
   }
  ref_FSub = ref_FSub/10.0;     
   yss = ys(lyear);
// Calculate winter subsistence harvest selectivity
// Assume subsistence harvest take all mature crab caught
   sselw = elem_prod(matc,selw); 
// Assume subsistence harvest discards all immature crab caught   
   dselw = elem_prod(imatc,selw); 
// Calculate legal crab caught by winter gear    
   lselw = elem_prod(lg,selw);   
// Calculate sublegal crab caught by winter gear (discards)  
   dlselw = elem_prod(slg,selw);      
// Calculate legal crab caught by Winter gear       
   sselc = elem_prod(sel(lyear),lg);   
// Calculate sublegal crab caught by Winter gear          
   dselc = elem_prod(sel(lyear),slg);  

   
//--------------------------------------------------------------------------
// Simulation for 100 years   
//--------------------------------------------------------------------------
  ref_Hts = 0.0;
  for (j = 1; j<= 101; j++)
   {
  ref_F = 0.01*j-0.01;
//--------------------------------------------------------------------------
// 8.1 Calculate Winter and Summer commercial Harvest rate: Fw, Fs  
//--------------------------------------------------------------------------	    	       
  par1 = 1-exp(-ref_F-0.42*Mn);
  par2 = 1-pwh*(1-exp(-0.42*Mn));
  fx = elem_div(1-pwh*par1,par2); 
  Fw = 1-fx;
  Fs = elem_prod(fx-exp(-ref_F),exp(-0.42*Mn)); 
// Initialize temp matrix
  ref_nps.initialize();
  ref_ops.initialize();
  ref_npw.initialize();
  ref_opw.initialize();
// Set First year
  ref_npw(1) = npw(lyear);
  ref_opw(1) = opw(lyear);
  ref_pw(1)= ref_npw(1)+ref_opw(1); 
	
//==========================================================================
// 8.7  Start 100 year Simulation 
//==========================================================================
    for (i=1;i<=100;i++)
     {
//--------------------------------------------------------------------------
// 8.1 Winter commercial catch and discards size composition 
//--------------------------------------------------------------------------	    	       
// Calculate number of newshell and oldshell in commercial catch 	  
	nwc = elem_prod(elem_prod(ref_npw(i),lselw),Fw);
	owc = elem_prod(elem_prod(ref_opw(i),lselw),Fw)
// Calculate number of newshell and oldshell discards in commercial  	  
	nwcd = elem_prod(elem_prod(ref_npw(i),dlselw),Fw);
	owcd = elem_prod(elem_prod(ref_opw(i),dlselw),Fw);
    ref_catch_w_c(i) = nwc + owc;
    ref_catch_w_cd(i) = nwcd + owcd;
	
//--------------------------------------------------------------------------
// 8.2 Winter subsistence catch and discards size composition 
//     Assume that Subsistence fishers took all length >= 3 crab  
//	   and discarded length 1 & 2 crab
//--------------------------------------------------------------------------
// Calculate number of newshell and oldshell in subsistence catch 
		nws = elem_prod(ref_npw(i),sselw)*FSub;
		ows = elem_prod(ref_opw(i),sselw)*FSub;
// Calculate number of newshell and oldshell in subsistence discards  		
		nwsd = elem_prod(ref_npw(i),dselw)*FSubd;
        owsd = elem_prod(ref_opw(i),dselw)*FSubd;
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
	nsc = elem_prod(elem_prod(ref_nps(i),sselc),Fs);
	osc = elem_prod(elem_prod(ref_ops(i),sselc),Fs); 
// Calculate number of newshell and oldshell discards in commercial  	  
	nscd = elem_prod(elem_prod(ref_nps(YrPass),dlselc),Fs);
	oscd = elem_prod(elem_prod(ref_ops(YrPass),dlselc),Fs);
// Total summer commercial catch 	  
	ref_catch_s(i) = nsc+osc;
	ref_catch_sd(i) = nscd+oscd;
//--------------------------------------------------------------------------
// 8.6  Calculate number of Summer Crab after summer fishery (nscaf)
//      nscaf = July 1st abundance*mortality till fishery 
//               - com catch and discards. 
//--------------------------------------------------------------------------
  nscaf = elem_prod(ref_nps(i)+ref_ops(i),exp(-ys*Mn))-ref_catch_s(i)-ref_catch_sd(i)*hm(1);	  
	
// Stop gap measure: abundance does not go below zero
    for (k=1;k<=na;k++)
     {	  
      if (nscaf(k) < 0.0) nscaf(k)= 0.001;  
	  }
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
	  //New shell crab molted + recruits
      ref_npw(i+1,j) = pp*exp(-(0.583-ref_ys)*Mn(j)) + rpp(j)*rec(lyear);    
	  //Old shell crab are unmolted 
      ref_opw(i+1,j) = nscaf(j)*(1.0-mp1(j))*exp(-(0.583-ref_ys)*Mn(j));    
     }
	
//==========================================================================
// 8.7  Calculate February 1st legal Crab and mature male biomass  for assessment
//==========================================================================

  ref_mmb(i) = sum(elem_prod(elem_prod(ref_npw(i)+ref_opw(i),matc),wm));  // Mature male biomass
  ref_legal(i) = sum(elem_prod(ref_npw(i)+ref_opw(i),lg));  // number of legal crab
  ref_legalb(i) = sum(elem_prod(ref_npw(i)+ref_opw(i),wmlg)); // legal biomass   
	
    ref_mbio(j) = ref_mmb(100);                    //lbs/R
//    ref_totc(j) = ref_catch_b(100);
    ref_legaln(j) = ref_legal(100);
    ref_Hr(j) = ref_catch_n(100)/ref_legaln(j);
   }
   }
  i35 = 0; i40 = 0;
  for (j = 1; j<= 101; j++)
   {
    if (i35 < 1.0)
// Find F35 F35 = 35% of reference mmb		
     if (ref_mbio(j) <= 0.35*ref_mbio(1))
      {
       f35 = 0.01*j-0.01;
       h35 = ref_Hr(j);
       i35 = 2.0;
      }
    if (i40 < 1.0)
     if (ref_mbio(j) <= 0.40*ref_mbio(1))
      {
       f40 = 0.01*j-0.01;
       h40 = ref_Hr(j);
       i40 = 2.0;
      }
   }

  ofstream report1("refp.out");
  report1 <<"Harvest rate in term of legal males as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_Hr<<endl;
  report1 <<"Total MMB (>93mm) Feb 01 as F = 0.00, 0.01, ... 1.0 (lbs/R)"<<endl;
  report1 << ref_mbio<<endl;
  report1 <<"Total catch as F = 0.00, 0.01, ... 1.0 (lbs/R)"<<endl;
  report1 << ref_totc<<endl;
  report1 <<"Total legals as F = 0.00, 0.01, ... 1.0 (crabs/R)"<<endl;
  report1 << ref_legaln<<endl;
  report1 <<"F35: "<<f35<<"  H35%: "<<h35<<endl;
  report1 <<"F40: "<<f40<<"  H40%: "<<h40<<endl;
  report1 <<"ref_Hwc = (mean of 10 years):  "<<ref_Hwc<<endl;
  cout<<"end of reference points"<<endl;

// ========================================================================
// ========================================================================

REPORT_SECTION
  cout << "Report Section" << endl;
  get_reference_points();
  