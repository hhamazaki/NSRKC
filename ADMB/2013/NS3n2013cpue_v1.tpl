//NS3n2013cpue.TPL
//This model is a NS3n20127c.TPL model adopted by the CPT in 2012
//The model incorporated cpue instead of effort


DATA_SECTION
  init_int ny                 //number of years
  init_int na                 //number of length groups
  init_int nyt                //number of years with trawl survey
  init_int nyp                //number of years with pot survey
  init_int nyw                //number of years with winter project
  init_int nyf                //number of years with pre-fishery survey
  init_int nyo                //number of years with observer's data
  init_int scy1                //year large summer commercial trawl fishery ended
  init_int scy2                //year summer commercial fishery escapement mechanism became requirement
  init_int scp                //year summer commrecial fishery CW>5 inch crab was accepted by buyers
  init_number slm             //mid-size of smallest length group (mm)
  init_number slt             //length interval (mm)
  init_number M               //instantaneous natural mortality
  init_number ms			  //instantaneous natural mortality multiplier for the last size class
  init_number slt6			  // trawl catchability multiplier for the last size class
  init_number qt			  // trawl catchability multiplier
  init_number spcv			  // summer pot survey CV
  init_number lamt            // NO LONGER USED weighting factor for trawl survey
  init_number lamp            // NO LONGER USED weighting factor for pot survey
  init_number lamc            //weighting factor for fis hing effort: summer fishery
  init_number maxss           //maximum effective sample size for length proportion
  init_number efn             //% of effective sample size due to multin.distribution
  init_vector lg(1,na)        //proportion of legals by length group
  init_matrix gr(1,na,1,na)   //growth matrix
  init_ivector it(1,nyt)      //year id for trawl survey data
  init_vector tt(1,nyt)       //total annual abundance from trawl survey
  init_vector pct(1,nyt)      //% of summer catch occurred before the mid point of survey
  init_vector yt(1,nyt)       //Mid point of trawl survey from July 1
  init_vector st(1,nyt)       //annual sample size from trawl survey
  init_vector cv1(1,nyt)       //annual cv of sublegals for trawl survey
  init_vector cv2(1,nyt)       //annual cv of legals for trawl survey
  init_vector cv3(1,nyt)       //annual cv of totals for trawl survey
  init_ivector ip(1,nyp)      //year id for pot survey data
  init_vector tp(1,nyp)       //total annual abundance from pot survey
  init_vector yp(1,nyp)       //Mid point of pot survey from July 1
  init_vector sp(1,nyp)       //annual sample size from pot survey
  init_ivector iw(1,nyw)      //year id for winter project
  init_vector sw(1,nyw)       //annual sample size from winter project
  init_int ipre               //year id for pre-fishery survey
  init_number spre            //annual sample size from pre-fishery survey
  init_vector sc(1,ny)        //annual sample size from summer fishery
  init_vector tc(1,ny)        //annual catch for summer fishery
  init_vector te(1,ny)        //annual effort for summer fishery
  init_vector stcpue(1,ny)    //annual standardized cpue for summer fishery
  init_vector secpue(1,ny)    //se of  annual cpue for summer fishery
  init_vector ys(1,ny)        //Mid point of summer fishery from July 1
  init_ivector io(1,nyo)      //year id for observer's data
  init_vector so(1,nyo)       //annual sample size from observer's data
  init_vector twc(1,ny-1)     //annual catch for winter fishery
  init_vector tsc(1,ny-1)     //annual catch for winter subsistence fishery
  init_matrix ont(1,na,1,nyt) //proportion of length group for new-shell, trawl survey
  init_matrix oot(1,na,1,nyt) //proportion of length group for old-shell, trawl survey
  init_matrix onp(1,na,1,nyp) //proportion of length group for new-shell, pot survey
  init_matrix oop(1,na,1,nyp) //proportion of length group for old-shell, pot survey
  init_matrix onw(1,na,1,nyw) //proportion of length group for new-shell, winter survey
  init_matrix oow(1,na,1,nyw) //proportion of length group for old-shell, winter survey
  init_vector onf(1,na)       //proportion of length group for new-shell, pre-fishery survey
  init_vector oof(1,na)       //proportion of length group for old-shell, pre-fishery survey
  init_matrix onc(1,na,1,ny)  //proportion of length group for new-shell, summer fishery
  init_matrix ooc(1,na,1,ny)  //proportion of length group for old-shell, summer fishery
  init_matrix ono(1,na,1,nyo) //proportion of length group for new-shell, observer's data
  init_matrix ooo(1,na,1,nyo) //proportion of length group for old-shell, observer's data
  init_vector wm(1,na)        //Mean weight
  init_number hm              //handling mortality rate for bycatch
  vector tt1(1,nyt)            //observed annual non-legal crab abundance: trawl survey
  vector tt2(1,nyt)            //observed annual legal crab abundance: trawl survey
  vector tp1(1,nyp)            //observed annual non-legal crab abundance: pot survey
  vector tp2(1,nyp)            //observed annual legal crab abundance: pot survey

  !! ad_comm::change_datafile_name("proj.ctl");
  init_int ProjType            // Set to 1 for Dana and 2 for Five-year forecasts
  init_int FirstYr             // First year for assessment
  int Nproj;                   // Projection length
  init_int NprojFmsy           // Projection length (for computing FMSY)
  init_int NprojSim            // Projection length (for actual simulations)
  !! Nproj = NprojFmsy;
  !! if (NprojSim > Nproj) Nproj = NprojSim;
  init_number alpha;           // Alpha - control rule
  init_number beta;            // Beta - control rule
  init_int Bmsy_yr1;           // Years for defining BMSY
  init_int Bmsy_yr2;           // Years for defining BMSY
  init_number CVEstWithin;     // CV of biomass estimate (captured by the assessment)
  init_number CVEstExtra;      // CV of biomass estimate (not captured by the assessment)
  init_number Buffer           // Buffer between OFL and ABC
  init_int BurnIn;             // Burn-in period
  init_int LastSim;            // Number of sims
  !! LastSim = LastSim + BurnIn;
  init_int SR_FORM;            // Stock-recruitment relationship
  init_number Steep;           // Steepness
  init_number SigmaR;          // Recruitment variation
  init_int Lag;                // Lag to recruitment
  init_int PrintAll;           // Print all results

  int mcmc_count;              // Counter for MCMC
  !! mcmc_count = 0;

  int fn_count;
  !! fn_count = 0;

  int MaxSim;
  !! MaxSim = 1000;

  int YrPass;                                    // Year pass (to calculate TAC)
  imatrix RecInd(1,MaxSim,ny,ny+Nproj);          // Index for selecting recruitment

  int SR_rel;                                    // Stock-recruitment relationship
  int PrintDiag;                                 // Test outputs

  !! cout << "Data Section Completed" << endl;
  !! cout << wm << endl;
  !! cout << hm << endl;

//==========================================================================

PARAMETER_SECTION
  init_bounded_number log_q1(-12.5,-8.5,1)
  init_bounded_number log_q2(-12.5,-10.0,1)
  init_bounded_number log_q3(-12.5,-10.0,1)
//  init_bounded_number log_initpop(2,15,1)
  init_bounded_vector log_initpop(1,na,0.1,15,1)
  init_bounded_number log_recscale(2,12,1)
  init_bounded_dev_vector log_relrec(1,ny-1,-15.,15.,1)
  init_bounded_number r1(0.5,0.9,1)             //proportion of recruits to length group 1
  init_bounded_number log_mo1(-5.5,-2.0,1)     //parameter for molting probability
  init_bounded_number log_mo2(0.55,5.00,1)       //length at 50% molting probability
  init_bounded_number log_st1(-5.5,-1.0,1)      //parameter for selectivity: NOAA trawl survey
  init_bounded_number log_st2(0.51,3.00,1)      //length at 50% selectivity: NOAA trawl survey
  init_bounded_number log_st3(-5.5,-1.0,1)      //parameter for selectivity: ADF&G trawl survey
  init_bounded_number log_st4(0.51,3.00,1)      //length at 50% selectivity: ADF&G trawl survey
  init_bounded_number log_sp1(-3.0,-1.00,1)      //parameter for selectivity: pot survey
  init_bounded_number log_sp2(3.9,5.50,1)      //length at 50% selectivity: pot survey
  init_bounded_number log_sw1(-3.0,1.00,1)      //parameter for selectivity: winter project
  init_bounded_number log_sw2(3.90,5.50,1)      //length at 50% selectivity: winter project
  init_bounded_number sw3(0.1,1.0,1)           //selectivity for length group NA: winter project
  init_bounded_number log_sc1(-3.0,-1.00,1)      //parameter for selectivity: summer fishery before 1993
  init_bounded_number log_sc2(3.9,5.5,1)      //length at 50% selectivity: summer fishery before 1993
  init_bounded_number log_sc3(-3.0,-1.00,1)      //parameter for selectivity: summer fishery after 1993
  init_bounded_number log_sc4(3.9,5.5,1)      //length at 50% selectivity: summer fishery after 1993
  init_bounded_number log_sc5(-3.0,-1.00,1)      //parameter for selectivity: summer fishery after 2008
  init_bounded_number log_sc6(3.9,5.50,1)      //length at 50% selectivity: summer fishery after 2008
  init_bounded_number sos0(0.50,1.50,-1)         //selectivity of old-shell crabs in summer commercial fishery
  init_bounded_number sow(0.300,1.70,-1)         //selectivity of old-shell crabs in winter project
 // init_bounded_number M(0.18,0.18,-1)         //Mortality
 // init_bounded_number ms(1.6,1.6,-1)         //Mortality multiplier for the last length group
 // init_bounded_number slt6(0.1,1.0,1)         //trawl selectivity of the last length group
  init_bounded_number p4(0.1,1,1)         //proporion of marketable legal for length 4 class since 2005
  vector log_rec(1,ny-1)
  matrix nps(1,na,1,ny+Nproj)                   //new-shell population abundance in summer
  matrix ops(1,na,1,ny+Nproj)                   //old-shell population abundance in summer
  matrix npw(1,na,1,ny+Nproj)                   //new-shell population abundance in winter
  matrix opw(1,na,1,ny+Nproj)                   //old-shell population abundance in winter
  vector rec(1,ny+Nproj)                        //annual recruits: starting in year 1.
  matrix ent(1,na,1,nyt)                        //estimated proportion: new shell,trawl survey
  matrix et0(1,na,1,nyt)                        //estimated trawl survey abundance
//  matrix eot(1,na,1,nyt)                        //estimated proportion: old shell,trawl survey
  matrix enp(1,na,1,nyp)                        //estimated proportion: new shell,pot survey
  matrix eop(1,na,1,nyp)                        //estimated proportion: old shell,pot survey
  matrix enw(1,na,1,nyw)                        //estimated proportion: new shell,winter project
  matrix eow(1,na,1,nyw)                        //estimated proportion: old shell,winter project
  vector enf(1,na)                              //estimated proportion: new shell,pre-fishery survey
  vector eof(1,na)                              //estimated proportion: old shell,pre-fishery survey
  matrix enc(1,na,1,ny)                         //estimated proportion: new shell,summer fishery
  matrix eoc(1,na,1,ny)                         //estimated proportion: old shell,summer fishery
  matrix eno(1,na,1,nyo)                        //estimated proportion: new shell,observer's data
  matrix eoo(1,na,1,nyo)                        //estimated proportion: old shell,observer's data
  matrix en0(1,na,1,ny-1)                       //estimated proportion: new shell,subsistence fishery
  matrix eo0(1,na,1,ny-1)                       //estimated proportion: old shell,subsistence fishery
  matrix en1(1,na,1,ny-1)                       //estimated proportion: new shell,winter fishery
  matrix eo1(1,na,1,ny-1)                       //estimated proportion: old shell,winter fishery
  vector ett(1,nyt)                             //estimated total annual abundance: trawl survey
  vector ett1(1,nyt)                            //estimated annual non-legal crab abundance: trawl survey
  vector ett2(1,nyt)                            //estimated annual legal crab abundance: trawl survey
  vector etp(1,nyp)                             //estimated total annual abundance: pot survey
  vector etp1(1,nyp)                            //estimated annual non-legal crab abundance: pot survey
  vector etp2(1,nyp)                            //estimated annual legal crab abundance: pot survey
  vector ete(1,ny)                              //estimated total annual effort: summer fishery
  vector ecpue(1,ny)                            //estimated cpue: summer fishery
  vector cvcpue(1,ny)                           //cv of cppue
  vector tb(1,ny)                               //estimated mean annual summer exploitable abundance
  vector mp1(1,na)                              //molting probability  1976-1981
  vector mp2(1,na)                              //molting probability  1982-1990
  vector mp3(1,na)                              //molting probability: 1991-1997
  vector mp0(1,na)                              //working vector: molting probability
  vector sel92(1,na)                            //selectivity of summer fishery before 1993
  vector sel93(1,na)                            //selectivity of summer fishery after 1992
  vector sel08(1,na)                            //selectivity of summer fishery after 2008
  matrix sel(1,na,1,ny)							//selectivity of summer fishery combined
  vector selt(1,na)                             //selectivity of trawl survey
  vector selp(1,na)                             //selectivity of pot survey
  vector selw(1,na)                             //working variable: selectivity of winter project.
  vector Mn(1,na)                               //working variable: natural mortality for new-shell crabs.
  vector Mo(1,na)                               //working variable: natural mortality for old-shell crabs.
  vector sos(1,na)                              //working variable: additional selectivity for old-shell.
  number q1                                     //catchability: summer fishery before 1993
  number q2                                     //catchability: summer fishery after 1992
  number q3                                     //catchability: summer fishery after 2008
  number pre                                    //working variable:
  vector bc(1,ny)                               //bycatch biomass.
  vector FSum(1,ny+Nproj)                       // Exploitation rate (Summer fishery)
  vector FWin(1,ny+Nproj)                       // Exploitation rate (Winter fishery)
  vector FSub(1,ny+Nproj)                       // Exploitation rate (Subsistance fishery)
  number FratSub;                               // Fratio for Subsistence fishery
  number FratWin;                               // Fratio for Winter fishery
  number Fmult;                                 // Multipler for F (for decision theory)
  number CatSum                                 // Catches (Summer)
  number CatSub                                 // Catches (Subsistence)
  number CatWin                                 // Catches (Winter)
  number FutCatSum;                             // Future catches (Summer)
  number FutCatSub;                             // Future catches (Subsistence)
  number FutCatWin;                             // Future catches (Winter)
  matrix npsEst(1,na,1,ny+Nproj)                //new-shell population abundance in summer
  matrix opsEst(1,na,1,ny+Nproj)                //old-shell population abundance in summer
  matrix npwEst(1,na,1,ny+Nproj)                //new-shell population abundance in winter
  matrix opwEst(1,na,1,ny+Nproj)                //old-shell population abundance in winter
  vector mmbEst(1,ny+Nproj)                     //annual mmb
  vector recEst(1,ny+Nproj)                     //annual recruits: starting in year 1.
  matrix npsSto(1,na,1,ny+Nproj)                //new-shell population abundance in summer
  matrix opsSto(1,na,1,ny+Nproj)                //old-shell population abundance in summer
  matrix npwSto(1,na,1,ny+Nproj)                //new-shell population abundance in winter
  matrix opwSto(1,na,1,ny+Nproj)                //old-shell population abundance in winter
  vector mmbSto(1,ny+Nproj)                     //annual mmb
  vector recSto(1,ny+Nproj)                     //annual recruits: starting in year 1.
  matrix EstInd(1,MaxSim,ny,ny+Nproj);          // Observation error for estimation
  matrix EstRec(1,MaxSim,ny,ny+Nproj);          // Process error for recruitment
  number MMB0;                                  // virgin MMB0
  number R0;                                    // virgin R0
  number PassCatch;                             // Catch passed
  number RecPass;                               // Passed recruitment
  number BmsyProx;                              // Bmsy proxy
  number fbest
  !! fbest = 1.0e+20;
  number tf1;									// Likelihood
  number tf2;									// Likelihood
  number tf3;									// Likelihood
  number tf4;									// Likelihood
  number tf5;									// Likelihood
  number tf6;									// Likelihood
  number tf7;									// Likelihood
  number tf8;									// Likelihood
  number tf9;									// Likelihood
  number tf10;									// Likelihood
  number tf11;									// Likelihood
  number tf12;									// Likelihood
  number tf13;									// Likelihood
  number tf14;									// Likelihood
  number tf15;									// Likelihood
  objective_function_value f
  sdreport_vector last_y(1,na)
  sdreport_vector legaln(1,ny)
  sdreport_vector legalb(1,ny)
  sdreport_vector mmb0(1,ny)
  vector legal(1,ny+Nproj)
  vector mmb(1,ny+Nproj)
  likeprof_number last_legal
  likeprof_number last_subl
  likeprof_number last_mmb
  likeprof_number last_ofl
  likeprof_number first_year
  number nc;
  number oc;
  vector first_y(1,na)
  !! cout << "Parameter Section Completed" << endl;

//==========================================================================

INITIALIZATION_SECTION
  log_q1         -11.5
  log_q2         -11.5
  log_q3         -11.05
  log_initpop      7.0
  r1               0.6
  log_mo1         -2.41358
  log_mo2          4.8211
  sw3              1.0
  log_sc1         -2.31267
  log_sc2          4.485
  log_sc3         -2.41358
  log_sc4          4.82115
  log_sc5         -1.73
  log_sc6          4.68
  sos0             1.0
  sow              1.00
 // M                0.18
 // ms               1.6
 // slt6             0.6
  p4               0.5

//==========================================================================

PRELIMINARY_CALCS_SECTION
  int i,j,Isim,iy;
  dvariable ErrComm;

  cout << "Starting preliminary calcs" << endl;
  Mo = M;
  Mn = M;
  Mo(na) = ms*M; Mn(na) = ms*M;
  sos = 1.0;
  selw = 1.0;
  selt = 1.0;
  selp = 1.0;
  en0.initialize();
  en1.initialize();
  eo0.initialize();
  eo1.initialize();
  //set the max effective sample size for length composition.
  // The effective sample size are
  // Trawl and summer pots: 50% of sample size or maxx*efn
  // Winter pots, summer commercial, and observers: 10% of sample size or 0.5*maxx*efn

  for (i=1;i<=nyt;i++)
   {
    st(i) *= 0.5;
    if (st(i) > maxss) st(i)=maxss*efn;             //trawl survey
   }
  for (i=1;i<=nyp;i++)
   {
    sp(i) *= 0.5;
    if (sp(i) > maxss) sp(i)=maxss*efn;             //pot survey
   }
  for (i=1;i<=nyw;i++)
   {
    sw(i) *= 0.10;
    if (sw(i) > 0.5*maxss) sw(i)=0.5*maxss*efn;             //winter project
//    if (sw(i) > maxss) sw(i)=maxss*efn;             //winter project
   }
//  pre = 0.0;                                        //do not use pre-fishery data (one year of data only)
//  if (pre > 0.5*maxss) pre=0.5*maxss*efn;                   //pre-fishery
//  if (pre > maxss) pre=maxss*efn;                   //pre-fishery

  for (i=1;i<=ny;i++)
   {
    sc(i) *= 0.10;
    if (sc(i) > 0.5*maxss) sc(i)=0.5*maxss*efn;             //summer fishery
//    if (sc(i) > maxss) sc(i)=maxss*efn;             //summer fishery
   }
  for (i=1;i<=nyo;i++)
   {
    so(i) *= 0.10;
    if (so(i) > 0.5*maxss) so(i)=0.5*maxss*efn;             //observer's data
//    if (so(i) > maxss) so(i)=maxss*efn;             //observer's data
   }
//combine new- and old-shell data for all trawl survey data & compute observed legal & non-legal abundance
  for (i=1;i<=nyt;i++)
   {
    for (j=3;j<=na;j++) ont(j,i) += oot(j,i);
	// tt1 non-legal abundance trawl survey
    tt1(i) = tt(i)*(ont(1,i)+ont(2,i)+ont(3,i));
    // tt2 legal abundance trawl survey
	tt2(i) = tt(i) - tt1(i);
   }

//combine compute observed legal & non-legal abundance
  for (i=1;i<=nyp;i++)
   {
	// tp1 non-legal abundance pot survey
    tp1(i) = tp(i)*(onp(1,i)+onp(2,i)+onp(3,i)+oop(3,i));
	// tp2 legal abundance pot survey
    tp2(i) = tp(i) - tp1(i);
   }


  // Generate ALL the random numbers
  random_number_generator r(9891);
  for (Isim=1;Isim<=MaxSim;Isim++)
   {
    // Recruitment
    for (iy=ny;iy<=ny+Nproj;iy++)
     RecInd(Isim,iy) = int(randu(r)*(ny-1))+1;
    // Process errr
    for (iy=ny;iy<=ny+Nproj;iy++)
     EstRec(Isim,iy) = mfexp(randn(r)*SigmaR-SigmaR*SigmaR/2.0);
    // Biomass estimation (common to all years)
    ErrComm = mfexp(randn(r)*CVEstExtra-square(CVEstExtra)/2.0);
    for (iy=ny;iy<=ny+Nproj;iy++)
     EstInd(Isim,iy) = ErrComm*mfexp(randn(r)*CVEstWithin-square(CVEstWithin)/2.0);
   }
// Calculate cv of cpue
  cvcpue = elem_div(secpue+1.e-3,stcpue+1.e-3);

  mcmc2 << "Year MMB MMN/MMB0 True-OFL Est-OFl ABC/Catch Recr" << endl;
  cout << "End preliminary calcs" << endl;


// ==========================================================================

PROCEDURE_SECTION
  dvariable SBPR;
  int i;

  fn_count += 1;

  convert_parameters_into_rates();
//   cout <<"OK for convert_parameters..."<<endl;

  get_first_year_abundance();
//  cout <<"OK for get_first_year..."<<endl;

  get_number_by_size();
//   cout <<"OK for get_number_by_size..."<<endl;

  get_proportion_and_effort();
//   cout <<"OK for get_proportion_and..."<<endl;

  evaluate_the_objective_function();

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
      else
       {
        SR_rel = 1;
        Fmult = 0;
        PrintDiag = 0;
        YrPass = ny;
        ProjConstF();
        SBPR = mmb(ny+Nproj-1)/rec(ny+Nproj-1);

        // Specify steepness as given above
        Get_Steepness();

        // Real stock-recruitment relationship
        SR_rel = SR_FORM;

        // Now do all those projections
        Nproj = NprojSim;
        if (R0 > 0)
         {
          mcmc2 << Steep << " " << MMB0 << " " << BmsyProx << endl;
          mcmc2 << " " << -1 << " " << MMB0 << " 1 0 0 0 " << R0 << endl;
          for (i=1;i<=ny-1;i++)
           mcmc2 << " " << FirstYr+i-1 << " " << mmb(i) << " " << mmb(i)/MMB0 << " -1 -1 " << tc(i)+twc(i)+tsc(i) << " " << rec(i) << endl;
          ProjAll();
         }
       }
     }
   }

//=============================================================================
//        Calculate molting probability and Pot selectivity
//=============================================================================
FUNCTION convert_parameters_into_rates
  int i, j;
// working variables: molting probability
  dvariable mo1,mo2;
// working variables: Selectivity
  dvariable sc1,sc2,sc3,sc4, sc5, sc6, st1,st2, st3, st4, sp1,sp2,sw1,sw2;
// local variable
  dvariable pp;

//catchability coefficient
   q1 = mfexp(log_q1);
   q2 = mfexp(log_q2);
   q3 = mfexp(log_q3);
//molting probability and selectivity for summer commercial fisheries
   mo1 = mfexp(log_mo1);
   mo2 = mfexp(log_mo2);
   sc1 = mfexp(log_sc1);
   sc2 = mfexp(log_sc2);
   sc3 = mfexp(log_sc3);
   sc4 = mfexp(log_sc4);
   sc5 = mfexp(log_sc5);
   sc6 = mfexp(log_sc6);
   st1 = mfexp(log_st1);
   st2 = mfexp(log_st2);
   st3 = mfexp(log_st3);
   st4 = mfexp(log_st4);
   sp1 = mfexp(log_sp1);
   sp2 = mfexp(log_sp2);
   sw1 = mfexp(log_sw1);
   sw2 = mfexp(log_sw2);

 // calculate molitng and catchability function
   for (j=1;j<=na;j++)
    {
     pp = (double(j)-1.0)*slt;
     pp += slm;
     mp1(j) = 1.0-1.0/(1.0+mfexp(-1.0*mo1*(pp-mo2)));
//     mp2(j) = 1.0-1.0/(1.0+mfexp(-1.0*mo1*(pp-mo2)));
//     mp3(j) = 1.0-1.0/(1.0+mfexp(-1.0*mo1*(pp-mo2)));
     sel92(j) = 1.0/(1.0+mfexp(-1.0*sc1*(pp-sc2)));
     sel93(j) = 1.0/(1.0+mfexp(-1.0*sc3*(pp-sc4)));
//	 sel08(j) = 1.0/(1.0+mfexp(-1.0*sc5*(pp-sc6)));
    }
// Adjust selectivity so that selectivity of 5th size class is 1	
   sel92 = sel92/sel92(na-1);
   if (sel92(na) > 1.0) sel92(na) = 1.0;
   sel93 = sel93/sel93(na-1);
   if (sel93(na) > 1.0) sel93(na) = 1.0;
//  sel08 = sel08/sel08(na-1);
//  if (sel08(na) > 1.0) sel08(na) = 1.0;

   for (i = 1; i<=ny; i++)
    {
     if(i<=scy1)
      for (j=1; j<=na; j++) sel(j,i) = sel92(j);
//	 else if(scy1+1 <= i <= scy2)
	  else
	  for (j=1; j<=na; j++) sel(j,i) = sel93(j);
//     else
//	  for (j=1; j<=na; j++) sel(j,i) = sel08(j);	
    }

//the first 2 length groups molt 100% and the last 2 length groups are 100% selected.
   mp1 = mp1/mp1(1);
   mp1(1) = 1.0;
//   mp2 = mp2/mp2(1);
//   mp2(1) = 1.0;
//   mp3 = mp3/mp3(1);
//   mp3(1) = 1.0;
   sos(na) = sos0;
   selw(1) = 1.0/(1.0+mfexp(-1.0*sw1*(78.5-sw2)));
   selw(2) = 1.0/(1.0+mfexp(-1.0*sw1*(88.5-sw2)));
   selw(na) = sw3;
   selt(1) = 1.0/(1.0+mfexp(-1.0*st1*(78.5-st2)));
   selt(2) = 1.0/(1.0+mfexp(-1.0*st1*(88.5-st2)));
 //  selt(3) = 1.0/(1.0+mfexp(-1.0*st1*(98.5-st2)));
 //  selt(4) = 1.0/(1.0+mfexp(-1.0*st1*(108.5-st2)));
   selt(na) = slt6;
   selp(1) = 1.0/(1.0+mfexp(-1.0*sp1*(78.5-sp2)));
   selp(2) = 1.0/(1.0+mfexp(-1.0*sp1*(88.5-sp2)));

// =========================================================================

FUNCTION get_first_year_abundance
  int i,j;
//  dvariable first_y(1,na);

  log_rec = log_relrec+log_recscale;
  first_y = mfexp(log_initpop);
  for (i=1;i<ny;i++) rec(i) = mfexp(log_rec(i));

   // recruits in the last year are assumed to be average of the most recent five years
   rec(ny-1) = 0.2*(rec(ny-2)+rec(ny-3)+rec(ny-4)+rec(ny-5)+rec(ny-6));
   //use the length proportion in the first year to estimate the abundance.
   tb(1) = 0.0;
   ops(1) = 0.0;
   ops(2) = 0.0;
   for (j=1;j<=na;j++)
    {
//   if (selt(j) < 0.1) selt(j) = 0.1;
//   if (j <= 2)
//    nps(j,1) = ont(j,1)*first_y/selt(j);
//   else
//    {
//     nps(j,1) = (ont(j,1)-oot(j,1))*first_y/selt(j);
//     ops(j,1) = oot(j,1)*first_y/selt(j);
//    }
     ops(j,1) = 0.0;
     nps(j,1) = first_y(j);
     tb(1) += (nps(j,1)+ops(j,1)*sos(j))*sel(j,1)*lg(j);
    }
   if (tb(1) < 0.001) tb(1) = 0.001;

   //estimated length proportion for summer commercial fishery.
   for (j=1;j<=na;j++)
    {
     enc(j,1) = nps(j,1)*sel(j,1)*lg(j)/tb(1);
     eoc(j,1) = ops(j,1)*sos(j)*sel(j,1)*lg(j)/tb(1);
    }

// =========================================================================

FUNCTION get_number_by_size
  int i,j,k;
  dvariable pp,tt1, ttt0;
  dvar_vector tt0(1,na);

// Estimate projected crab abundance
  for (i=1;i<ny;i++)
  {
//    if (i < 7)
      mp0 = mp1;
//    else
//     if (i < 16)
//      mp0 = mp2;
//     else
//      mp0 = mp3;

// Calculate survivors from the summer fishery and discareds
// ttt0: total number of legal crab  bc: discards biomass
    ttt0 = 0.0;
    bc(i) = 0.0;	
    for (k=1;k<=na;k++) {ttt0 += (nps(k,i)+ops(k,i))*lg(k);}	
// if projected abundance is lower than actual catch then projected cach should be adjusted to 1.2 times of actual harvest
    if (ttt0 < tc(i)) { ttt0 = tc(i)*1.2;}
// if projected abundance is lower than 0.001 thaen  projected cach should be adjusted to 0.001
    else if (ttt0 < 0.001) ttt0 = 0.001;

    for (k=1;k<=na;k++)
     {
      tt0(k) = ((nps(k,i)+ops(k,i))*exp(-ys(i)*Mn(k))-tc(i)*(enc(k,i)+eoc(k,i))-(tc(i)/ttt0)*(nps(k,i)+ops(k,i))*sel(k,i)*(1.0-lg(k))*hm)*exp(-(0.583-ys(i))*Mn(k));
      if (tt0(k) < 0.0) tt0(k)= 0.0;
      bc(i) += (tc(i)/ttt0)*(nps(k,i)+ops(k,i))*sel(k,i)*(1.0-lg(k))*hm*wm(k);
     }
	// Commercial accept only > CW 5 inch since 2005 and sub comercial size crabs are assumed discarded
	if (i >= scp) bc(i) += (tc(i)/ttt0)*(nps(4,i)+ops(4,i))*sel(4,i)*(lg(k)*(1-p4))*hm*wm(4);
	
    // Dynamics
    for (j=1;j<=na;j++)
     {
      pp = 0;
      for (k=1;k<=j;k++) pp += gr(k,j)*tt0(k)*mp0(k);
      npw(j,i) = pp;
// Unobserved juvenile crab abundance was estimated
      if (j==1) npw(j,i) += r1*rec(i);
      if (j==2) npw(j,i) += (1.0-r1)*rec(i);
      opw(j,i) = tt0(j)*(1.0-mp0(j));
      if (opw(j,i)<0.0) opw(j,i) = 0.0001;
     }

    //propotions for winter fishery and subsistence fishery
    pp = 0.0; tt1 = 0.0;
    for (j=3;j<=na;j++)
     {
      pp += (npw(j,i)+opw(j,i)*sow)*selw(j);
      tt1 += (npw(j,i)+opw(j,i)*sow)*selw(j)*lg(j);
     }
    if (pp <= 0.0) pp = 0.000001;
    if (tt1 <= 0.0) tt1 = 0.000001;
    for (j=3;j<=na;j++)
     {
      en0(j,i) = npw(j,i)*selw(j)/pp;
      eo0(j,i) = opw(j,i)*sow*selw(j)/pp;
      en1(j,i) = npw(j,i)*selw(j)*lg(j)/tt1;
      eo1(j,i) = opw(j,i)*sow*selw(j)*lg(j)/tt1;
     }

    // Extract exploitation rate
    FSum(i) = tc(i)/tb(i);
    FSub(i) = tsc(i)/pp;
    FWin(i) = twc(i)/tt1;

    //update to abundance in next summer (new shell and old shell)
    for (j=1;j<=na;j++)
	{
     nps(j,i+1) = (npw(j,i)-en0(j,i)*tsc(i)-en1(j,i)*twc(i))*exp(-0.417*Mn(j));
     ops(j,i+1) = (opw(j,i)-eo0(j,i)*tsc(i)-eo1(j,i)*twc(i))*exp(-0.417*Mo(j));
	}

    //proportion for summer fishery. Higher selectivity for new-shell than for old-shell.
	// tb: mean exploitable leagal abundance
	tb(i+1) = 0.0;
    for (j=1;j<=na;j++)
     tb(i+1) += (nps(j,i+1)+ops(j,i+1)*sos(j))*sel(j,i+1)*lg(j);	
// Commercial accept only > CW 5 inch since 2005  	
	if (i >= scp) tb(i+1) += (p4-1)*(nps(4,i+1)+ops(4,i+1)*sos(4))*sel(4,i+1)*lg(4);
	
	// mean exploitable leagal abundance does not go neative
    if (tb(i+1) < 0.001) tb(i+1) = 0.001;
	// Calculate proprotion of newshell and oldshell by comm fiish 	
	 for (j=1;j<=na;j++)
     {
	  enc(j,i+1) = nps(j,i+1)*sel(j,i+1)*lg(j)/tb(i+1);
      eoc(j,i+1) = ops(j,i+1)*sos(j)*sel(j,i+1)*lg(j)/tb(i+1);
     }
	// Since 2005 Commercial accept only > CW 5 inch crab, which affects only 4th crab size.	 	
	if (i >= scp)
		{
		enc(4,i+1) = nps(4,i+1)*sel(4,i+1)*lg(4)*p4/tb(i+1);
		eoc(4,i+1) = ops(4,i+1)*sos(4)*sel(4,i+1)*lg(4)*p4/tb(i+1);
		}
   }

  // if at the last_phase, calculated abundances in last year & legal abundance for all years.
  last_y = column(nps,ny);
  for (j=3; j<=na; j++) last_y(j) += ops(j,ny);

  for (i=1; i<=ny; i++)
   {
    legal(i) = 0.0; legalb(i) = 0.0;
    for (j=1; j<=na; j++)
    {
       legal(i) += (nps(j,i)+ops(j,i))*lg(j);
       legalb(i) += (nps(j,i)+ops(j,i))*lg(j)*wm(j);
    }
    mmb(i) = 0.0;
    for (j=3; j<=na; j++) mmb(i) += (nps(j,i)+ops(j,i))*wm(j);
    legaln(i) = legal(i);
    mmb0(i) = mmb(i);
   }
  // last_legal = abundance*proportion of legal*catchability*average weight
//  last_legal = last_y(3)*lg(3)*wm(3)*sel08(3)+last_y(4)*lg(4)*wm(4)*sel08(4)+last_y(5)*lg(5)*wm(5)*sel08(5)+last_y(6)*lg(6)*wm(6)*sel08(6);
  last_legal = last_y(3)*lg(3)*wm(3)*sel93(3)+last_y(4)*lg(4)*p4*wm(4)*sel93(4)+last_y(5)*lg(5)*wm(5)*sel93(5)+last_y(6)*lg(6)*wm(6)*sel93(6);
  last_subl = last_y(1)*(1-lg(1))*wm(1)*sel93(1)+last_y(2)*(1-lg(2))*wm(2)*sel93(2)+last_y(3)*(1-lg(3))*wm(3)*sel93(3)+last_y(4)*(1-lg(4))*wm(4)*sel93(4)+last_y(5)*(1-lg(5))*wm(5)*sel93(5)+last_y(6)*(1-lg(6))*wm(6)*sel93(6);
  last_ofl = last_legal*(1.0-exp(-M));
  last_mmb   = mmb(ny);

// ==========================================================================

FUNCTION get_proportion_and_effort
  int i,j;
  dvariable tf8,tf10,pp;

  ett.initialize();
  ett1.initialize();
  etp.initialize();
  etp1.initialize();
  eop.initialize();
  eof.initialize();
  eow.initialize();
  eoo.initialize();
  ecpue.initialize();


  //proprotion and total abundance for trawl survey
  for (i=1;i<=nyt;i++)
   {
    if (yt(i) > ys(it(i)))
     {
      tf8 = ys(it(i));              //time lag from July 1 to fishery
      tf10 = yt(i) - ys(it(i));     //time lag from fishery to survey
     }
    else
     {
      tf8 = yt(i);
      tf10 = 0.0;
     }
    for (j=1;j<=na;j++)
     {
      et0(j,i) = ((nps(j,it(i))+ops(j,it(i)))*exp(-tf8*Mn(j))-(enc(j,it(i))+eoc(j,it(i)))*pct(i)*tc(it(i)))*exp(-tf10*Mn(j))*selt(j);
      if (et0(j,i) < 0.0) et0(j,i) = 0.0;
      ett(i) += et0(j,i);
      if (j<4) ett1(i) += et0(j,i);
     }
    if (ett(i) <= 0.0) ett(i) = 0.000001;
    if (ett1(i) <= 0.0) ett1(i) = 0.000001;
    ett2(i) = ett(i) - ett1(i);
    for (j=1;j<=na;j++)
     {
      ent(j,i) = et0(j,i)/ett(i);
      if (ent(j,i) < 0.0) ent(j,i) = 0.0;
     }
   }

  //proprotion and total abundance for pot survey
  for (i=1;i<=nyp;i++)
   {
    for (j=1;j<=na;j++)
     {
      etp(i) += (nps(j,ip(i))+ops(j,ip(i)))*exp(-yp(i)*Mn(j))*selp(j);
      if (j < 4) etp1(i) += (nps(j,ip(i))+ops(j,ip(i)))*exp(-yp(i)*Mn(j))*selp(j);
     }
    if (etp(i) <= 0.0) etp(i) = 0.000001;
    if (etp1(i) <= 0.0) etp1(i) = 0.000001;
    etp2(i) = etp(i) - etp1(i);
    for (j=1;j<=na;j++)
     {
      enp(j,i) = nps(j,ip(i))*exp(-yp(i)*Mn(j))*selp(j)/etp(i);
      if (j > 1) eop(j,i) = ops(j,ip(i))*exp(-yp(i)*Mn(j))*selp(j)/etp(i);
     }
   }

  //proprotion for winter project, higher selectivity for old-shell than for new-shell.
  for (i=1;i<=nyw;i++)
   {
    pp = 0.0;
    for (j=1; j<=na; j++)
     pp += (npw(j,iw(i))+opw(j,iw(i))*sow)*selw(j);
    if (pp <= 0.0) pp = 0.000001;
    for (j=1;j<=na;j++)
     {
      enw(j,i) = npw(j,iw(i))*selw(j)/pp;
      if (j > 1) eow(j,i) = opw(j,iw(i))*selw(j)*sow/pp;
     }
   }

  //proprotion for observer's data
  for (i=1;i<=nyo;i++)
   {
    pp = 0.0;
    for (j=1;j<=na;j++)
     pp += (nps(j,io(i))+ops(j,io(i)))*sel(j,io(i))*(1.0-lg(j));
    if (pp <= 0.0) pp = 0.000001;
    for (j=1;j<=na;j++)
     {
      eno(j,i) = nps(j,io(i))*sel(j,io(i))*(1.0-lg(j))/pp;
      if (j > 1) eoo(j,i) = ops(j,io(i))*sel(j,io(i))*(1.0-lg(j))/pp;
     }
   }

  //proprotion for pre-fishery survey
  pp = 0.0;
  for (j=1;j<=na;j++)
   pp += (nps(j,ipre)+ops(j,ipre))*sel(j,ipre);
  if (pp <= 0.0) pp = 0.000001;
  for (j=1;j<=na;j++)
   {
    enf(j) = nps(j,ipre)*sel(j,ipre)/pp;
    if (j > 1) eof(j) = ops(j,ipre)*sel(j,ipre)/pp;
   }

//fishing effort and cpue for summer fishery
  for (i=1;i<=ny;i++)
  {
    tb(i) = tb(i) - 0.5*tc(i);  //exploitable abundance at the middle of the season
    if (tb(i) < 0.001) tb(i) = 0.001;
    if (i <= scy1)
	 {
	 ecpue(i) = q1*tb(i);
     ete(i) = tc(i)/(q1*tb(i));     //before 1993
	 }
//	else if (scy1 < i < scy2)
	else
	 {
	 ete(i) = tc(i)/(q2*tb(i));     //after 1993
	 ecpue(i) = q2*tb(i);	
	 }
//	else
//	 {	
//	 ete(i) = tc(i)/(q3*tb(i));     //after 2008
//	 ecpue(i) = q3*tb(i);
//	 }
  }

//fishing cpue: remove data
  for (i=1;i<=ny;i++)
   {
   if (stcpue(i) <= 0.0) {ecpue(i) = 0.0;}
   }

	 //Calculate summer fishery 3 year average cpue


// =======================================================================

FUNCTION evaluate_the_objective_function

//weight each year estimate by 1/(2*variance) - use cv of biomass in sqrt(log(cv^2+1)) as sd of log(biomass)
//  tf1 = norm2(elem_div((log(tt1+1.e-3)-log(ett1+1.e-3)),
//               sqrt(2)*sqrt(log(elem_prod(cv1,cv1)+1.0))));
//  f += tf1;
//  tf2 = norm2(elem_div((log(tt2+1.e-3)-log(ett2+1.e-3)),
//               sqrt(2)*sqrt(log(elem_prod(cv2,cv2)+1.0))));

// Log likelihodd for trawl survey
  tf2 = 0.5*norm2(elem_div((log(tt+1.e-3)-log(qt*(ett1+ett2)+1.e-3)),sqrt(log(elem_prod(cv3,cv3)+1.0))));
  f = tf2;
// Log likelihodd crab pot suervey
  tf3 = 0.5*norm2(log(tp+1.e-3)-log(qt*(etp1+etp2)+1.e-3))/log(spcv*spcv+1.0);
  f += tf3;

  // Log likelihodd legal crab pot suervey
//  tf4 = norm2(log(tp2+1.e-3)-log(etp2+1.e-3))/(2.0*log(0.27*0.27+1.0));
//  f += tf4;
// Log likelihood effort in summer fishery
//   tf5 = norm2(log(te+1.e-3)-log(ete+1.e-3));
//  f += lamc*tf5;

// Log likelihood standard cpue in summer fishery
   tf5 = 0.5*norm2(elem_div((log(stcpue+1.e-3)-log(1000*ecpue+1.e-3)),sqrt(log(elem_prod(cvcpue,cvcpue)+1.0))));
   f += lamc*tf5;

//Log likelhihood for trawl survey multinomial proportion
  tf6 = sum(elem_prod(st,colsum(elem_prod(ont,log(ent+1.e-3)))));
  f -=tf6;
//Log likelhihood for pot survey multinomial proportion
  tf7 = sum(elem_prod(sp,colsum(elem_prod(onp+oop,log(enp+eop+1.e-3)))));
  f -= tf7;
//Log likelhihood for winter multinomial proportion
  tf9 = sum(elem_prod(sw,colsum(elem_prod(onw+oow,log(enw+eow+1.e-3)))));
  f -= tf9;

//Log likelhihood size proportion for summer fishery survey
  tf11 = sum(elem_prod(sc,colsum(elem_prod(onc+ooc,log(enc+eoc+1.e-3)))));
  f -= tf11;

  tf13 = 0.01*norm2(log_relrec);                         //deviation in recruits.
  f += tf13;

 //Log likelhihood for pre-fishery survey propotion
//  tf14 = pre*(sum(elem_prod(onf+oof,log(enf+eof+1.e-3)))+sum(elem_prod(onf+oof,log(enf+eof+1.e-3))));
//  f -= tf14;
//Log likelhihood size proportion for observer survey
  tf15 = sum(elem_prod(so,colsum(elem_prod(ono+ooo,log(eno+eoo+1.e-3)))));
  f -= tf15;

  if (f < fbest)
   {
    fbest = f;
    cout <<"Tr_imma Tr-Mat Pot1  Pot2  effort  Tr_pr  Pot_pr  Win_pr  Sum_pr "<<endl;
    cout <<tf1<<" "<<tf2<<" "<<tf3<<" "<<tf4<<" "<<tf5<<" "<<tf6<<" "<<tf7<<" "<<tf9<<" "<<tf11<<endl;
    cout <<"Recruit deviation Observer's data total F"<<endl;
    cout <<tf13<<" "<<tf15<<" "<<f<<endl;
    cout<<"Legal in ny = "<<legal(ny)/1000.0<<" mmb ny = "<<mmb(ny)/1000.0<<endl;
   }

// ==========================================================================

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

// ==========================================================================

FUNCTION Find_OFL
  dvariable Bmsy,Fmsy;
  int ii,iy;

  // Define Bmsy
  Bmsy = 0;
  for (iy=Bmsy_yr1;iy<=Bmsy_yr2;iy++) Bmsy += mmb(iy);
  Bmsy /= float(Bmsy_yr2-Bmsy_yr1+1);

  // Specify recruitment
  RecPass = 0;
  for (iy=YrPass-5;iy<=YrPass-1;iy++) RecPass += rec(iy);
  RecPass /= 5;

  Fmsy = M;

  Fmult = Fmsy;
  Proj1Yr();
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

// ===========================================================================

FUNCTION Proj1Yr
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
   Total += (nps(Len,YrPass)+ops(Len,YrPass)*sos(Len))*sel(Len,ny)*lg(Len);
  if (Total < 0.01) Total = 0.001;
  CatSum = 0;
  for (Len=1;Len<=na;Len++)
   {
    encf(Len) = nps(Len,YrPass)*sel(Len,ny)*lg(Len)/Total;
    CatNum = nps(Len,YrPass)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
    CatSum += CatNum;
    TotalCatch += CatNum*wm(Len);
    eocf(Len) = ops(Len,YrPass)*sos(Len)*sel(Len,ny)*lg(Len)/Total;
    CatNum = ops(Len,YrPass)*sos(Len)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
    CatSum += CatNum;
    TotalCatch += CatNum*wm(Len);
   }

  // Survivors from the summer fishery
  for (Len=1;Len<=na;Len++)
   {
    tt0(Len) = ((nps(Len,YrPass)+ops(Len,YrPass))*exp(-yss*Mn(Len))-CatSum*(encf(Len)+eocf(Len)))*exp(-(0.583-yss)*Mn(Len));
    if (tt0(Len) < 0.0) tt0(Len)= 0.0;
   }

  // Upodate dynamics (after summer fishery)
  for (Len=1;Len<=na;Len++)
   {
    pp = 0;
    for (Len1=1;Len1<=Len;Len1++) pp += gr(Len1,Len)*tt0(Len1)*mp0(Len1);
    npw(Len,YrPass) = pp;
    if (Len==1) npw(1,YrPass) += r1*RecPass;
    if (Len==2) npw(1,YrPass) += (1.0-r1)*RecPass;
    opw(Len,YrPass) = tt0(Len)*(1.0-mp0(Len));
    if (opw(Len,YrPass)<0.0) opw(Len,YrPass) = 0.0001;
   }

  //propotions for winter fishery and subsistence fishery
  pp = 0.0; tt1 = 0.0;
  for (Len=3;Len<=na;Len++)
   {
    pp += (npw(Len,YrPass)+opw(Len,YrPass)*sow)*selw(Len);
    tt1 += (npw(Len,YrPass)+opw(Len,YrPass)*sow)*selw(Len)*lg(Len);
   }
  if (pp <= 0.0) pp = 0.000001;
  if (tt1 <= 0.0) tt1 = 0.000001;

  CatSub = 0; CatWin = 0;
  en0f.initialize();
  eo1f.initialize();
  for (Len=3;Len<=na;Len++)
   {
    en0f(Len) = npw(Len,YrPass)*selw(Len)/pp;
    CatNum = npw(Len,YrPass)*selw(Len)*Fmult*FratSub;
    CatSub += CatNum;
    TotalCatch += CatNum*wm(Len);
    eo0f(Len) = opw(Len,YrPass)*sow*selw(Len)/pp;
    CatNum = opw(Len,YrPass)*sow*selw(Len)*Fmult*FratSub;
    CatSub += CatNum;
    TotalCatch += CatNum*wm(Len);
    en1f(Len) = npw(Len,YrPass)*selw(Len)*lg(Len)/tt1;
    CatNum = npw(Len,YrPass)*selw(Len)*lg(Len)*Fmult*FratWin;
    CatWin += CatNum;
    TotalCatch += CatNum*wm(Len);
    eo1f(Len) = opw(Len,YrPass)*sow*selw(Len)*lg(Len)/tt1;
    CatNum = opw(Len,YrPass)*sow*selw(Len)*lg(Len)*Fmult*FratWin;
    CatWin += CatNum;
    TotalCatch += CatNum*wm(Len);
   }

  //update to abundance in next summer (new shell)
  for (Len=1;Len<=na;Len++)
   nps(Len,YrPass+1) = (npw(Len,YrPass)-en0f(Len)*CatSub-en1f(Len)*CatWin)*exp(-0.417*Mn(Len));
  //update to abundance in next summer (old shell)
  for (Len=1;Len<=na;Len++)
   ops(Len,YrPass+1) = (opw(Len,YrPass)-eo0f(Len)*CatSub-eo1f(Len)*CatWin)*exp(-0.417*Mo(Len));

  // Compute legal and MMB
  legal(YrPass+1) = 0.0;
  for (Len=1; Len<=na; Len++) legal(YrPass+1) += (nps(Len,YrPass+1)+ops(Len,YrPass+1))*lg(Len);
  mmb(YrPass+1) = 0.0;
  for (Len=3; Len<=na; Len++) mmb(YrPass+1) += (nps(Len,YrPass+1)+ops(Len,YrPass+1))*wm(Len);

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
   Total += (npsEst(Len,YrPass)+opsEst(Len,YrPass)*sos(Len))*sel(Len,ny)*lg(Len);
  if (Total < 0.01) Total = 0.001;
  CatSum = 0;
  for (Len=1;Len<=na;Len++)
   {
    encf(Len) = npsEst(Len,YrPass)*sel(Len,ny)*lg(Len)/Total;
    CatNum = npsEst(Len,YrPass)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
    CatSum += CatNum;
    TotalCatch += CatNum*wm(Len);
    eocf(Len) = opsEst(Len,YrPass)*sos(Len)*sel(Len,ny)*lg(Len)/Total;
    CatNum = opsEst(Len,YrPass)*sos(Len)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
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
       Total += (nps(Len,FutYr)+ops(Len,FutYr)*sos(Len))*sel(Len,ny)*lg(Len);
      if (Total < 0.01) Total = 0.001;
      FutCatSum = 0;
      for (Len=1;Len<=na;Len++)
       {
        encf(Len) = nps(Len,FutYr)*sel(Len,ny)*lg(Len)/Total;
        CatNum = nps(Len,FutYr)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*FutMort;
        FutCatSum += CatNum;
        PredCatchW += CatNum*wm(Len);
        eocf(Len) = ops(Len,FutYr)*sos(Len)*sel(Len,ny)*lg(Len)/Total;
        CatNum = ops(Len,FutYr)*sos(Len)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*FutMort;
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
     Total += (nps(Len,FutYr)+ops(Len,FutYr)*sos(Len))*sel(Len,ny)*lg(Len);
    if (Total < 0.01) Total = 0.001;
    CatSum = 0;
    for (Len=1;Len<=na;Len++)
     {
      encf(Len) = nps(Len,FutYr)*sel(Len,ny)*lg(Len)/Total;
      CatNum = nps(Len,FutYr)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
      CatSum += CatNum;
      TotalCatch += CatNum*wm(Len);
      eocf(Len) = ops(Len,FutYr)*sos(Len)*sel(Len,ny)*lg(Len)/Total;
      CatNum = ops(Len,FutYr)*sos(Len)*sel(Len,ny)*lg(Len)*exp(-yss*Mn(Len))*Fmult;
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
  dvar_matrix ref_catch_w_n(1,na,1,100);
  dvar_matrix ref_catch_w_o(1,na,1,100);
  dvar_matrix ref_catch(1,na,1,100);
  dvar_matrix ref_ps(1,na,1,100);
  dvar_matrix ref_nps(1,na,1,100);
  dvar_matrix ref_ops(1,na,1,100);
  dvar_matrix ref_npw(1,na,1,100);
  dvar_matrix ref_opw(1,na,1,100);
  dvar_vector ref_catch_b(1,100);
  dvar_vector ref_catch_n(1,100);
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
  for (i = ny-10; i<= ny-1; i++)
   {
    tt0 = 0.0;
    for (j=3;j<=na;j++) tt0 += (npw(j,i)+opw(j,i)*sow)*selw(j)*lg(j);
    if (tt0 <= 0.0) tt0 = 0.000001;
    ref_Htw(i) = (twc(i)+tsc(i))/tt0;
    ref_Hwc += ref_Htw(i);
   }
  ref_Hwc = ref_Hwc/10.0;        //use average harvest rate of the most recent 10 years as winter and subsistence catch
  ref_Hts = 0.0;
  for (j = 1; j<= 101; j++)
   {
    ref_F = 0.01*j-0.01;
    if (j > 1) ref_Hts = 1.0-mfexp(-1.0*ref_F);
    ref_catch.initialize();
    ref_nps.initialize();
    ref_ops.initialize();

    //initial year
    for (k=1; k<=na; k++)
     {
      ref_nps(k,1) = nps(k,ny);
      ref_ops(k,1) = ops(k,ny);
      ref_ps(k,1)  = ref_nps(k,1)  + ref_ops(k,1);
     }

    //numbers at length
    for (i=1;i<100;i++)
     {
      for (k=1;k<=na;k++)
       {
        ref_catch_s(k,i) = (ref_nps(k,i)+ref_ops(k,i)*sos(k))*sel(k,ny-1)*lg(k)*ref_Hts;
        pp = 0;
        for (kk=1;kk<=k;kk++)
        {
          tt0 = (ref_ps(kk,i)*mfexp(-ref_ys*Mn(kk))-ref_catch_s(kk,i)-ref_Hts*ref_ps(kk,i)*sel(kk,ny-1)*(1.0-lg(kk))*hm)*mfexp(-(0.583-ref_ys)*Mn(kk));
          if (tt0<0.0) tt0 = 0.0;
          pp += gr(kk,k)*tt0*mp0(kk);
        }
        ref_npw(k,i) = pp;
        if (k==1) ref_npw(k,i) += r1*1000.0;
        if (k==2) ref_npw(k,i) += (1.0-r1)*1000.0;
        ref_opw(k,i) = (ref_ps(k,i)*mfexp(-ref_ys*Mn(k))-ref_catch_s(k,i)-ref_Hts*ref_ps(k,i)*sel(k,ny-1)*(1.0-lg(k))*hm)*mfexp(-(0.583-ref_ys)*Mn(k))*(1.0-mp0(k));
        if (ref_opw(k,i)<0.0) ref_opw(k,i) = 0.0001;
       }
      //update to abundance in next summer
      for (k=1;k<=na;k++)
       {
        if (j == 1)                            //all fishing is closed.
         {
          ref_catch_w_n(k,i) = 0.0;
          ref_catch_w_o(k,i) = 0.0;
         }
        else
         {
          ref_catch_w_n(k,i) = ref_npw(k,i)*selw(k)*lg(k)*ref_Hwc;
          ref_catch_w_o(k,i) = ref_opw(k,i)*sow*selw(k)*lg(k)*ref_Hwc;
         }
        ref_catch(k,i) = ref_catch_s(k,i) + ref_catch_w_n(k,i) + ref_catch_w_o(k,i);
       }
      for (k=1;k<=na;k++)
       ref_ops(k,i+1) = (ref_opw(k,i)-ref_catch_w_o(k,i))*mfexp(-0.417*Mo(k));
      for (k=1;k<=na;k++)
       {
        ref_nps(k,i+1) = (ref_npw(k,i)-ref_catch_w_n(k,i))*mfexp(-0.417*Mn(k));
        ref_ps(k,i+1)  = ref_nps(k,i+1)  + ref_ops(k,i+1);
       }
      ref_legal(i) = 0.0;
      ref_legal_b(i) = 0.0;
      ref_mmb(i) = 0.0;
      ref_catch_b(i) = 0.0;
      ref_catch_n(i) = 0.0;
      for (k=1; k<=na; k++)
       {
        ref_legal(i) += ref_ps(k,i)*lg(k);
        ref_legal_b(i) += ref_ps(k,i)*lg(k)*wm(k);
        ref_catch_b(i) += ref_catch(k,i)*wm(k);
        ref_catch_n(i) += ref_catch(k,i);
        if (k > 2) ref_mmb(i) += ref_ps(k,i)*wm(k);
       }
     }
    ref_mbio(j) = ref_mmb(99)/1000.0;                    //lbs/R
    ref_totc(j) = ref_catch_b(99)/1000.0;
    ref_legaln(j) = ref_legal(99)/1000.0;
    ref_Hr(j) = (ref_catch_n(99)/1000.0)/ref_legaln(j);
   }
  i35 = 0; i40 = 0;
  for (j = 1; j<= 101; j++)
   {
    if (i35 < 1.0)
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
  report1 <<"Harvest rate in term of legal males 7/1 as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_Hr<<endl;
  report1 <<"Total MMB (>93mm) 7/1 as F = 0.00, 0.01, ... 1.0 (lbs/R)"<<endl;
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

REPORT_SECTION
  cout << "Report Section" << endl;

  get_reference_points();

  FratSub = (FSub(ny-5)+FSub(ny-4)+FSub(ny-3)+FSub(ny-2)+FSub(ny-1))/(FSum(ny-5)+FSum(ny-4)+FSum(ny-3)+FSum(ny-2)+FSum(ny-1));
  FratWin = (FWin(ny-5)+FWin(ny-4)+FWin(ny-3)+FWin(ny-2)+FWin(ny-1))/(FSum(ny-5)+FSum(ny-4)+FSum(ny-3)+FSum(ny-2)+FSum(ny-1));

  report << "nps" << endl << nps << endl;  // Modeled new shell summer
  report << "ops" << endl << ops << endl;  // Modeled old shell summer
//  report << "npw" << endl << npw << endl;  // Modeled new shell winter
//  report << "opw" << endl << opw << endl;  // Modeled oldshell shell winter
  report << "ett" << endl << ett << endl;  // Estimated trawl abundance
  report << "etp" << endl << etp << endl;  // Estimated Pot abundance
  report << "ete" << endl << ete << endl;  // Estimated Summer fishery effort
  report << "ecpue" << endl << ecpue << endl;  // Estimated Summer fishery cpue
//  report << "ecp3" << endl << ecp3 << endl;  // Estimated Summer fishery cpue
  report << "ent" << endl << ent << endl;  // Estimated trawl size proportion
  report << "enp" << endl << enp << endl;  // Estimated pot newshell size proportion
  report << "eop" << endl << eop << endl;  // Estimated pot oldshell size proportion
  report << "enw" << endl << enw << endl;  // Estimated winter survey newshell size proportion
  report << "eow" << endl << eow << endl;  // Estimated winter survey newshell size proportion
  report << "enf" << endl << enf << endl;  // Estimated pre-fishery newshell size proportion
  report << "eof" << endl << eof << endl;  // Estimated pre-fishery oldshell size proportion
  report << "enc" << endl << enc << endl;  // Estimated summer fishery newshell size proportion
  report << "eoc" << endl << eoc << endl;  // Estimated summer fishery oldshell size proportion
  report << "eno" << endl << eno << endl;  // Estimated observer newshell size proportion
  report << "eoo" << endl << eoo << endl;  // Estimated observer oldshell size proportion
  report << "Rec" << endl << rec << endl;  // Estimated Recruits abundance
  report << "Legal" << endl << legal << endl;  // Estimated legal abundance
  report << "MMB" << endl << mmb << endl;  // Estimated mmb abundance
  report << "ett1" << endl << ett1 << endl;  // Estimated non-legal trawl survey abundance
  report << "ett2" << endl << ett2 << endl;  // Estimated legal trawl survey abundance
  report << "etp1" << endl << etp1 << endl;  // Estimated non-legal Pot survey abundance
  report << "etp2" << endl << etp2 << endl;  // Estimated legal Pot survey abundance
  report << "bc" << endl << bc << endl;  // Estimated discards biomass
  report << "f" << endl << f << endl;  // Total likelihood
  // Individual likelihood
  report << "tf1" << endl << tf1 << " " << tf2 << " " << tf3 << " " << tf4 << " " << tf5 << " " << tf6 << " " << tf7 << " " << tf8 << " " << tf9 << " " << tf10 << " " << tf11 << " " << tf12 << " " << tf13 << " " << tf14 << " " << tf15 << endl;

GLOBALS_SECTION
  #include <math.h>
  #include <admodel.h>
  #include <time.h>

  ofstream mcmc1,mcmc2;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

// ===========================================================================

TOP_OF_MAIN_SECTION
  arrmblsize = 10000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);

  time(&start);
  mcmc1.open("Pstar.Out");
  mcmc2.open("FiveY.Out");

// ===========================================================================

FINAL_SECTION

 // Output summary stuff
 time(&finish);
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;


