//===============================================================================================================
// NS3n2016_Feb01_6c.TPL
// This model is identical with NNS3n2016_Feb01.TPL.  
// Model codes were cleaned and simplifiled, all unnecessary data were removed 
// This version eliminate q2 and sel93() of commercial catch
// q3 and sel08() represents fishery 1993 after
// p4 was assumed 1.0 (i.e. change of commercial catch size did not affect commercial catch) 
// Incorporated discards of Winter commercial and subsistence 
//===============================================================================================================

//===============================================================================================================
// 1.0  Data Entry 
//===============================================================================================================
DATA_SECTION
  init_int fyear              //first year
  init_int lyear              //last year  
  init_int na                 //number of length groups
  init_int nyt                //number of years with trawl survey
  init_int nyw                //number of years with winter project
  init_int nyo                //number of years with observer's data
  init_ivector scy(1,2)           //summer commercial fishery event indices: 1)year large trawl fishery ended, 2) escape mechanism installed
  init_int scp                //year summer commrecial fishery CW>5 inch crab was accepted by buyers
  init_number slm             //mid-size of smallest length group (mm)
  init_number slt             //length interval (mm)
  init_number M2              //instantaneous natural mortality
  init_number ms6			  //instantaneous natural mortality multiplier for the last size class
  init_vector maxs(1,2)       //maximum effective sample size for length proportion
  init_vector efn(1,2)           //% of effective sample size due to multin.distribution
  init_vector hm(1,2)              //handling mortality rate for bycatch summer and winter   
  init_vector lg(1,na)        //proportion of legals by length group
  init_vector wm(1,na)        //Mean weight 
  init_ivector it(1,nyt)      //year id for trawl survey data
  init_vector tt(1,nyt)       //total annual abundance from trawl survey
  init_vector pct(1,nyt)      //% of summer catch occurred before the mid point of survey
  init_vector yt(1,nyt)       //Mid point of trawl survey from July 1
  init_vector cv(1,nyt)       //annual cv of totals for trawl survey
  init_matrix nont(1,na,1,nyt) //number of length group for new-shell, trawl survey
  init_matrix noot(1,na,1,nyt) //number of length group for old-shell, trawl survey  
  init_ivector iw(1,nyw)      //year id for winter project
  init_vector wcpue(1,nyw)   //Winter pot survey cpue  
  init_matrix nonw(1,na,1,nyw) //number of length group for new-shell, winter survey
  init_matrix noow(1,na,1,nyw) //number of length group for old-shell, winter survey
  init_vector tc(fyear,lyear)        //annual catch for summer fishery
  init_vector te(fyear,lyear)        //annual effort for summer fishery
  init_vector stcpue(fyear,lyear)    //annual standardized cpue for summer fishery
  init_vector secpue(fyear,lyear)    //se of  annual cpue for summer fishery
  init_vector ys(fyear,lyear)        //Mid point of summer fishery from July 1
  init_matrix nonc(1,na,fyear,lyear)  //proportion of length group for new-shell, summer fishery
  init_matrix nooc(1,na,fyear,lyear)  //proportion of length group for old-shell, summer fishery
  init_vector twc(fyear,lyear)     //annual catch for winter fishery
  init_vector tws(fyear,lyear)     //annual retained catch for winter subsistence fishery
  init_vector twst(fyear,lyear)     //annual total catch for winter subsistence fishery
  init_ivector io(1,nyo)      //year id for observer's data
  init_matrix nono(1,na,1,nyo) //number of length group for new-shell, observer survey
  init_matrix nooo(1,na,1,nyo) //number of length group for old-shell, observer survey
  init_imatrix tag_recov1(1,21,1,8)   // tag recovery data  1976-1992 
  init_imatrix tag_recov2(1,21,1,8)   // tag recovery data  1993-present    
  init_number SDRec           // SD of recruitment 
  init_number SDW             // SD of winter cpue    
  init_number log_initp       // initial population control
  init_int initp_phase      // determine phase
  init_int qtno_phase         // Survey Q phase control vector (NOAA) 
  init_int M_phase            // M estimate
  init_int ms_phase           // ms estimate
  init_int lamc               // likelifhood for cpue
  init_int lamw               // likelihood for winter cpue
  init_int lawp               // likelihood weight for winter pot comp
  init_number latag              // likelihood weight for tagging data
//  init_number sw6
  !! cout << "Data Section Completed" << endl;
  !! cout << "lg " << lg << endl;

 
//===============================================================================================================
// 2.0  Define Parameters 
//===============================================================================================================
PARAMETER_SECTION
  init_bounded_vector log_q(1,2,-20.5,20.0,1)     //Commercial firherey Catchability 
  init_bounded_number log_qw(-10.5,20.0,-1)    // Winter sruvery catchability  
  init_bounded_number log_initpop(2.0,15.0,1) //Initial population size  
  init_bounded_number log_recscale(2.0,12.0,1)
  init_bounded_dev_vector log_relrec(fyear,lyear-1,-40.0,40.0,1)
  init_bounded_vector flnp(1,na-1,-10.0,10.0,1)     //Initial year length proportion (log scale)
  init_bounded_number r1(0.1,0.9,1)             //proportion of recruits to length group 1
  init_bounded_number log_mo1(-10.5,-1.0,1)     //parameter for molting probability
  init_bounded_vector log_st1(1,1,-15.0,1.0,1)      //parameter for selectivity: trawl NOAA & ADF&G survey 
  init_bounded_number log_sw1(-15.0,1.0,1)      //parameter for selectivity: winter project
  init_bounded_number sw3(0.1,1.0001,1)           //selectivity for length group NA: winter project
  init_bounded_vector log_sc1(1,1,-15.0,1.0,1)      //parameter for selectivity: summer fishery before 1993, after 1993
  init_bounded_number advar(0.0,6.0,1)         //selectivity of old-shell crabs in winter project
  init_bounded_number qtno(0.1,1.0,1)           // Catchability of NOAA trawl survey
  init_bounded_number M(0.02,0.34,M_phase)           //Mortality  Does not use
  init_bounded_number ms(1.0,5.0,ms_phase)         //Mortality multiplier for the last length group
//  init_bounded_number p4(1.0,1.0,-1)             //proporion of marketable legal for length 4 class since 2005
  init_bounded_number sigma(0.,30.,2)            // growth matrix standard deviation       
  init_bounded_vector ig(1,2,0.,20.,1)       //  growth matrix mean increment  
  
// Delived parameters
  vector st(1,nyt)       //annual sample size from trawl survey
  vector sw(1,nyw)       //annual sample size from winter project    
  vector sc(fyear,lyear)       //annual sample size from summer comfishery  
  vector so(1,nyo)       //annual sample size from observer's data  
  matrix ont(1,na,1,nyt) //proportion of length group for new-shell, trawl survey
  matrix oot(1,na,1,nyt) //proportion of length group for old-shell, trawl survey  
  matrix onw(1,na,1,nyw) //proportion of length group for new-shell, winter survey
  matrix oow(1,na,1,nyw) //proportion of length group for old-shell, winter survey
  matrix onc(1,na,fyear,lyear)  //proportion of length group for new-shell, summer fishery
  matrix ooc(1,na,fyear,lyear)  //proportion of length group for old-shell, summer fishery
  matrix ono(1,na,1,nyo) //proportion of length group for new-shell, observer survey
  matrix ooo(1,na,1,nyo) //proportion of length group for old-shell, observer survey
  vector log_rec(fyear,lyear-1)                        //Annual recruitment log scale 
  matrix nps(1,na,fyear,lyear+1)                   //new-shell population abundance: summer
  matrix ops(1,na,fyear,lyear+1)                   //old-shell population abundance: summer
  matrix npw(1,na,fyear,lyear+1)                   //new-shell population abundance: winter
  matrix opw(1,na,fyear,lyear+1)                   //old-shell population abundance: winter
  matrix enc(1,na,fyear,lyear)                     //Modeled new-shell length proportion: Summer commercial catch
  matrix eoc(1,na,fyear,lyear)                     //Modeled old-shell length proportion: Summer commercial catch  
  vector ecpue(fyear,lyear)                        //Modeled estimated cpue: summer fishery
  vector cvcpue(fyear,lyear)                       //cv of cppue   
  matrix enwp(1,na,fyear,lyear)                      //Modeled new-shell length proportion: Winter subsistence catch
  matrix eowp(1,na,fyear,lyear)                      //Modeled old-shell length proportion: Winter subsistence catch
  matrix enwpd(1,na,fyear,lyear)                     //Modeled new-shell length proportion: Winter subsistence discards
  matrix eowpd(1,na,fyear,lyear)                     //Modeled old-shell length proportion: Winter subsistence discards
  matrix enwc(1,na,fyear,lyear)                      //Modeled new-shell length proportion: Winter commercial catch  
  matrix eowc(1,na,fyear,lyear)                      //Modeled old-shell length proportion: Winter commercial catch
  vector twsd(fyear,lyear)                             //Winter subsitence catch discarded 
  vector tb(fyear,lyear)                               //estimated mean annual summer exploitable abundance
  vector tw(fyear,lyear)                           	//estimated mean annual winter exploitable abundance
  vector twp(fyear,lyear)                           	//estimated mean annual winter exploitable abundance for pot survey 
  vector rec(fyear,lyear)                        //annual recruits: starting in year 1.
  matrix ent(1,na,1,nyt)                        //new-shell length propotion: summer trawl survey
  matrix eot(1,na,1,nyt)                        //old-shell length propotion: summer trawl survey  
  matrix ent0(1,na,1,nyt)                       //Modeled trawl survey newshell abundance
  matrix eot0(1,na,1,nyt)                       //Modeled trawl survey oldshell abundance  
  matrix enw(1,na,1,nyw)                        //new-shell length proportion: Winter pot survey
  matrix eow(1,na,1,nyw)                        //old-shell length proportion: Winter pot survey  
  matrix eno(1,na,1,nyo)                        //new-shell length proportion: Summer commercial observer discards
  matrix eoo(1,na,1,nyo)                        //old-shell length proportion: Summer commercial observer discards
  vector ettq(1,nyt)                            //estimated catchability adjusted total annual abundance: trawl survey
  vector ett(1,nyt)                             //estimated total annual abundance: trawl survey
  vector ett1(1,nyt)                            //estimated annual non-legal crab abundance: trawl survey
  vector ett2(1,nyt)                            //estimated annual legal crab abundance: trawl survey
  vector ewcpue(1,nyw)                          // estimated cpue for pot survey
  vector mp1(1,na)                              //molting probability
  vector mp0(1,na)                              //working vector: molting probability
//Commerical Fishery Selectivity Valiables
  matrix selc1(1,2,1,na)                        //selectivity of summer fishery 
  matrix sel(1,na,fyear,lyear)                  //annual selectivity of summer fishery
  vector st1(1,1)
  vector sc1(1,1)
//Trawl Survey Selectivity Valiables  
  matrix selt1(1,2,1,na)                            //selectivity of trawl NOAA ADDFGsurvey
  matrix selt(1,na,fyear,lyear)                 //selectivity of combined trawl survey
  vector selw(1,na)                             //working variable: selectivity of winter project.
  vector Mn(1,na)                               //working variable: natural mortality for new-shell crabs.
  vector q(1,2)                                     //catchability: summer fishery 1977-1992
  number qw
  vector FSum(fyear,lyear)                       // Exploitation rate (Summer fishery)
  vector FWin(fyear,lyear)                       // Exploitation rate (Winter fishery)
  vector FSub(fyear,lyear)                       // Exploitation rate (Subsistance fishery)
  vector bc(fyear,lyear)                         // Summer commercial diccards.
  vector bcw(fyear,lyear)                        // winter commercial discards.
  vector T_var(fyear,lyear);                     // additional variance for cppue 
  vector expn(1,na-1);
  vector expo(1,na-2);
  vector npp(1,na);
  vector opp(1,na-1);
  vector tt1(1,nyt)            //observed annual non-legal crab abundance: trawl survey
  vector tt2(1,nyt)            //observed annual legal crab abundance: trawl survey
  vector tf(1,13);                              // Likelihood vector
  
//==========================================================================
//Tagging related parameters   
//==========================================================================  
  matrix tagrecap1(1,na,1,na)   //tagging freq before 1993 for Year 1 
  matrix tagrecap2(1,na,1,na)   //tagging freq before 1993 for Year 2
  matrix tagrecap3(1,na,1,na)   //tagging freq before 1993 for Year 3
  matrix tagrecap12(1,na,1,na)   //tagging freq after 1993 Year 1
  matrix tagrecap22(1,na,1,na)   //tagging freq after 1993 Year 2
  matrix tagrecap32(1,na,1,na)   //tagging freq after 1993 Year 3
  matrix ptagrecap1(1,na,1,na)   //prob tagging before 1993 Year 1 
  matrix ptagrecap2(1,na,1,na)   //prob tagging before 1993 Year 2
  matrix ptagrecap3(1,na,1,na)   //prob tagging before 1993 Year 3
  matrix ptagrecap12(1,na,1,na)   //prob tagging after 1993 Year 1
  matrix ptagrecap22(1,na,1,na)   //prob tagging after 1993 Year 2
  matrix ptagrecap32(1,na,1,na)   //prob tagging after 1993for Year 3
  matrix tgr(1,na,1,na)          // True transition growth matrix
  matrix mgr1(1,na,1,na)          //molting probability adjsteted matrix Year 1
  matrix mgr2(1,na,1,na)          //molting probability adjsteted matrix Year 2
  matrix mgr3(1,na,1,na)          //molting probability adjsteted matrix Year 3
  matrix egr1(1,na,1,na)    // estimated growth matrix for sel1 Year 1
  matrix egr2(1,na,1,na)   //estimated growth matrix for sel1 Year 2
  matrix egr3(1,na,1,na)   //estimated growth matrix for sel1 Year 3
  matrix egr12(1,na,1,na)   //estimated growth matrix for sel2 Year 1
  matrix egr22(1,na,1,na)   //estimated growth matrix for sel2 Year 2
  matrix egr32(1,na,1,na)   //estimated growth matrix for sel2 Year 3 
  vector ef(1,6)            // objective function vector for tagging data
  
//===============================================================================================================

  objective_function_value f
  sdreport_vector last_y(1,na)
  sdreport_vector legaln(fyear,lyear)
  sdreport_vector legalb(fyear,lyear)
  sdreport_vector mmb0(fyear,lyear+1)
  vector legal(fyear,lyear+1)
  vector mmb(fyear,lyear+1)
  sdreport_number last_legal  
  number last_subl
  number bmsy
  sdreport_number last_mmb
  likeprof_number last_ofl
  !! cout << "Parameter Section Completed" << endl;

//===============================================================================================================
INITIALIZATION_SECTION
  log_q         -6.5
  r1             0.6
  log_mo1       -2.41358
  log_st1         -2.5
  log_sc1         -2.31267
  log_sw1         -2.31267
  log_qw         -6.5
//  p4 1.0
  advar 0.5
  qtno 1.0
  sigma  5.0
  ig 10.0
//===============================================================================================================

//===============================================================================================================
// 4.0  Preliminary Calculation
//===============================================================================================================
PRELIMINARY_CALCS_SECTION
  int i,j;
  dvariable tt0,n0; // Calculated working variable does not need defferentiated
//  cout << "Starting preliminary calcs" << endl;
  M = M2;
  ms = ms6;
  Mn = M;
  Mn(na) = ms*M;
  tt0.initialize();
  n0.initialize();

//===============================================================================================================
// 4.1  Calculate annual sample size and length proportions     
//===============================================================================================================
  st = colsum(nont)+colsum(noot);       //annual sample size from trawl survey
  sw = colsum(nonw)+colsum(noow);       //annual sample size from winter pot survey
  sc = colsum(nonc)+colsum(nooc);       //annual sample size from commercial fishery survey
  so = colsum(nono)+colsum(nooo);       //annual sample size from oserver survey

  for (i=1;i<=na;i++)
    {
     ont(i) = elem_div(nont(i),st);
     oot(i) = elem_div(noot(i),st);
     onw(i) = elem_div(nonw(i),sw);
     oow(i) = elem_div(noow(i),sw);
     onc(i) = elem_div(nonc(i),sc);
     ooc(i) = elem_div(nooc(i),sc);
	 onc(i,fyear) = 0;  // no fishery in 1976
	 ooc(i,fyear) = 0;  // no fishery in 1976
	 onc(i,1991) = 0;   // no fishery in 1991
	 ooc(i,1991) = 0;	// no fishery in 1991 
     ono(i) = elem_div(nono(i),so);
     ooo(i) = elem_div(nooo(i),so); 
	}
  	
//===============================================================================================================
// 4.2  Set Maximum effective sample size: 
//      Trawl & Summer pots:  50% of actual sample size or maxss
//      Winter pots, Summer commercial, Observer:  10% of actual sample size or maxsc
//===============================================================================================================

// Calculate maximum effective smple sieze for Trawl survey
  for (i=1;i<=nyt;i++)
   {
    st(i) *= efn(1);
    if (st(i) > maxs(1)) st(i)=maxs(1);
   }
 cout << st << endl;
// Calculate maximum effective smple sieze for Winter pot survey
  for (i=1;i<=nyw;i++)
   {
    sw(i) *= efn(2);
    if (sw(i) > maxs(2)) sw(i)=maxs(2);
   }
// Calculate maximum effective smple sieze for Summer commercial 
  for (i=fyear;i<=lyear;i++)
   {
    sc(i) *= efn(2);
    if (sc(i) > maxs(2)) sc(i)=maxs(2);
   }
// Calculate maximum effective smple sieze for Observer
  for (i=1;i<=nyo;i++)
   {
    so(i) *= efn(2);
    if (so(i) > maxs(2)) so(i)=maxs(2);
   }

//===============================================================================================================
// 4.3  Calculate cv of cpue index      
//===============================================================================================================
  cvcpue = elem_div(secpue+1.e-3,stcpue+1.e-3);

//===============================================================================================================
// 4.4  Calculate Winter subsistence discards (twsd)
//      For early years, total number of crab caught is not available. 
//      Estimae discards based on average proporion of discards   
//===============================================================================================================
// Discards twsd is a subtraction of retained catach (tws) from total catch (twsd) 
  twsd = twst - tws;
// Calculate Average discards to retained ratio
  for (i=fyear;i<=lyear;i++)
  {
  if (twsd(i) > 0)
        {
          tt0 += twsd(i)/tws(i);  // Sum proprotion of discards
          n0 += 1;                // Number of sample 
        }
    }
// For early unknown discards, estimate by average discards rate. 		  
  for (i=fyear;i<=lyear;i++)
    {
    if (twsd(i)< 0)
        {
          twsd(i) = tws(i)*(tt0/n0);
        }   
    }
 cout << twsd << endl;
   
//===============================================================================================================
// 4.6  Create tagrecovery matrices       
//===============================================================================================================
    tagrecap1.initialize();
    tagrecap2.initialize();
    tagrecap3.initialize();
    tagrecap12.initialize();
    tagrecap22.initialize();
    tagrecap32.initialize();  
    ptagrecap1.initialize();
    ptagrecap2.initialize();
    ptagrecap3.initialize();
    ptagrecap12.initialize();
    ptagrecap22.initialize();
    ptagrecap32.initialize();           
          
// Create tag-recoverey frequency matrix by recovery years          
  for(i=1;i<=21;i++)
  {
        tagrecap1(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,3);
        tagrecap2(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,4);
        tagrecap3(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,5);
        tagrecap12(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,3);
        tagrecap22(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,4);       
        tagrecap32(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,5);  
  }
  
// Create tag-recoverey probability matrix by recovery years   
  
  for(i=1;i<=na;i++)
  {  
  ptagrecap1(i)=tagrecap1(i)/rowsum(tagrecap1)(i);
  ptagrecap2(i)=tagrecap2(i)/rowsum(tagrecap2)(i);
  ptagrecap3(i)=tagrecap3(i)/rowsum(tagrecap3)(i);
  ptagrecap12(i)=tagrecap12(i)/rowsum(tagrecap12)(i);
  ptagrecap22(i)=tagrecap22(i)/rowsum(tagrecap22)(i);  
  ptagrecap32(i)=tagrecap32(i)/rowsum(tagrecap32)(i);
  }
// Change to zero with no tag recovery data  
  ptagrecap1(1) = 0.0;
  ptagrecap12(1) = 0.0;  
  ptagrecap22(6) = 0.0;
  ptagrecap32(5) = 0.0;  
  ptagrecap32(6) = 0.0;    
 
 cout << "End preliminary calcs" << endl; 

 
//===============================================================================================================
RUNTIME_SECTION
  convergence_criteria 1E-6
  maximum_function_evaluations 10000
//===============================================================================================================

//===============================================================================================================
// 5.0  Procedure       
//===============================================================================================================
PROCEDURE_SECTION

  convert_parameters_into_rates();
//  cout <<"OK for convert_parameters..."<<endl;

  get_first_year_abundance();
//  cout <<"OK for get_first_year..."<<endl;

  growth_matrix();
//  cout <<"OK growth matrix..."<<endl;
  
  get_number_by_size();
//  cout <<"OK for get_number_by_size..."<<endl;

  get_proportion_and_effort();
//  cout <<"OK for get_proportion_and..."<<endl;

  evaluate_the_objective_function();
  
   
//===============================================================================================================
// 6.0  Function:  Convert parameters into rates      
//===============================================================================================================
FUNCTION convert_parameters_into_rates
  int i, j;
// working variables: molting probability
  dvariable mo1,sw1;    
// working variables: Selectivity 
//  dvariable sc1,sc3,sc5, st1,st3,sw1;   
  double pp, L,Lm;                                       
//catchability coefficient
   q = mfexp(log_q);
   qw = mfexp(log_qw);   
//molting probability and selectivity for summer commercial fisheries
   mo1 = mfexp(log_mo1);
// Summer Commercial fisheris selectivity parameters    
   sc1 = mfexp(log_sc1);
// Summer Trawl Survey selectivity parameters       
   st1 = mfexp(log_st1); 
   sw1 = mfexp(log_sw1);
   
//===============================================================================================================
// 6.1  Set molting and net selectivity funuctions      
//===============================================================================================================
// Lm is the maximum mean length of crab
    Lm = slm+slt*(na-1);
 // calculate molting and selectivity functions   
   for (j=1;j<=na;j++)
    {
     L = (double(j)-1.0)*slt;
     L += slm;
     mp1(j) = 1.0-1.0/(1.0+mfexp(mo1*(slm-L)+log(1.0/0.001-1.0)));  // molting probability : reverse logistic function//
     selc1(1,j) = 1.0/(1.0+mfexp(sc1(1)*(Lm-L)+log(1.0/0.999-1.0)));  // commercial pot selectivity 1976-1992: logistic function//
     selc1(2,j) = 1.0/(1.0+mfexp(sc1(1)*(Lm-L)+log(1.0/0.999-1.0)));  // commercial pot selectivity 1993 - : logistic function//
//     selw(j) = 1.0/(1.0+mfexp(sw1*(Lm-L-slt)+log(1.0/0.999-1.0)));   // winter pot selectivity 1981-2011: logistic function//
     selw(j) = 1.0-1.0/(1.0+mfexp(sw1*(slm-L)+log(1.0/0.001-1.0)));  // winter pot selectivity 1981-2011: reverse logistic function//
     selt1(1,j) = 1.0/(1.0+mfexp(st1(1)*(Lm-L)+log(1.0/0.999-1.0)));   // NOAA summer trawl net selectivity 1976-1992: logistic function//
     selt1(2,j) = 1.0/(1.0+mfexp(st1(1)*(Lm-L)+log(1.0/0.999-1.0)));   // ADFG summer trawl net selectivity: 1996- logistic function//     
    }

   selw(1) = sw3;   
   
//===============================================================================================================
// 6.2  Adjust molting and selectivity functions      
//===============================================================================================================
  
// Assign Trawl survey seletivity function in different time periods. 
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(1))    // period 1976 - 1992 
      for (j=1; j<=na; j++) selt(j,i) = selt1(1,j);   //NOAA trawl survey 
      else                          
      for (j=1; j<=na; j++) selt(j,i) = selt1(1,j);   //ADFG trawl survey
    }


// Assign Summer commrecial catch seletivity function in different time periods. 
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(1))    // period 1976 - 1992 
      for (j=1; j<=na; j++) sel(j,i) = selc1(1,j); 
      else                          // period 2008 - present
      for (j=1; j<=na; j++) sel(j,i) = selc1(2,j);
    }

//===============================================================================================================
// 7.0  Function:  Estimate the First Year (Feb 1st) abundance and length composition
//      Both abundance and length composition is estimated     
//===============================================================================================================
FUNCTION get_first_year_abundance
  int i, j;
  dvariable first_y;

// Assign first year abundace: 
  first_y = mfexp(log_initpop);
  
// Estimate first year length proportion    
   for (j=1;j<=(na-1);j++) expn(j) = mfexp(flnp(j));
   npp(1,na-1) = expn/(1+sum(expn));
   npp(na) = 1-sum(npp(1,na-1));
   
// Assign first year abundance on Feb 1st.  Assume no oldshell crabs     
   for (j=1;j<=na;j++) 
	{
		npw(j,fyear) = npp(j)*first_y;
		opw(j,fyear) = 0.0;
	}
  
// Assign Recruits   
  log_rec = log_relrec+log_recscale;  

  for (i=fyear;i<=(lyear-1);i++) rec(i) = mfexp(log_rec(i));  
// recruits of the last year are assumed to be average of the most recent five years
   rec(lyear) = 0.2*(rec(lyear-1)+rec(lyear-2)+rec(lyear-3)+rec(lyear-4)+rec(lyear-5));
   
//=========================================================================================
// 8.0  Function:  Gorwth matrix
//     Estimate Growth transition matrix from tagging data   
//========================================================================================= 
  
FUNCTION growth_matrix
   int i, j; 
   dvariable t1, mu, fa, fb;
  tgr.initialize();
  mgr1.initialize();
  egr1.initialize();
  egr2.initialize();
  egr3.initialize();
  egr12.initialize();
  egr22.initialize();
  egr32.initialize();
  ef.initialize();
// assign normal probability for each length class   
   for (i=1;i<=(na-1);i++)
   {
    mu = slm+ig(1)+ig(2)*i;
// Assume that crab does not shrink    
    for (j=i;j<=na;j++)
     {
     fa = slm+slt*(j-1.5);
     fb = slm+slt*(j-0.5);
     tgr(i,j) = cumd_norm((fb-mu)/sigma)-cumd_norm((fa-mu)/sigma);
     }
   }
   tgr(na,na) = 1;
   
// Normalized the proabaiblity  
 
   for (i=1;i<=na;i++)
   {
    tgr(i) = tgr(i)/rowsum(tgr)(i);
   }
   // growth increment of the last length class is 1.0
 
// Include molting probability 
 
   for (i=1;i<=na;i++)
   {
    mgr1(i) = tgr(i)*mp1(i);
   }

// Add (1-mi) to the model :  mgr is the moling probability adjusted growth-matrix
 
   for (i=1;i<=na;i++)
   {
    mgr1(i,i) += (1-mp1(i));
   }
    
 // Calculate expected matrix for year 2-3
   mgr2 = mgr1*mgr1;    //estimated growth matrix Year 2
   mgr3 = mgr2*mgr1;   //estimated growth matrix Year 3
// 

//============================================================================================================
//    Adjustment for selectivity
//============================================================================================================
  
// multiply by fishery selectivity
   for (i=1;i<=na;i++)
    {
    egr1(i) = elem_prod(mgr1(i),selc1(1));
    egr2(i) = elem_prod(mgr2(i),selc1(1));
    egr3(i) = elem_prod(mgr3(i),selc1(1));
// Selectivity for 1993-present
    egr12(i) = elem_prod(mgr1(i),selc1(2));
    egr22(i) = elem_prod(mgr2(i),selc1(2));
    egr32(i) = elem_prod(mgr3(i),selc1(2));    
    }

// Normalize  
   for (i=1;i<=na;i++)
    {
    egr1(i) = egr1(i)/rowsum(egr1)(i);
    egr2(i) = egr2(i)/rowsum(egr2)(i);
    egr3(i) = egr3(i)/rowsum(egr3)(i);
    egr12(i) = egr12(i)/rowsum(egr12)(i);
    egr22(i) = egr22(i)/rowsum(egr22)(i);
    egr32(i) = egr32(i)/rowsum(egr32)(i);   
    }
//  cout << egr1 << endl;
//  cout << egr12 << endl;
  
 

//===============================================================================================================
// 8.0  Calculate Population dynamics by length class  
//===============================================================================================================
FUNCTION get_number_by_size
  int i,j,k;
  dvariable pp, pp1, pp2, tt1, ttt0, ttt1;
  dvar_vector tt0(1,na);    // Abundance after summer fishery
  dvariable TotalCatch, CatNum;    
  enc.initialize();
  eoc.initialize();
  eno.initialize();
  eoo.initialize();
  enwp.initialize();
  eowp.initialize();  
  enwpd.initialize();
  eowpd.initialize();
  enwc.initialize();
  eowc.initialize();  
  tt0.initialize();

// Estimate projected crab abundance 
  for (i=fyear;i<=lyear;i++)
  {
// Assign molting probability   
      mp0 = mp1;
 	
//===============================================================================================================
// 8.1 Winter commercial catch and discards size composition 
//===============================================================================================================
//propotions for winter fishery and subsistence fishery
	tt1 = 0.0;  //commercial catch by length class 
	tw(i) = 0.0;
	twp(i) = 0.0;
// Calculate population selective to winter commercial fisheries 
    for (j=1;j<=na;j++)
     {
      tw(i) += (npw(j,i)+opw(j,i))*selw(j)*lg(j); // For commercial Catch 	  
     }

// Calculate population selective to winter pot survey
    for (j=1;j<=na;j++)
     {
      twp(i) += (npw(j,i)+opw(j,i))*selw(j); // For commercial Catch 	  
     }
          
// Calculate exploitable abundance at the middle of the season
     twp(i) = twp(i) - 0.5*(twc(i)+tws(i));  

// Calculate proportion of new shell and old shell in winter commercial catch   	 
    for (j=1;j<=na;j++)
     {
      enwc(j,i) = npw(j,i)*selw(j)*lg(j)/tw(i);
      eowc(j,i) = opw(j,i)*selw(j)*lg(j)/tw(i);
     }
// Calculate Legal abundance in winter 	 
	for (j=1;j<=na;j++)
     {
      tt1 += (npw(j,i)+opw(j,i))*lg(j); 	  
     }
	 
//===============================================================================================================
// 8.2 Winter subsistence catch and discards size composition 
//     Assume that Subsistence fishers took all length >= 3 crab  and discarded length 1 & 2 crab
//===============================================================================================================

// Calculate adult population cathable to winter sfisheries 
	pp1 =0.0;
    for (j=3;j<=na;j++)
     {
      pp1 += (npw(j,i)+opw(j,i))*selw(j);  // For subsistence catch   
     }

// Calculate proportiojn of new shell and old shell in winter subsistence catch   	 
    for (j=3;j<=na;j++)
     {
	  enwp(j,i) = npw(j,i)*selw(j)/pp1;
      eowp(j,i) = opw(j,i)*selw(j)/pp1;
     }

// Calculate juvenile population cathable to winter sfisheries 
	pp2 = 0.0;
    for (j=1;j<=2;j++)
     {
      pp2 += (npw(j,i)+opw(j,i))*selw(j);  // For subsistence catch   
     }
	 
// Calculate proportiojn of new shell and old shell discarded.  	 
    for (j=1;j<=2;j++)
     {
	  enwpd(j,i) = npw(j,i)*selw(j)/pp2;
      eowpd(j,i) = opw(j,i)*selw(j)/pp2;	  
	  }

    // Extract exploitation rate
    FSum(i) = tc(i)/tb(i);
    FSub(i) = tws(i)/pp1;
    FWin(i) = twc(i)/tt1;

//===============================================================================================================
// 8.3 Calculate abundance for summber July 1st:
//     Summer population is a survivor of winter fisheries
//===============================================================================================================
  bcw(i) = 0.0;
    for (j=1;j<=na;j++)
	{
     nps(j,i) = (npw(j,i)-enwp(j,i)*tws(i)-enwc(j,i)*twc(i)-(twc(i)/tt1)*(npw(j,i))*selw(j)*(1-lg(j))*hm(2)-twsd(i)*enwpd(j,i)*hm(2))*exp(-0.417*Mn(j));  
     ops(j,i) = (opw(j,i)-eowp(j,i)*tws(i)-eowc(j,i)*twc(i)-(twc(i)/tt1)*(opw(j,i))*selw(j)*(1-lg(j))*hm(2)-twsd(i)*eowpd(j,i)*hm(2))*exp(-0.417*Mn(j));	
	 if (nps(j,i) < 0.0) nps(j,i) = 0.001; // summer abundance should not go bellow zero
	 if (ops(j,i) < 0.0) ops(j,i) = 0.001; // summer abundance should not go bellow zero	
     bcw(i) += (twc(i)/tt1)*(npw(j,i)+opw(j,i))*selw(j)*(1-lg(j))*hm(2); 
	}
  
//===============================================================================================================
// 8.4  Calculate Summber Crab length proportion catchable to commercial fishery 
//===============================================================================================================
 
	tb(i) = 0.0; 
    for (j=1;j<=na;j++)
	{
// 	    if (j==4&i>= scp) // Since 2005 Commerical buyers accest only CW > 5.0 crab instead of CW > 4.5	 
//		{
//	    tb(i) += (nps(j,i)+ops(j,i))*sel(j,i)*lg(j)*p4;  
//		}
//	    else
//	    {
        tb(i) += (nps(j,i)+ops(j,i))*sel(j,i)*lg(j); //total summer crab abundance available to 
//	    }
	}		
	// mean exploitable leagal abundance does not go negative
    if (tb(i) < 0.001) tb(i) = 0.001;
	// Calculate proprotion of newshell and oldshell by comm fiish 	
	for (j=1;j<=na;j++)
     {
        enc(j,fyear) = 0.0;
        eoc(j,fyear) = 0.0;
// 	    if (j==4&i>= scp) // Since 2005 Commerical buyers accest only CW > 5.0 crab instead of CW > 4.5	 
//		{
//		enc(j,i) = nps(j,i)*sel(j,i)*lg(j)*p4/tb(i);
//		eoc(j,i) = ops(j,i)*sel(j,i)*lg(j)*p4/tb(i);
//		}	 
//	     else
//		{
	 	enc(j,i) = nps(j,i)*sel(j,i)*lg(j)/tb(i);
		eoc(j,i) = ops(j,i)*sel(j,i)*lg(j)/tb(i);
//		}
	}

//===============================================================================================================
// 8.5  Calculate Crab aundance after summer fishery and discards 
//===============================================================================================================
    ttt0 = 0.0;   // total number of legal crab on July 1st
    bc(i) = 0.0;  // discards biomass
	
// Calculate total number of legal crab on July 1st (ttt0)
    for (k=1;k<=na;k++) {ttt0 += (nps(k,i)+ops(k,i))*lg(k);}	
// if projected abundance is lower than actual catch then projected abundance should be adjusted to 1.2 times of actual harvest 
    if (ttt0 < tc(i)) {ttt0 = tc(i)*1.2;}
// if projected abundance is lower than 0.001 thaen  projected cach should be adjusted to 0.001
    else if (ttt0 < 0.001) ttt0 = 0.001;
	
// Crab abundance right after summer commercial fishery (tt0) and bycatch (bc)
    for (k=1;k<=na;k++)
     {	  
      tt0(k) = (nps(k,i)+ops(k,i))*exp(-ys(i)*Mn(k))-tc(i)*(enc(k,i)+eoc(k,i))-(tc(i)/ttt0)*(nps(k,i)+ops(k,i))*sel(k,i)*(1.0-lg(k))*hm(1);	  
      if (tt0(k) < 0.0) tt0(k)= 0.001;  // Stock gap measure: abundance of each length class should not go below zero
      bc(i) += (tc(i)/ttt0)*(nps(k,i)+ops(k,i))*sel(k,i)*(1.0-lg(k))*hm(1)*wm(k);  // Bycatch biomass 
// Commercial accept only > CW 5 inch since 2005 and sub comercial size crabs are assumed discarded
//	if (i >= scp) bc(i) += (tc(i)/ttt0)*(nps(4,i)+ops(4,i))*sel(4,i)*(lg(k)*(1-p4))*hm(1)*wm(4);
	  }
     
//===============================================================================================================
// 8.6  Calculate Crab abundance on Feb 1st 
//===============================================================================================================
     
// Calculate New Shell popululation abundance by length class:
    for (j=1;j<=na;j++)
     {
      pp = 0.0; 
      for (k=1;k<=j;k++) pp += tgr(k,j)*tt0(k)*mp0(k); //Each crab molts right after fishery ended      
      npw(j,i+1) = pp*exp(-(0.583-ys(i))*Mn(j));  //New shell crab on Feb 1st 
// Add Recruits to length 1 and 2
      if (j==1) npw(j,i+1) += r1*rec(i);
      if (j==2) npw(j,i+1) += (1.0-r1)*rec(i);
      opw(j,i+1) = tt0(j)*(1.0-mp0(j))*exp(-(0.583-ys(i))*Mn(j));   //Old shell crab are unmolted crab  
     }
   }

  
//===============================================================================================================
// 8.7  Calculate February 1st legal Crab and mature male biomass  for assessment
//===============================================================================================================

 for (i=fyear; i<=lyear; i++)
   {
    legal(i) = 0.0; legalb(i) = 0.0;
    for (j=1; j<=na; j++)
    {
       legal(i) += (npw(j,i)+opw(j,i))*lg(j);
       legalb(i) += (npw(j,i)+opw(j,i))*lg(j)*wm(j);
    }
    mmb(i) = 0.0;
    
    for (j=3; j<=na; j++) mmb(i) += (npw(j,i)+opw(j,i))*wm(j);
    mmb0(i) = mmb(i);    
    legaln(i) = legal(i);
	}
   
  last_y = column(npw,lyear+1);
  for (j=1; j<=na; j++) last_y(j) += opw(j,lyear+1);

// Calculate Jul 1st projected abundance (Assume no winter fishery)
  last_legal = 0.0;
  for (j=1; j<=na; j++) last_legal += last_y(j)*lg(j)*wm(j)*selc1(2,j);
  
  last_subl = 0.0;
  for (j=1; j<=na; j++) last_subl += last_y(j)*(1-lg(j))*wm(j)*selc1(2,j);

  last_mmb = 0.0;
  for (j=3; j<=na; j++) last_mmb += (npw(j,lyear+1)+opw(j,lyear+1))*wm(j);
  mmb(lyear+1) = last_mmb;
  mmb0(lyear+1) = last_mmb;
  bmsy = (sum(mmb))/(lyear-fyear+1);
  last_ofl = last_legal*(1.0-exp(-M*(last_mmb/bmsy-0.1)/0.9));

  
//===============================================================================================================
// 9.0  produce model estimated observation data  
//===============================================================================================================
FUNCTION get_proportion_and_effort
  int i,j;
  dvariable bf,af,pp;
  ett.initialize();
  ett1.initialize();
  eow.initialize();
  eoo.initialize();
  ecpue.initialize();


//===============================================================================================================
// 9.1  Predicted Trawl survey abundance 
//      Survey abundance is adjusted by fishery and mortality 
//      When mid-point of survey date is later than that of commercial fishery. 
//         Remove commercial catch from survey abundance. 
//===============================================================================================================
  for (i=1;i<=nyt;i++)
   {
    if (yt(i) > ys(it(i))) // Mid-point of survey date is later than that of commercial fishery
     {
      bf = ys(it(i));              //time lag from July 1 to fishery
      af = yt(i) - ys(it(i));     //time lag from fishery to survey
     }
    else
     {
      bf = yt(i);
      af = 0.0;
     }
    for (j=1;j<=na;j++)
     {
      ent0(j,i) = ((nps(j,it(i)))*exp(-bf*Mn(j))-(enc(j,it(i)))*pct(i)*tc(it(i)))*exp(-af*Mn(j))*selt(j,it(i));
      eot0(j,i) = ((ops(j,it(i)))*exp(-bf*Mn(j))-(eoc(j,it(i)))*pct(i)*tc(it(i)))*exp(-af*Mn(j))*selt(j,it(i)); 
	  if (ent0(j,i) < 0.0) ent0(j,i) = 0.0;
	  if (eot0(j,i) < 0.0) eot0(j,i) = 0.0;	  
      ett(i) += ent0(j,i)+eot0(j,i);
      if (j<4) ett1(i) += ent0(j,i)+eot0(j,i);
     }
    if (ett(i) <= 0.0) ett(i) = 0.00001;
        ett2(i) = ett(i) - ett1(i);
    for (j=1;j<=na;j++)
     {
      ent(j,i) = ent0(j,i)/ett(i);
      eot(j,i) = eot0(j,i)/ett(i);	  
      if (ent(j,i) < 0.0) ent(j,i) = 0.0;
	  if (eot(j,i) < 0.0) eot(j,i) = 0.0;
     }

// Estimate Survey Q for NOAA trawls
  if(it(i) < scy(1)) {ettq(i) = qtno*ett(i);}
    else
    ettq(i) = (ett(i));
   }

//===============================================================================================================
// 9.2  Predicted Winter Pot survey length proportion 
//===============================================================================================================
  for (i=1;i<=nyw;i++)
   {
    pp = 0.0;
    for (j=1; j<=na; j++)
     pp += (npw(j,iw(i))+opw(j,iw(i)))*selw(j);
    if (pp <= 0.0) pp = 0.000001;
    for (j=1;j<=na;j++)
     {
      enw(j,i) = npw(j,iw(i))*selw(j)/pp;
      eow(j,i) = opw(j,iw(i))*selw(j)/pp;
     }
// Calculate winter survey cpue     
      ewcpue(i) = qw*twp(iw(i)); 	 
   }

//===============================================================================================================
// 9.3  Predicted Observer survey length proportion  
//===============================================================================================================
  for (i=1;i<=nyo;i++)
   {
    pp = 0.0;
    for (j=1;j<=na;j++)
     pp += (nps(j,io(i))+ops(j,io(i)))*sel(j,io(i))*(1.0-lg(j));
    if (pp <= 0.0) pp = 0.000001;
    for (j=1;j<=na;j++)
     {
      eno(j,i) = nps(j,io(i))*sel(j,io(i))*(1.0-lg(j))/pp;
      eoo(j,i) = ops(j,io(i))*sel(j,io(i))*(1.0-lg(j))/pp;
     }
   }

//===============================================================================================================
// 9.4  Predicted CPUE
//===============================================================================================================

  for (i=fyear;i<=lyear;i++)
  {
    tb(i) = tb(i) - 0.5*tc(i);  //exploitable abundance at the middle of the season
    if (tb(i) < 0.001) tb(i) = 0.001;
    if (i <= scy(1))
	 {
	 ecpue(i) = q(1)*tb(i);
	 }
    else
	 {	
	 ecpue(i) = q(2)*tb(i);
	 }
  }
  
//fishing cpue: remove data 
  for (i=fyear;i<=lyear;i++)
   {
   if (stcpue(i) <= 0.0) {ecpue(i) = 0.0;}
   }
   

//===============================================================================================================
// 10.0  Likelihood Calculation  
//===============================================================================================================
FUNCTION evaluate_the_objective_function

// Log likelihood for trawl survey 
  tf(1) = 0.5*norm2(elem_div((log(tt+1.e-3)-log(ettq+1.e-3)),sqrt(log(elem_prod(cv,cv)+1.0))));

// Log likelihood standard cpue in winter fishery 
//s  tf(2) = lamw*norm2(log(wcpue+1.e-3)-log(ewcpue+1.e-3))/(2*SDW*SDW); 
   
// Log likelihood standard cpue in summer fishery 
  T_var = sqrt(log(elem_prod(cvcpue,cvcpue)+1.0)+advar); 
  tf(3) = lamc*(sum(log(T_var))+0.5*(norm2(elem_div((log(stcpue+1.e-3)-log(ecpue+1.e-3)),T_var))));  

//Log likelhihood for trawl survey multinomial proportion     
  tf(4) = -(sum(elem_prod(st,colsum(elem_prod(ont,log(ent+1.e-3))))) - sum(elem_prod(st,colsum(elem_prod(ont,log(ont+1.e-3))))));   
//Log likelhihood for trawl survey multinomial proportion     
  tf(5) = -(sum(elem_prod(st,colsum(elem_prod(oot,log(eot+1.e-3))))) - sum(elem_prod(st,colsum(elem_prod(oot,log(oot+1.e-3))))));   
   
//Log likelhihood for winter multinomial proportion   
  tf(6) = -lawp*(sum(elem_prod(sw,colsum(elem_prod(onw,log(enw+1.e-3))))) - sum(elem_prod(sw,colsum(elem_prod(onw,log(onw+1.e-3))))));     
  
//Log likelhihood for winter multinomial proportion   
  tf(7) = -lawp*(sum(elem_prod(sw,colsum(elem_prod(oow,log(eow+1.e-3))))) - sum(elem_prod(sw,colsum(elem_prod(oow,log(oow+1.e-3))))));       
 
//Log likelhihood size proportion for summer fishery survey    
  tf(8) = -(sum(elem_prod(sc,colsum(elem_prod(onc,log(enc+1.e-3))))) - sum(elem_prod(sc,colsum(elem_prod(onc,log(onc+1.e-3))))));      

  //Log likelhihood size proportion for summer fishery survey    
  tf(9) = -(sum(elem_prod(sc,colsum(elem_prod(ooc,log(eoc+1.e-3))))) - sum(elem_prod(sc,colsum(elem_prod(ooc,log(ooc+1.e-3))))));      
  
//Log likelhihood size proportion for observer survey    
  tf(10) = -(sum(elem_prod(so,colsum(elem_prod(ono,log(eno+1.e-3))))) - sum(elem_prod(so,colsum(elem_prod(ono,log(ono+1.e-3))))));   

//Log likelhihood size proportion for observer survey    
  tf(11) = -(sum(elem_prod(so,colsum(elem_prod(ooo,log(eoo+1.e-3))))) - sum(elem_prod(so,colsum(elem_prod(ooo,log(ooo+1.e-3))))));   
  
//deviation in recruits.  
  tf(12) = norm2(log_relrec)/(2*SDRec*SDRec);                            
  
// tag recovery likelihood
  
  ef(1) = -(sum(elem_prod(rowsum(tagrecap1),rowsum(elem_prod(ptagrecap1,log(egr1+1.e-3))))-elem_prod(rowsum(tagrecap1),rowsum(elem_prod(ptagrecap1,log(ptagrecap1+1.e-3)))))); 
  ef(2) = -(sum(elem_prod(rowsum(tagrecap2),rowsum(elem_prod(ptagrecap2,log(egr2+1.e-3))))-elem_prod(rowsum(tagrecap2),rowsum(elem_prod(ptagrecap2,log(ptagrecap2+1.e-3))))));  
  ef(3) = -(sum(elem_prod(rowsum(tagrecap3),rowsum(elem_prod(ptagrecap3,log(egr3+1.e-3))))-elem_prod(rowsum(tagrecap3),rowsum(elem_prod(ptagrecap3,log(ptagrecap3+1.e-3))))));  
  ef(4) = -(sum(elem_prod(rowsum(tagrecap12),rowsum(elem_prod(ptagrecap12,log(egr12+1.e-3))))-elem_prod(rowsum(tagrecap12),rowsum(elem_prod(ptagrecap12,log(ptagrecap12+1.e-3)))))); 
  ef(5) = -(sum(elem_prod(rowsum(tagrecap22),rowsum(elem_prod(ptagrecap22,log(egr22+1.e-3))))-elem_prod(rowsum(tagrecap22),rowsum(elem_prod(ptagrecap22,log(ptagrecap22+1.e-3))))));  
  ef(6) = -(sum(elem_prod(rowsum(tagrecap32),rowsum(elem_prod(ptagrecap32,log(egr32+1.e-3))))-elem_prod(rowsum(tagrecap32),rowsum(elem_prod(ptagrecap32,log(ptagrecap32+1.e-3))))));  
  
  tf(13) = latag*sum(ef);
  
// Total Likelihood
  f += sum(tf);   

//===============================================================================================================
REPORT_SECTION
  cout << "Report Section" << endl;
  report << "nps" << endl << nps << endl;  // Modeled new shell summer 
  report << "ops" << endl << ops << endl;  // Modeled old shell summer
  report << "npw" << endl << npw << endl;  // Modeled new shell winter
  report << "opw" << endl << opw << endl;  // Modeled oldshell shell winter  
  report << "ett" << endl << ett << endl;  // Estimated trawl abundance
  report << "ecpue" << endl << ecpue << endl;  // Estimated Summer fishery cpue
  report << "ewcpue" << endl << ewcpue << endl;  // Estimated Winter survey cpue
  report << "ent" << endl << ent << endl;  // Estimated trawl newshell size proportion 
  report << "eot" << endl << eot << endl;  // Estimated trawl oldshell size proportion   
  report << "enw" << endl << enw << endl;  // Estimated winter survey newshell size proportion
  report << "eow" << endl << eow << endl;  // Estimated winter survey oldshell size proportion
  report << "enc" << endl << enc << endl;  // Estimated summer fishery newshell size proportion
  report << "eoc" << endl << eoc << endl;  // Estimated summer fishery oldshell size proportion
  report << "eno" << endl << eno << endl;  // Estimated observer newshell size proportion
  report << "eoo" << endl << eoo << endl;  // Estimated observer oldshell size proportion
  report << "Rec" << endl << rec << endl;  // Estimated Recruits abundance
  report << "Legal" << endl << legal << endl;  // Estimated legal abundance
  report << "Legalb" << endl << legalb << endl;  // Estimated legal abundance  
  report << "MMB" << endl << mmb << endl;  // Estimated mmb abundance
  report << "ett1" << endl << ett1 << endl;  // Estimated non-legal trawl survey abundance
  report << "ett2" << endl << ett2 << endl;  // Estimated legal trawl survey abundance
  report << "bc" << endl << bc << endl;  // Estimated Summer discards biomass
  report << "bcw" << endl << bcw << endl; // Estimated Winter discards
  report << "f" << endl << f << endl;  // Total likelihood
  // Individual likelihood
  report << "tf" << endl << tf << endl;
  report << "ef" << endl << ef << endl;
  report << "tgr" << endl << tgr << endl;
  report << "mgr" << endl << mgr1 << endl;
  report << "egr1" << endl << egr1 << endl;
  report << "egr2" << endl << egr2 << endl;
  report << "egr3" << endl << egr3 << endl;
  report << "egr12" << endl << egr12 << endl;
  report << "egr22" << endl << egr22 << endl;
  report << "egr32" << endl << egr32 << endl;
 
//===============================================================================================================
GLOBALS_SECTION
  #include <math.h>
  #include <admodel.h>
  #include <time.h>
  ofstream mcmc1,mcmc2;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

//===============================================================================================================
TOP_OF_MAIN_SECTION
  arrmblsize = 10000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  time(&start);

//===============================================================================================================
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