//==========================================================================
// NS3n2017_Feb1_tr_new.TPL
// This model is identical with NNS3n2017_Feb1.TPL.   
// Transpose format column:Lengths - row:years to column:years - row:lenghts   
// Model codes were cleaned and simplifiled 
// This version eliminate q2 and sel93() of commercial catch
// Incorporated discards of Winter commercial and subsistence 
//  
//==========================================================================

//==========================================================================
// 1.0  Data Entry 
//==========================================================================
DATA_SECTION
//--------------------------------------------------------------------------
// Control Data
// fyear	->	 first year
// lyear	->	 last year
// na		->	 number of length classes
// ra		->	 number of recruit length classes
// sla		->	 number of sublegal length classes
// nyr		->	 number of years with trawl survey
// nyw		->	 number of years with winter survey
// nyo		->	 number of years with observer survey
// nysp		->	 number of years with Sprint tagging survey
// ntag1	->	 number of tag data dimention 1976-1993
// ntag2	->	 number of tag data dimention 1993-
//--------------------------------------------------------------------------
  init_int fyear              
  init_int lyear               
  init_int na                 
  init_int ra				  
  init_int sla                
  init_int nyt               
  init_int nyw                
  init_int nyo 
  init_int nysp  
  init_int ntag1     		  
  init_int ntag2              
//--------------------------------------------------------------------------  
// summer commercial fishery event indices: 
// scy		->	 
//		1)year large trawl fishery ended, 
//		2) escape mechanism installed, 
// 		3) year only  CW>5 inch crab was accepted by buyers
// qyear	-> Year noaa survey catchability changed
//--------------------------------------------------------------------------
  init_ivector scy(1,3)       
  init_int qyear              
//--------------------------------------------------------------------------  
// Crab biology  
// slm	->	median length (mm) of smalllest length class 
// sl5	->	length interval (mm)
// M2	->	assumed instantaneous natural mortality  
// msn  ->	natural mortality multiplier for the last two size classes
//--------------------------------------------------------------------------
  init_number slm             //median length of smallest group (mm)
  init_number slt             //length interval (mm)
  init_number M2              //instantaneous natural mortality
  init_vector msn(1,2)	
//--------------------------------------------------------------------------  
// Effective maximum sample size for multinomial length proportion
// maxs	->  
// (1) Commercial length, (2) Survey length
// efn  -> % of effective sample size
//--------------------------------------------------------------------------  
  init_vector maxs(1,2)       
  init_vector efn(1,2)   
//--------------------------------------------------------------------------  
// Other assumptions 
// hm	->	discard mortality: (1) summer, (2) winter 
// lg	->	proportion of legal by length class  
// wm	->	mean weight (lb) by length class   
//--------------------------------------------------------------------------   
  init_vector hm(1,2)         
  init_vector lg(1,na)        
  init_vector wm(1,na)        
//--------------------------------------------------------------------------
// Trawl Survey Data
// 1) it	->   survey year
// 2) tt	->   survey abundance
// 3) pct	->   % of summer catch occurred before the mid point of survey
// 4) yt	->   Mid point of trawl survey from July 1
// 5) cv	->   survey cv
// 6) nont	->   number of sample by length newshell
// 7) noot	->   number of sample by length oldwshell
//--------------------------------------------------------------------------
  init_ivector it(1,nyt)      
  init_vector tt(1,nyt)
  init_vector pct(1,nyt)
  init_vector yt(1,nyt)
  init_vector cv(1,nyt)
  init_matrix nont(1,na,1,nyt)
  init_matrix noot(1,na,1,nyt)

//--------------------------------------------------------------------------
// Winter Pot Survey Data
// 1) iw	->   survey year
// 2) wcpue	->   survey cpue
// 3) nonw	->   number of sample by length newshell
// 4) noow	->   number of sample by length oldwshell
//--------------------------------------------------------------------------
  init_ivector iw(1,nyw)
  init_vector wcpue(1,nyw)
  init_matrix nonw(1,na,1,nyw)
  init_matrix noow(1,na,1,nyw)

//--------------------------------------------------------------------------
// Summer Commercial Catch Data
// 1) tc	->   Annual catch 
// 2) te	->   Annual effort (potlift) Not Use 
// 3) stcpue	->   Standardized CPUE estimated outside of the model 
// 4) secpue	->   Standardized CPUE SE estimated outside of the model
// 5) ys	->   Mid point of summer fishery from July 1
// 6) nonc	->   number of sample by length newshell
// 7) nooc	->   number of sample by length oldwshell
//--------------------------------------------------------------------------
  init_vector tc(fyear,lyear)
  init_vector te(fyear,lyear)
  init_vector stcpue(fyear,lyear)
  init_vector secpue(fyear,lyear)
  init_vector ys(fyear,lyear)
  init_matrix nonc(1,na,fyear,lyear)
  init_matrix nooc(1,na,fyear,lyear)
//--------------------------------------------------------------------------
// Winter Catch Data
// 1) twc	->   Annual Cmmercial catch 
// 2) tws	->   Annual Subsistence Catch Retained
// 3) twst	->   Annual Subsistence Catch Total 
//--------------------------------------------------------------------------
  init_vector twc(fyear,lyear)
  init_vector tws(fyear,lyear)
  init_vector twst(fyear,lyear)
//--------------------------------------------------------------------------
// Summer Commercial Observer dicarded Catch data
// 1) io	->   Year of Observer data 
// 2) nono	->   number of sample by length newshell
// 3) nooo	->   number of sample by length olfshell
//--------------------------------------------------------------------------
  init_ivector io(1,nyo)
  init_matrix nono(1,na,1,nyo)
  init_matrix nooo(1,na,1,nyo)
//--------------------------------------------------------------------------
// Spring Tagging Suivey data
// 1) isp	->   Year of tagging data 
// 2) nonsp	->   number of sample by length newshell
// 3) noosp	->   number of sample by length olfshell
//--------------------------------------------------------------------------
  init_ivector isp(1,nysp)
  init_matrix nonsp(1,na,1,nysp)
  init_matrix noosp(1,na,1,nysp)  
//--------------------------------------------------------------------------
// Tag recovery data data
// 1) tag_recov1	->  tag recovery 1976-1992 
// 2) tag_recov2    ->	tag recovery 1993-present
//--------------------------------------------------------------------------
  init_imatrix tag_recov1(1,ntag1,1,5)
  init_imatrix tag_recov2(1,ntag2,1,5)
//--------------------------------------------------------------------------
// Misc Control data 
// SDRec	-> SD of recruitment 
// SDW		-> SD of winter cpue  
// pwh		-> Proportion of Winter Crab harvest 
//--------------------------------------------------------------------------
  init_number SDRec           // SD of recruitment 
  init_number SDW             // SD of winter cpue    
  init_int M_phase            // M estimat
  init_int rmol_phase
  init_int lamc               // likelifhood for cpue
  init_int lawp               // likelihood weight for winter pot comp
  init_number latag              // likelihood weight for tagging data          
  init_int nst
  init_int nsc
  init_int ssc
  init_vector bse(1,2)
  init_number pwh

  !! cout << "Data Section Completed" << endl;
  !! cout << "lg " << lg << endl;
  !! cout << "twc " << endl << twc << endl;
  !! cout << "ssc " << endl << ssc << endl;  
  !! cout << "bse " << endl << bse << endl;  
  !! cout << "pwh " << endl << pwh << endl;

//==========================================================================
// 2.0  Define Parameters 
//==========================================================================
PARAMETER_SECTION
//--------------------------------------------------------------------------
// Model estimated parameters
// 1) log_q			->	ln Summer Comfish catchability (two peirods)
// 2) log_qw		->	ln Winter Survey catchability   (Not use)
// 3) log_initpop	->	ln Initial population size 
// 4) log_rescale	->	ln mean recruit 
// 5) log_relrec	->	ln recruit SD
// 6) flnp			->	ln Initial length prop
// 7) rlnp			->	ln Initial recruit length prop
// 8) log_amol		->	ln molt prob parameter alpha
// 9) log_bmol		->	ln molt prob parameter beta
// 10) log_amol_dev ->	Random walk molinting paramter alpha
// 11) log_bmol_dev  ->	Random walk molinting paramter beta
// 10) log_ast 		->	ln trawl selectivity logistic alpha
// 11) log_bst 		->	ln trawl selectivity logistic beta
// 12) log_asw		->	ln Winter pot selectivity logistic alpha
// 13) log_bsw		->	ln Winter pot selectivity logistic beta
// 14) sw3 			->	Winber pot selectivity recruit
// 12) log_asc 		->	ln Summer com selectivity logistid alpha 
// 13) log_bsc 		->	ln Summer com selectivity logistid beta 
// 15) advar 		->	Additional variance for st CPUE
// 16) qtno 		->	1976-1991 NOAA Trawl Survey Catchability 
// 17) M 			->	Natural mortality (Estunate)
// 18) ig			->	linear growth parameters 
// 19) ms1			->	Length specific mortality muliplier
//==========================================================================
  init_bounded_vector log_q(1,2,-20.5,10.0,1)
  init_bounded_number log_qw(-10.5,20.0,-1)  
  init_bounded_number log_initpop(2.0,15.0,1)  
  init_bounded_number log_recscale(2.0,12.0,1)
  init_bounded_dev_vector log_relrec(fyear,lyear-1,-40.0,40.0,1)
  init_bounded_vector flnp(1,na-1,0.0,10.0,1)
  init_bounded_vector rlnp(1,ra-1,0.0,10.0,1)
  init_bounded_number log_amol(-5.0,-1.0,1)     
  init_bounded_number log_bmol(1.0,5.5,1)
  init_bounded_dev_vector log_amol_dev(fyear,lyear,-10.0,10.0,rmol_phase)      
  init_bounded_dev_vector log_bmol_dev(fyear,lyear,-10.0,10.0,rmol_phase)     
  init_bounded_vector log_ast(1,nst,-5,1.0,1) 
  init_bounded_vector log_bst(1,nst,0.0,6.0,1)      
  init_bounded_number log_asw(-5.5,-1.0,1)
  init_bounded_number log_bsw(1.5,6.0,1) 
  init_bounded_vector sw3(1,ra,0.,1.0,1)                      
  init_bounded_vector log_asc(1,nsc,-5.0,1.0,-ssc)      
  init_bounded_vector log_bsc(1,nsc,0.0,6.0,-ssc)      
  init_bounded_number advar(0.0,6.0,2) 
  init_bounded_number qtno(0.1,1.0,1)     
  init_bounded_number M(0.02,1.0,M_phase)           
  init_bounded_number sigma(0.,30.,2)    
  init_bounded_vector ig(1,2,0.,20.,1)
  init_bounded_vector ms1(1,na,0.5,5.0,1) 

//--------------------------------------------------------------------------
// Derived Working Parameters: Total ample sizes 
// 1) st		->	Trawl survey
// 2) sw		->	Winter Pot 
// 3) sc		->	Summer Commercial  
// 4) so		->	Summer Commercial Observer discards  
// 5) ssp		->	Spring Pot survey 
// 6) efst		->	Effective sample size Trawl survey
// 7) efsw		->	Effective sample size Winter Pot 
// 8) efsc		->	Effective sample size Summer Commercial  
// 9) efso		->	Effective sample size Summer Commercial Observer discards  
  
//--------------------------------------------------------------------------
  vector st(1,nyt)       
  vector sw(1,nyw)        
  vector sc(fyear,lyear)         
  vector so(1,nyo)  
  vector ssp(1,nysp)  
  vector efst(1,nyt)       
  vector efsw(1,nyw)        
  vector efsc(fyear,lyear)         
  vector efso(1,nyo) 
  vector efsp(1,nysp)  
  
//--------------------------------------------------------------------------
// Derived Observed and Modeled Parameters: Length proportion
// 1) ont oot		->	Observed Trawl survey newshell oldshell
// 2) onw oow		->	Observed Winter Pot survey newshell oldshell
// 3) onc ooc		->	Observed Summer Commercial newshell oldshell
// 4) ono ooo		->	Observed Summer Com discards newshell oldshell
// 5) ent eot		->	Modeled Trawl survey newshell oldshell
// 6) enw eow		->	Modeled Winter Pot survey newshell oldshell
// 7) enc eoc		->	Modeled Summer Commercial newshell oldshell
// 8) eno eoo		->	Modeled Summer Com discards newshell oldshell
// 9) enpw eopw		->	Modeled Winter subsisntence newshell oldshell
// 10) enpwd eopwd	->	Modeled Winter sub discards newshell oldshell
// 11) ensp eosp	->	Modeled Spring Pot survey newshell oldshell
//--------------------------------------------------------------------------
  matrix ont(1,na,1,nyt) 
  matrix oot(1,na,1,nyt) 
  matrix onw(1,na,1,nyw)
  matrix oow(1,na,1,nyw)
  matrix onc(1,na,fyear,lyear)
  matrix ooc(1,na,fyear,lyear)
  matrix ono(1,na,1,nyo)
  matrix ooo(1,na,1,nyo)
  matrix onsp(1,na,1,nysp)                        
  matrix oosp(1,na,1,nysp)   
  matrix ent(1,nyt,1,na)                        
  matrix eot(1,nyt,1,na)                         
  matrix enw(1,nyw,1,na)                        
  matrix eow(1,nyw,1,na)
  matrix enc(fyear,lyear,1,na)
  matrix eoc(fyear,lyear,1,na)
  matrix eno(1,nyo,1,na)
  matrix eoo(1,nyo,1,na)
  matrix enwp(fyear,lyear,1,na) 
  matrix eowp(fyear,lyear,1,na)
  matrix enwpd(fyear,lyear,1,na) 
  matrix eowpd(fyear,lyear,1,na)
  matrix ensp(1,nysp,1,na)                        
  matrix eosp(1,nysp,1,na)  

 //--------------------------------------------------------------------------
// Modeled Abundance, catch, cpue
// 1) nps ops		->	Summer abundance newshell oldshell
// 2) npw opw		->	Winter abundance newshell oldshell
// 3) ent0 eot0		->	Trawl survey abundance newshell oldshell
// 4) ett			->  Trawl survey abundance
// 4) ettq			->  Survey q adjusted Trawl survey abundance
// 4) enwc eowc		->	Winter com catch newshell oldshell
// 5) twsd  		->	Winter sub catch discarded (Observed)
// 6) tb			->  Summer total exploitable abundance
// 7) twp			->  Winter total exploitable abundance
// 8) log_rec		->  ln recruitment 
// 9) rec			->  Recruitment 
// 10) ecpue		->  Summer commercial cpue
// 11) cvcpue		->  CV of cpue (data)
// 12) ewcpue 		->  Winter Pot cpue (do not use)
// 13) bc			->  Summer commercial discards 
// 14) bcw			->  Winter commercial discards
//-------------------------------------------------------------------------- 
  matrix nps(fyear,lyear+1,1,na)    
  matrix ops(fyear,lyear+1,1,na)
  matrix npw(fyear,lyear+1,1,na)
  matrix opw(fyear,lyear+1,1,na)
  matrix ent0(1,nyt,1,na)
  matrix eot0(1,nyt,1,na) 
  vector ett(1,nyt)
  vector ettq(1,nyt)
  matrix enwc(fyear,lyear,1,na)
  matrix eowc(fyear,lyear,1,na)
  vector twsd(fyear,lyear)
  vector tb(fyear,lyear)
//  vector twp(fyear,lyear)
  vector log_rec(fyear,lyear-1)
  vector rec(fyear,lyear)  
  vector ecpue(fyear,lyear)  
  vector cvcpue(fyear,lyear)                       
//  vector ewcpue(1,nyw)
  vector bc(fyear,lyear) 
  vector bcw(fyear,lyear)
//--------------------------------------------------------------------------
// Molting - Selectivity Probability vectors 
// 1) mlen 			->	Mide length for each length class
// 2) molp			->	Annual Molting Probability 
// 3) mp1			->	Mean molting Probability 
// 4) selc1			->  Summer commercial selectivity 
// 5) selc			->  Annual Summer commercial selectivity
// 6) selt1			->  Trawl survey selectivity
// 7) selw			->  winter survey.Commercial selectivity
//--------------------------------------------------------------------------
  vector mlen(1,na)                             
  matrix molp(fyear,lyear,1,na)                 
  vector mp1(1,na)                              
  matrix selc1(1,2,1,na)                       
  matrix sel(fyear,lyear,1,na)                 
  matrix selt1(1,2,1,na)
  matrix selt(fyear,lyear,1,na)
  vector selw(1,na)
//--------------------------------------------------------------------------
// Crab length specific vectors 
// 1) matc 			->	Mature crab
// 2) imat			->	Immature crab 
// 3) Mn			->	Length specific nautral mortality 
// 4) expn			->	First year length proportion: all are newshell
// 5) expr			->  Proportion of recruits
// 6) npp			->	First year length proportion: all are newshell
// 7) rpp			->  Proportion of recruits
// 8) slg			->  Proprtion of sublegal
// 9) wmlg			->  lb.*proportion of legal
// 10) wmslg			->  lb.*proportion of sublegal
// 11) q				->	Commercial Catchability Coeffient (1 and 2)
// 12) qw			->  Winter survey catchability
//--------------------------------------------------------------------------
  vector matc(1,na)                            
  vector imatc(1,na)                           
  vector Mn(1,na)
  vector expn(1,na-1) 
  vector expr(1,ra-1) 
  vector npp(1,na)
  vector rpp(1,na)
  vector slg(1,na)
  vector wmlg(1,na)
  vector wmslg(1,na)
  vector q(1,2)
  //  number qw  
//--------------------------------------------------------------------------
// Likelihood vectors  
// 1) T_var 		->	Additional variance for com cpue
// 2) tf			->	Likelihood vectors  
// 3) offset		->	Length prop liklihood offset  
// 4) toffser		->	Tagging length prop offset
//--------------------------------------------------------------------------  
  vector T_var(fyear,lyear);
  vector tf(1,17);           
  vector offset(1,10);  
  vector toffset(1,6);
//==========================================================================
// Tagging related matrix   
// 1) tag1	 		->	tagging freq before 1993 for Year 1 
// 2) tag2	 		->	tagging freq before 1993 for Year 2  
// 3) tag3	 		->	tagging freq before 1993 for Year 3
// 4) tag12	 		->	tagging freq after 1993 for Year 1 
// 5) tag22	 		->	tagging freq after 1993 for Year 2  
// 6) tag32	 		->	tagging freq after 1993 for Year 3
// 7) ptag1	 		->	prob tagging freq before 1993 for Year 1 
// 8) ptag2	 		->	prob tagging freq before 1993 for Year 2  
// 9) ptag3	 		->	prob tagging freq before 1993 for Year 3
// 10) ptag12	 	->	prob tagging freq after 1993 for Year 1 
// 11) ptag22	 	->	prob tagging freq after 1993 for Year 2  
// 12) ptag32	 	->	prob tagging freq after 1993 for Year 3
// 13) tgr	 		->	Transition growth matrix 
// 14) mgr1	 		->	molting probability adjsteted matrix Year 1
// 15) mgr2	 		->	molting probability adjsteted matrix Year 2
// 16) mgr3	 		->	molting probability adjsteted matrix Year 3
// 17) egr1	 		->	estimated groth matrix for sel1 Year 1
// 18) egr2	 		->	estimated groth matrix for sel1 Year 2
// 19) egr3	 		->	estimated groth matrix for sel1 Year 3
// 20) egr12	 	->	estimated groth matrix for sel2 Year 1
// 21) egr22	 	->	estimated groth matrix for sel2 Year 2
// 22) egr32		->	estimated groth matrix for sel2 Year 3
// 23) ef			->	objective likelihood  
//--------------------------------------------------------------------------  
  matrix tag1(1,na,1,na)
  matrix tag2(1,na,1,na)
  matrix tag3(1,na,1,na)
  matrix tag12(1,na,1,na)
  matrix tag22(1,na,1,na)
  matrix tag32(1,na,1,na)
  matrix ptag1(1,na,1,na) 
  matrix ptag2(1,na,1,na)
  matrix ptag3(1,na,1,na)
  matrix ptag12(1,na,1,na)
  matrix ptag22(1,na,1,na)
  matrix ptag32(1,na,1,na)
  matrix tgr(1,na,1,na)
  matrix mgr1(1,na,1,na)
  matrix mgr2(1,na,1,na)
  matrix mgr3(1,na,1,na)
  matrix egr1(1,na,1,na)
  matrix egr2(1,na,1,na)
  matrix egr3(1,na,1,na)
  matrix egr12(1,na,1,na)
  matrix egr22(1,na,1,na)
  matrix egr32(1,na,1,na)
  vector ef(1,6)         
//==========================================================================

  objective_function_value f
  sdreport_vector last_y(1,na)
  vector last_legalb(1,na)  
  vector last_sublb(1,na)
  sdreport_vector legaln(fyear,lyear+1)
  sdreport_vector legalb(fyear,lyear+1)
  sdreport_vector mmb(fyear,lyear+1)
  sdreport_number last_ofl
  sdreport_number temp
  sdreport_number last_subofl
  sdreport_number bmsy
  sdreport_number last_mmb

  !! cout << "Parameter Section Completed" << endl;

//==========================================================================
INITIALIZATION_SECTION
  log_q         -6.5
  log_amol     -2.4
  log_bmol       4.5
  log_ast         -2.5
  log_bst         1.0
  log_asc         -2.5
  log_bsc         5.5  
  log_asw         -2.5
  log_bsw         4.5
  log_qw         -6.5
  log_amol_dev   0
  log_bmol_dev   0  
  advar 0.5
  qtno 1.0
  sigma  5.0
  ig 10.0
//==========================================================================

//==========================================================================
// 4.0  Preliminary Calculation
//==========================================================================
PRELIMINARY_CALCS_SECTION
  int i,j;
  dvariable tt0,n0; // Calculated working variables nondefferentiated
  
//  cout << "Starting preliminary calcs" << endl;
  M = M2;
  tt0.initialize();
  n0.initialize();
  matc.initialize();

//--------------------------------------------------------------------------
// 4.1   Fishery selectivity weight conversions     
//--------------------------------------------------------------------------
//  Calculate mid length  
  for (i=1;i<=na;i++) mlen(i) = slm + (double(i)-1.0)*slt;
//  Calculate proporion of sublegal
  slg = -lg+1.0;
//  Calculate lb.proprtion lega 
  wmlg  = elem_prod(lg,wm);       
//  Calculate lb.proprtion siblegal 
  wmslg = elem_prod(slg,wm);      
  
  for (i=1;i<=na;i++)
    {
	if (i>ra) matc(i) = 1.0;          // mature crab is 1.0
	}
  imatc = -matc+1.0;                  // inmature crab  

//--------------------------------------------------------------------------
// 4.2  Calculate annual sample size and length proportions     
// 1) st		->	annual sample size trawl survey
// 2) sw		->	annual sample size winter pot survey 
// 3) sc		->	annual sample size Summer Commercial Catch  
// 4) so		->	annual sample size Summer Commercial discards  
// 5) ssp		->	annual sample size Spring pot survey  
//--------------------------------------------------------------------------
  st = colsum(nont)+colsum(noot);
  sw = colsum(nonw)+colsum(noow);
  sc = colsum(nonc)+colsum(nooc);
  so = colsum(nono)+colsum(nooo);
  ssp = colsum(nonsp)+colsum(noosp);
  
// Calculate Length proportions: newshell oldshell
    
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
     onsp(i) = elem_div(nonsp(i),ssp);
     oosp(i) = elem_div(noosp(i),ssp); 	 
	}
	
  
//--------------------------------------------------------------------------
// 4.3  Set Maximum effective sample size: 
//    Trawl & Summer pots:  50% of actual sample size or maxss
//    Winter pots, Spring Pots, Summer commercial, Observer: 
//	  10% of actual sample size or maxsc
//--------------------------------------------------------------------------

// Calculate maximum effective smple sieze for Trawl survey
// efn(1) = 0.5
  efst = st*efn(1);  
  for (i=1;i<=nyt;i++)
   {
    if (efst(i) > maxs(1)) efst(i) = maxs(1);
   }
   
// Calculate maximum effective smple sieze for Winter pot survey
// efn(2) = 0.1
   efsw = sw*efn(2);
  for (i=1;i<=nyw;i++)
   {
     if (efsw(i) > maxs(2)) efsw(i) = maxs(2);
   }

// Calculate maximum effective smple sieze for Spring pot survey
// efn(2) = 0.1
   efsp = ssp*efn(2);
  for (i=1;i<=nysp;i++)
   {
     if (efsp(i) > maxs(2)) efsp(i) = maxs(2);
   }
   
// Calculate maximum effective smple sieze for Summer commercial 
  efsc = sc*efn(2);
  for (i=fyear;i<=lyear;i++)
   {
    if (efsc(i) > maxs(2)) efsc(i) = maxs(2);
   }
   
// Calculate maximum effective smple sieze for Observer
   efso = so*efn(2);  
  for (i=1;i<=nyo;i++)
   {
    if (efso(i) > maxs(2)) efso(i)=maxs(2);
   }

//--------------------------------------------------------------------------
// 4.4  Calculate cv of cpue index      
//--------------------------------------------------------------------------  
  cvcpue = elem_div(secpue+1.e-3,stcpue+1.e-3);

//--------------------------------------------------------------------------
// 4.5  Calculate Winter subsistence discards (twsd)
//      For early years, total number of crab caught is not available. 
//      Estimae discards based on average proporion of discards   
//--------------------------------------------------------------------------
// Discard is a total catch minus retained catch 
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
   
//--------------------------------------------------------------------------
// 4.6  Create tagrecovery matrices       
//--------------------------------------------------------------------------
    tag1.initialize();
    tag2.initialize();
    tag3.initialize();
    tag12.initialize();
    tag22.initialize();
    tag32.initialize();  
    ptag1.initialize();
    ptag2.initialize();
    ptag3.initialize();
    ptag12.initialize();
    ptag22.initialize();
    ptag32.initialize();           
          
// Create tag-recoverey frequency matrix by recovery years          
  for(i=1;i<=ntag1;i++)
  {
        tag1(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,3);
        tag2(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,4);
        tag3(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,5); 
  }
 
  for(i=1;i<=ntag2;i++)
  {
        tag12(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,3);
        tag22(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,4);
        tag32(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,5); 
  }
  
// Create tag-recoverey probability matrix by recovery years   
  if(nsc == 1)
  {
  tag1 = tag1+tag12;
  tag2 = tag2+tag22;
  tag3 = tag3+tag32;
  }
  
  for(i=1;i<=na;i++)
  {  
  ptag1(i)=tag1(i)/(rowsum(tag1)(i)+1.e-3);
  ptag2(i)=tag2(i)/(rowsum(tag2)(i)+1.e-3);
  ptag3(i)=tag3(i)/(rowsum(tag3)(i)+1.e-3);
  ptag12(i)=tag12(i)/(rowsum(tag12)(i)+1.e-3);
  ptag22(i)=tag22(i)/(rowsum(tag22)(i)+1.e-3);  
  ptag32(i)=tag32(i)/(rowsum(tag32)(i)+1.e-3);
  }


//--------------------------------------------------------------------------
// 4.7  Calculate Length proporiton offset used for final likelihood      
//--------------------------------------------------------------------------
  offset(1) = sum(elem_prod(efst,colsum(elem_prod(ont,log(ont+1.e-3)))));
  offset(2) = sum(elem_prod(efst,colsum(elem_prod(oot,log(oot+1.e-3)))));
  offset(3) = sum(elem_prod(efsw,colsum(elem_prod(onw,log(onw+1.e-3)))));
  offset(4) = sum(elem_prod(efsw,colsum(elem_prod(oow,log(oow+1.e-3)))));
  offset(5) = sum(elem_prod(efsc,colsum(elem_prod(onc,log(onc+1.e-3)))));
  offset(6) = sum(elem_prod(efsc,colsum(elem_prod(ooc,log(ooc+1.e-3)))));
  offset(7) = sum(elem_prod(efso,colsum(elem_prod(ono,log(ono+1.e-3)))));
  offset(8) = sum(elem_prod(efso,colsum(elem_prod(ooo,log(ooo+1.e-3)))));
  offset(9) = sum(elem_prod(efsp,colsum(elem_prod(onsp,log(onsp+1.e-3)))));
  offset(10) = sum(elem_prod(efsp,colsum(elem_prod(oosp,log(oosp+1.e-3)))));
    
//--------------------------------------------------------------------------
// 4.8  Calculate Tagging Length offsets       
//--------------------------------------------------------------------------
  toffset(1) = sum(elem_prod(rowsum(tag1),rowsum(elem_prod(ptag1,log(ptag1+1.e-3)))));
  toffset(2) = sum(elem_prod(rowsum(tag2),rowsum(elem_prod(ptag2,log(ptag2+1.e-3)))));
  toffset(3) = sum(elem_prod(rowsum(tag3),rowsum(elem_prod(ptag3,log(ptag3+1.e-3)))));
  toffset(4) = sum(elem_prod(rowsum(tag12),rowsum(elem_prod(ptag12,log(ptag12+1.e-3)))));
  toffset(5) = sum(elem_prod(rowsum(tag22),rowsum(elem_prod(ptag22,log(ptag22+1.e-3)))));
  toffset(6) = sum(elem_prod(rowsum(tag32),rowsum(elem_prod(ptag32,log(ptag32+1.e-3)))));

 cout << "End preliminary calcs" << endl; 

//==========================================================================
RUNTIME_SECTION
  convergence_criteria 1E-6
  maximum_function_evaluations 10000
//==========================================================================

//==========================================================================
// 5.0  Procedure Section      
//==========================================================================
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
//  cout <<"OK for objective function ..."<<endl;  
   
//==========================================================================
// 6.0  Function:  Convert parameters into rates      
//==========================================================================
FUNCTION convert_parameters_into_rates
  int i, j;
// working variables: molting probability
  dvariable mol,asw,bsw,alpha,beta;    
  dvar_vector ast(1,nst),bst(1,nst), asc(1,nsc),bsc(1,nsc);
// Initialize length proportions.   
  npp.initialize();				
  rpp.initialize();				
//catchability coefficient
   q = mfexp(log_q);
//   qw = mfexp(log_qw);   
// Summer Commercial fisheris selectivity parameters    
   asc = mfexp(log_asc);
   bsc = mfexp(log_bsc);
// Summer Trawl Survey selectivity parameters 
   ast = mfexp(log_ast);
   bst = mfexp(log_bst); 
   asw = mfexp(log_asw);
   bsw = mfexp(log_bsw);
   alpha = mfexp(log_amol);
   beta = mfexp(log_bmol);
//--------------------------------------------------------------------------
// 6.1  Set molting and selectivity funuctions      
//--------------------------------------------------------------------------
// Molting Parameter Radom Walk implementaiton   
  mp1 = 1.0/(1.0+mfexp(alpha*(mlen - beta)));
    for (i=fyear;i<=lyear;i++)	
    {	  
	   alpha *= exp(log_amol_dev(i));
	   beta *= exp(log_bmol_dev(i));	 	
       molp(i) = 1.0/(1.0+mfexp(alpha*(mlen - beta)));			
    }
// Summmer Commercial selectivity 
// commercial pot selectivity 1976-1992: logistic function
   selc1(1) = 1.0/(1.0+mfexp(-asc(1)*(mlen - bsc(1)))); 
   selc1(2) = 1.0/(1.0+mfexp(-asc(nsc)*(mlen - bsc(nsc))));  
   if(bse(1) >= 1){    // period 1976 - 1992 
   selc1(1) = 1.0/(1.0+mfexp(asc(1)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));  
   selc1(2) = 1.0/(1.0+mfexp(asc(nsc)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));   
   }   
   selc1(1,na) = 1.0;
   selc1(2,na) = 1.0;
// Winter Pot Survey Selectivity  	
   selw = 1.0/(1.0+mfexp(asw*(mlen-bsw)));  
// Summer Trawl survey Selectivity 
// NOAA summer trawl net selectivity 1976-1992	 
//     selt1(1) = 1.0/(1.0+mfexp(-ast(1)*(mlen-bst(1))));   
// ADFG summer trawl net selectivity: 1996- 	 
//     selt1(2) = 1.0/(1.0+mfexp(-ast(nst)*(mlen-bst(nst))));   
   selt1(1) = 1.0/(1.0+mfexp(ast(1)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));   
   selt1(2) = 1.0/(1.0+mfexp(ast(nst)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));        
   selt1(1,na) = 1.0;
   selt1(2,na) = 1.0;
   
// Directly estimate small length class 
  for (j=1;j<=msn(2);j++) selw(j) = sw3(j);

// Define natural mortality  
// Set mortality be 1.0 for some length crab
  for (j=1;j<=(na-msn(1));j++) ms1(j) = 1.0;
    if(bse(2)>=1) ms1(7)=ms1(8);    
	Mn = M*ms1;
	    
  
// Assign Trawl survey seletivity function in different time periods. 
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(1))    // period 1976 - 1992 
      selt(i) = selt1(1);   //NOAA trawl survey 
      else                          
      selt(i) = selt1(2);   //ADFG trawl survey
    }

// Assign Summer commrecial catch seletivity function in different time periods. 
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(2))    // period 1976 - 2007
        sel(i) = selc1(1); 
	  else            // period 2008 - present
        sel(i) = selc1(2);
    }
  
// First year length proportion    
   for (j=1;j<=(na-1);j++) expn(j) = mfexp(flnp(j));
   npp(1,na-1) = expn/(1+sum(expn));
   npp(na) = 1-sum(npp(1,na-1));
   
// Recruits length proportion   
   for (j=1;j<=(ra-1);j++) expr(j) = mfexp(rlnp(j));
   rpp(1,ra-1) = expr/(1+sum(expr));
   rpp(ra) = 1-sum(rpp(1,ra-1));
   
//==========================================================================
// 7.0  Function: Estimate the First Year (Feb 1st) abundance and 
//		length composition
//      Both abundance and length composition is estimated     
//==========================================================================
FUNCTION get_first_year_abundance
  int i, j;
  dvariable first_y;

// Assign first year abundace: 
  first_y = mfexp(log_initpop);
  
// Assign first year abundance on Feb 1st.  Assume no oldshell crabs     
   for (j=1;j<=na;j++) 
	{
		npw(fyear,j) = npp(j)*first_y;
		opw(fyear,j) = 0.0;
	}
  
// Assign Recruits   
  log_rec = log_relrec+log_recscale;  

  for (i=fyear;i<=(lyear-1);i++) rec(i) = mfexp(log_rec(i));  
// recruits of the last year are assumed to be average of the most recent five years
   rec(lyear) = 0.2*(rec(lyear-1)+rec(lyear-2)+rec(lyear-3)+rec(lyear-4)+rec(lyear-5));
   
//==========================================================================
// 8.0  Function:  Gorwth matrix
//     Estimate Growth transition matrix from tagging data   
//========================================================================== 
  
FUNCTION growth_matrix
  int i, j; 
  dvariable t1, mu, fa, fb;
// Initialize matrix
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
   mu = slm;
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
// growth increment of the last length class is 1.0
   tgr(na,na) = 1;
   
// Normalized the growth proabaiblity  
   for (i=1;i<=na;i++)
   {
    tgr(i) = tgr(i)/rowsum(tgr)(i);
   }
   
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

//==========================================================================
// 8.0  Calculate Population dynamics by length class  
//==========================================================================
FUNCTION get_number_by_size
  int i,j,k;
  dvariable pp, wpm, wpim, tt1, ttt0, ttt1, sb, tw;
  dvar_vector nscaf(1,na), tsc(1,na),tscd(1,na), tt0(1,na);    // Abundance after summer fishery
  dvar_vector sselw(1,na),dselw(1,na),lselw(1,na),dlselw(1,na);
  dvar_vector sselc(1,na),dselc(1,na); 
  dvar_vector nwc(1,na),owc(1,na),nwcd(1,na),owcd(1,na);
  dvar_vector nws(1,na),ows(1,na),nwsd(1,na),owsd(1,na); 
  dvar_vector par1(1,na),par2(1,na), FOFL(1,na);   
  dvariable TotalCatch, CatNum, cofl,fl;    
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
  last_ofl.initialize();
  last_subofl.initialize();
// 
// Calculate winter subsistence harvest selectivity
// Assume subsistence harvets take all mature crab caught
   sselw = elem_prod(matc,selw); 
// Assume subsistence harvest discards all immature crab caught   
   dselw = elem_prod(imatc,selw); 
// Calculate legal crab caught by winter gear    
   lselw = elem_prod(lg,selw);   
// Calculate sublegal crab caught by winter gear (discards)  
   dlselw = elem_prod(slg,selw);      

//==========================================================================
//==========================================================================
// 8.0  Population Dynamics 
//==========================================================================
//==========================================================================
  for (i=fyear;i<=lyear;i++)
  {
//--------------------------------------------------------------------------
// 8.1  Winter Fishery 
//      Winter Newshell: npw, Oldshell: opw
//--------------------------------------------------------------------------

//==========================================================================
// 8.1 Winter commercial catch and discards size composition 
//==========================================================================
	    	       
// Calculate Total Winter Legal Crab Catchable to commercial fisheries 
      tw = sum(elem_prod(npw(i)+opw(i),lselw));	  
// Calculate proportion of newshell and oldshell in commercial catch   	  
      enwc(i) = elem_prod(npw(i),lselw)/tw;
      eowc(i) = elem_prod(opw(i),lselw)/tw;
// Calculate number of newshell and oldshell in commercial catch   	  
      nwc = twc(i)*enwc(i);
      owc = twc(i)*eowc(i);
// Calculate number of newshell and oldshell discards in commercial  	  
      nwcd = (twc(i)/tw)*elem_prod(npw(i),dlselw);
      owcd = (twc(i)/tw)*elem_prod(npw(i),dlselw);
// Calculate total number of winter commercial discards mortality
	  bcw(i) = sum(nwcd+owcd)*hm(2); 	  
// Calculate total Winter Crab (Leal + Sublegal) Catchable to commercial
//      twp(i) = sum(elem_prod(npw(i)+opw(i),selw)); //   
   
// Calculate exploitable abundance at the middle of the season
//      twp(i) = twp(i) - 0.5*(twc(i)+tws(i));  

//==========================================================================
// 8.2 Winter subsistence catch and discards size composition 
//     Assume that Subsistence fishers took all length >= 3 crab  
//	   and discarded length 1 & 2 crab
//==========================================================================

// Calculate total mature population catchable to winter fisheries 
     wpm = sum(elem_prod(npw(i)+opw(i),sselw));
// Calculate total immature population cathable to winter sfisheries 
     wpim = sum(elem_prod(npw(i)+opw(i),dselw));
// Calculate proportiojn of newshell and oldshell in winter subsistence  
	 enwp(i) = elem_prod(npw(i),sselw)/wpm;
     eowp(i) = elem_prod(opw(i),sselw)/wpm;
// Calculate number of newshell and oldshell in winter subsistence 
	 nws = tws(i)*enwp(i);
     ows = tws(i)*eowp(i); 
// Calculate proportiojn of new and old shell discards in winter subsistence  
	 enwpd(i) = elem_prod(npw(i),dselw)/wpim;
     eowpd(i) = elem_prod(opw(i),dselw)/wpim;
// Calculate number of new and old shell discards in winter subsistence  
	 nwsd = twsd(i)*enwpd(i);
     owsd = twsd(i)*eowpd(i); 	 
	 
//--------------------------------------------------------------------------
// 8.3  Summer Population: July 1st
//--------------------------------------------------------------------------

//==========================================================================
// 8.4 Calculate abundance for summber July 1st:
//     Summer population is survivors of winter fisheriem 
//     and natural mortality
//==========================================================================
    nps(i) = elem_prod(npw(i)-nwc-nws-nwcd*hm(2)-nwsd*hm(2),exp(-0.417*Mn));  
    ops(i) = elem_prod(opw(i)-owc-ows-owcd*hm(2)-owsd*hm(2),exp(-0.417*Mn));  

// Make sure that abundance does not bellow zero
    for (j=1;j<=na;j++)
	{
	 if (nps(i,j) < 0.0) nps(i,j) = 0.001; 
	 if (ops(i,j) < 0.0) ops(i,j) = 0.001; 	
	}

//==========================================================================
// 8.4  Calculate Commerical catch  
//==========================================================================
// Calculate summer crab legal and sublegal crab selectivity
   sselc = elem_prod(sel(i),lg);   
   dselc = elem_prod(sel(i),slg);     
// Total number of crab catchable to commercial fishery    
   tb(i) = sum(elem_prod(nps(i)+ops(i),sselc));	

// Calculate proprotion of newshell and oldshell in commercial		
	enc(i) = elem_prod(nps(i),sselc)/tb(i);
	eoc(i) = elem_prod(ops(i),sselc)/tb(i);
// Put zero to years of no crab fisheries
	enc(fyear) = 0.0;
    eoc(fyear) = 0.0;
	enc(1991) = 0.0;
    eoc(1991) = 0.0;
// Calculate Summer commercial legal and sublegal catch by length;
	tsc = tc(i)*(enc(i)+eoc(i));
    tscd = (tc(i)/tb(i))*elem_prod(nps(i)+ops(i),dselc);
// Commercial discards biomass 	
    bc(i) = hm(1)*sum(elem_prod(tscd,wm));     

//==========================================================================
//==========================================================================
// 8.5  Winter Population: Feb 1st next year
//==========================================================================
//==========================================================================
	
//==========================================================================
// 8.5  Calculate number of Summer Crab after summer fishery (nscaf)
//      nscaf = July 1st abundance*mortality till fishery 
//               - com catch and discards. 
//==========================================================================
  		
    nscaf = elem_prod(nps(i)+ops(i),exp(-ys(i)*Mn))-tsc-hm(1)*tscd;	  
	
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
      npw(i+1,j) = pp*exp(-(0.583-ys(i))*Mn(j)) + rpp(j)*rec(i);    
	  //Old shell crab are unmolted 
      opw(i+1,j) = nscaf(j)*(1.0-molp(i,j))*exp(-(0.583-ys(i))*Mn(j));    
     }
   }


//==========================================================================
// 8.7  Calculate February 1st legal Crab and mature male biomass  for assessment
//==========================================================================

 for (i=fyear; i<=lyear+1; i++)
   {
    mmb(i) = sum(elem_prod(elem_prod(npw(i)+opw(i),matc),wm));  // Mature male biomass
    legaln(i) = sum(elem_prod(npw(i)+opw(i),lg));  // number of legal crab
    legalb(i) = sum(elem_prod(npw(i)+opw(i),wmlg)); // legal biomass   
	}
   
  last_y = npw(lyear+1)+opw(lyear+1);

// Calculate legal b exploitable to fisheries 
  last_legalb = elem_prod(elem_prod(last_y,wmlg),selc1(2));
  last_sublb = elem_prod(elem_prod(last_y,wmslg),selc1(2));  
  last_mmb= mmb(lyear+1);
  bmsy = (sum(mmb) - mmb(1976)-mmb(1977)-mmb(1978)-mmb(1979))/(lyear-fyear-2);
  cofl = last_mmb/bmsy;
  if (cofl < 0.25){
	  fl = 0;  
        }
	   else if (cofl <=1){ fl = (cofl-0.1)/0.9;
	   }
	   else{fl = 1;
	}
// Calculate FOFL	
  FOFL = fl*0.18;
  par1 = 1-exp(-FOFL-0.42*Mn);
  par2 = 1-exp(-0.42*Mn);
  temp = sum(elem_prod(last_legalb,par1-elem_prod(par2,elem_div(1-pwh*par1,1-pwh*par2)))); 
  for(k=1;k<=na;k++)
	{		
  last_ofl += last_legalb(k)*(1-exp(-FOFL(k)-0.42*Mn(k))-(1-exp(-0.42*Mn(k)))*((1-pwh*(1-exp(-FOFL(k)-0.42*Mn(k))))/(1-pwh*(1-exp(-0.42*Mn(k))))));
  last_subofl += last_sublb(k)*(1-exp(-M-0.42*Mn(k))-(1-exp(-0.42*Mn(k)))*((1-pwh*(1-exp(-M-0.42*Mn(k))))/(1-pwh*(1-exp(-0.42*Mn(k))))));
	} 
  
//==========================================================================
// 9.0  produce model estimated observation data  
//==========================================================================
FUNCTION get_proportion_and_effort
  int i,j;
  dvariable bf,af,pp;
  ett.initialize();
  eow.initialize();
  eoo.initialize();
  ecpue.initialize();

//==========================================================================
// 9.1  Predicted Trawl survey abundance 
//      Survey abundance is adjusted by fishery and mortality 
//      When mid-point of survey date is later than that of commercial fishery. 
//         Remove commercial catch from survey abundance. 
//==========================================================================
  for (i=1;i<=nyt;i++)
   {
 
   if (yt(i) > ys(it(i))) // Mid-point of survey date is later than that of commercial fishery
     {
      bf = ys(it(i));             //time lag from July 1 to fishery
      af = yt(i) - ys(it(i));     //time lag from fishery to survey
     }
    else
     {
      bf = yt(i);
      af = 0.0;
     }
// Calculate trawl survey abundance by Shell, length
    ent0(i) = elem_prod(elem_prod((elem_prod(nps(it(i)),exp(-bf*Mn))-pct(i)*tc(it(i))*enc(it(i))),exp(-af*Mn)),selt(it(i)));
    eot0(i) = elem_prod(elem_prod((elem_prod(ops(it(i)),exp(-bf*Mn))-pct(i)*tc(it(i))*eoc(it(i))),exp(-af*Mn)),selt(it(i)));
// Stop gap measure	 
	for (j=1;j<=na;j++)
     {
	  if (ent0(i,j) < 0.0) ent0(i,j) = 0.0;
	  if (eot0(i,j) < 0.0) eot0(i,j) = 0.0;	  
     }
// Calculate trawl survey abundance 
	 ett(i) = sum(ent0(i)+eot0(i));
// Stop gap measuer
    if (ett(i) <= 0.0) ett(i) = 0.00001;

// Calculate proprion of new and old shell crab by trawl survey	
    ent(i) = ent0(i)/ett(i);
	eot(i) = eot0(i)/ett(i);
 
// Estimate Survey Q for NOAA trawls
  if(it(i) < qyear) {ettq(i) = qtno*ett(i);}
    else
    ettq(i) = (ett(i));
   }

//==========================================================================
// 9.2  Predicted Winter Pot survey length proportion 
//==========================================================================
  for (i=1;i<=nyw;i++)
   {
// Calculate total winter pot catchable 
    pp = sum(elem_prod(npw(iw(i))+opw(iw(i)),selw));
    if (pp <= 0.0) pp = 0.000001;
// Calculate proportion of pot survey shell length  
	enw(i) = elem_prod(npw(iw(i)),selw)/pp;
	eow(i) = elem_prod(opw(iw(i)),selw)/pp;
// Calculate winter survey cpue     
//      ewcpue(i) = qw*twp(iw(i)); 	 Not Use 
   }

//==========================================================================
// 9.2  Predicted Spting Pot survey length proportion 
//==========================================================================
  if(lawp>=1){
  for (i=1;i<=nysp;i++)
   {   
// Calculate total winter pot catchable 
    pp = sum(elem_prod(nps(isp(i))+nps(isp(i)),selw));
    if (pp <= 0.0) pp = 0.000001;
// Calculate proportion of pot survey shell length  
	ensp(i) = elem_prod(nps(isp(i)),selw)/pp;
	eosp(i) = elem_prod(ops(isp(i)),selw)/pp;
   }
  }
//==========================================================================
// 9.3  Predicted Observer survey length proportion  
//==========================================================================
  for (i=1;i<=nyo;i++)
   {
// Calculate total summer com sublegal   
	pp = sum(elem_prod(elem_prod(nps(io(i))+ops(io(i)),sel(io(i))),slg));
    if (pp <= 0.0) pp = 0.000001;
// Calculate proportion of sublegal shell, length 
    eno(i) = elem_prod(elem_prod(nps(io(i)),sel(io(i))),slg)/pp;
    eoo(i) = elem_prod(elem_prod(ops(io(i)),sel(io(i))),slg)/pp;
   }

//==========================================================================
// 9.4  Predicted CPUE
//==========================================================================

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
// if stcpue is zero (no fishery)  ecpue = 0;
	if (stcpue(i) <= 0.0) ecpue(i) = 0.0;
  }
  

//==========================================================================
// 10.0  Likelihood Calculation  
//==========================================================================
FUNCTION evaluate_the_objective_function

  int i,j,k;
// Log likelihood for trawl survey 
  tf(1) = 0.5*norm2(elem_div((log(tt+1.e-3)-log(ettq+1.e-3)),sqrt(log(elem_prod(cv,cv)+1.0))));

// Log likelihood standard cpue in winter fishery 
//  tf(2) = norm2(log(wcpue+1.e-3)-log(ewcpue+1.e-3))/(2*SDW*SDW); 
   
// Log likelihood standard cpue in summer fishery 
  T_var = sqrt(log(elem_prod(cvcpue,cvcpue)+1.0)+advar); 
  tf(3) = lamc*(sum(log(T_var))+0.5*(norm2(elem_div((log(stcpue+1.e-3)-log(ecpue+1.e-3)),T_var))));  

//Log likelhihood for trawl survey multinomial proportion     
  tf(4) = -(sum(elem_prod(efst,rowsum(elem_prod(trans(ont),log(ent+1.e-3))))) - offset(1))      
          -(sum(elem_prod(efst,rowsum(elem_prod(trans(oot),log(eot+1.e-3))))) - offset(2));   
   
//Log likelhihood for winter multinomial proportion   
  tf(5) = -(sum(elem_prod(efsw,rowsum(elem_prod(trans(onw),log(enw+1.e-3))))) - offset(3))     
          -(sum(elem_prod(efsw,rowsum(elem_prod(trans(oow),log(eow+1.e-3))))) - offset(4));       
 
//Log likelhihood size proportion for summer fishery survey    
  tf(6) = -(sum(elem_prod(efsc,rowsum(elem_prod(trans(onc),log(enc+1.e-3))))) - offset(5))         
          -(sum(elem_prod(efsc,rowsum(elem_prod(trans(ooc),log(eoc+1.e-3))))) - offset(6));      
  
//Log likelhihood size proportion for observer survey    
  tf(7) = -(sum(elem_prod(efso,rowsum(elem_prod(trans(ono),log(eno+1.e-3))))) - offset(7)) 
          -(sum(elem_prod(efso,rowsum(elem_prod(trans(ooo),log(eoo+1.e-3))))) - offset(8));   

//Log likelhihood for Spring survey  multinomial proportion   
  if(lawp>=1){
  tf(8) = -(sum(elem_prod(efsp,rowsum(elem_prod(trans(onsp),log(ensp+1.e-3))))) - offset(9));     
          -(sum(elem_prod(efsp,rowsum(elem_prod(trans(oosp),log(eosp+1.e-3))))) - offset(10));       
	}	  
//deviation in recruits.  
  tf(9) = norm2(log_relrec)/(2*SDRec*SDRec);                            
//deviation in molting parameter.  
  tf(10) = norm2(log_amol_dev)/(2*SDRec*SDRec);     
//deviation in molting parameter.  
  tf(11) = norm2(log_bmol_dev)/(2*SDRec*SDRec);     

  
// tag recovery likelihood
  
  ef(1) = -(sum(elem_prod(rowsum(tag1),rowsum(elem_prod(ptag1,log(egr1+1.e-3)))))- toffset(1)); 
  ef(2) = -(sum(elem_prod(rowsum(tag2),rowsum(elem_prod(ptag2,log(egr2+1.e-3)))))-toffset(2));  
  ef(3) = -(sum(elem_prod(rowsum(tag3),rowsum(elem_prod(ptag3,log(egr3+1.e-3)))))-toffset(3));  
  ef(4) = -(nsc-1)*(sum(elem_prod(rowsum(tag12),rowsum(elem_prod(ptag12,log(egr12+1.e-3)))))-toffset(4)); 
  ef(5) = -(nsc-1)*(sum(elem_prod(rowsum(tag22),rowsum(elem_prod(ptag22,log(egr22+1.e-3)))))-toffset(5));  
  ef(6) = -(nsc-1)*(sum(elem_prod(rowsum(tag32),rowsum(elem_prod(ptag32,log(egr32+1.e-3)))))-toffset(6));  
  
  tf(12) = latag*sum(ef);

  
// Total Likelihood
  f = sum(tf);   

//==========================================================================
REPORT_SECTION
  cout << "Report Section" << endl;
  report << "nps" << endl << trans(nps) << endl;  // Modeled new shell summer 
  report << "ops" << endl << trans(ops) << endl;  // Modeled old shell summer
  report << "npw" << endl << trans(npw) << endl;  // Modeled new shell winter
  report << "opw" << endl << trans(opw) << endl;  // Modeled oldshell shell winter  
  report << "ett" << endl << ett << endl;  // Estimated trawl abundance
  report << "ecpue" << endl << ecpue << endl;  // Estimated Summer fishery cpue
//  report << "ewcpue" << endl << ewcpue << endl;  // Estimated Winter survey cpue
  report << "ent" << endl << trans(ent) << endl;  // Estimated trawl newshell size proportion 
  report << "eot" << endl << trans(eot) << endl;  // Estimated trawl oldshell size proportion   
  report << "enw" << endl << trans(enw) << endl;  // Estimated winter survey newshell size proportion
  report << "eow" << endl << trans(eow) << endl;  // Estimated winter survey oldshell size proportion
  report << "enc" << endl << trans(enc) << endl;  // Estimated summer fishery newshell size proportion
  report << "eoc" << endl << trans(eoc) << endl;  // Estimated summer fishery oldshell size proportion
  report << "eno" << endl << trans(eno) << endl;  // Estimated observer newshell size proportion
  report << "eoo" << endl << trans(eoo) << endl;  // Estimated observer oldshell size proportion
  report << "ensp" << endl << trans(ensp) << endl;  // Estimated spring Pot survey newshell size proportion
  report << "eosp" << endl << trans(eosp) << endl;  // Estimated spring Pot survey oldshell size proportion
  report << "Rec" << endl << rec << endl;  // Estimated Recruits abundance
  report << "Legal" << endl << legaln << endl;  // Estimated legal abundance
  report << "Legalb" << endl << legalb << endl;  // Estimated legal abundance  
  report << "MMB" << endl << mmb << endl;  // Estimated mmb abundance
  report << "bc" << endl << bc << endl;  // Estimated Summer discards biomass
  report << "bcw" << endl << bcw << endl; // Estimated Winter discards
  report << "npp" << endl << npp << endl;
  report << "rpp" << endl << rpp << endl;  
  report << "f" << endl << f << endl;  // Total likelihood
  // Individual likelihood
  report << "tf" << endl << tf << endl;
  report << "ef" << endl << ef << endl;
// selectivity
  report << "selc1" << endl << selc1 << endl;
  report << "selt1" << endl << selt1 << endl; 
  report << "selw" << endl << selw << endl;
  report << "M" << endl << Mn << endl;  
  report << "tgr" << endl << tgr << endl;
  report << "mgr" << endl << mgr1 << endl; 
  report << "egr1" << endl << egr1 << endl;
  report << "egr2" << endl << egr2 << endl;
  report << "egr3" << endl << egr3 << endl;
  report << "egr12" << endl << egr12 << endl;
  report << "egr22" << endl << egr22 << endl;
  report << "egr32" << endl << egr32 << endl;
  report << "molp" << endl << trans(molp) <<endl;
  report << "selt" << endl << trans(selt) <<endl; 
  report << "selc" << endl << trans(sel) <<endl; 

  
 
//==========================================================================
GLOBALS_SECTION
  #include <math.h>
  #include <admodel.h>
  #include <time.h>
  ofstream mcmc1,mcmc2;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

//==========================================================================
TOP_OF_MAIN_SECTION
  arrmblsize = 10000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  time(&start);

//==========================================================================
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