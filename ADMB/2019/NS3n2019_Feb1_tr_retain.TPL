//==========================================================================
// NS3n2019_Feb1_tr_new.TPL
// This model is identical with NNS3n2018_Feb1.TPL.   
// Transpose format column:Lengths - row:years to column:years - row:lengths   
// Model codes were cleaned and simplified 
// This version eliminate q2 and sel93() of commercial catch
// Incorporated discards of Winter commercial and subsistence 
// Included Tier 3 F35 calculation  
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
// nyod		->	 number of years with observer discard survey
// nyot		->	 number of years with observer total catch survey
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
  init_int nycw  
  init_int nyod
  init_int nyot  
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
// 1) iod	->   Year of Observer data 
// 2) nonod	->   number of sample by length newshell
// 3) noood	->   number of sample by length olfshell
//--------------------------------------------------------------------------
  init_ivector iod(1,nyod)
  init_matrix nonod(1,na,1,nyod)
  init_matrix noood(1,na,1,nyod)
//--------------------------------------------------------------------------
// Summer Commercial Observer total  Catch data
// 1) iot	->   Year of Observer data 
// 2) nonot	->   number of sample by length newshell
// 3) nooot	->   number of sample by length olfshell
//--------------------------------------------------------------------------
  init_ivector iot(1,nyot)
  init_matrix nonot(1,na,1,nyot)
  init_matrix nooot(1,na,1,nyot)  
//--------------------------------------------------------------------------
// Winter Commercial Retained Catch data
// 1) icw	->   Year of  data 
// 2) noncw	->   number of sample by length newshell
// 3) noocw ->   number of sample by length olfshell
//--------------------------------------------------------------------------
  init_ivector icw(1,nycw)
  init_matrix noncw(1,na,1,nycw)
  init_matrix noocw(1,na,1,nycw)    
  
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
  init_vector SD(1,2)   // SD of recruitment & Winter CPUE   
  init_ivector base(1,6)  // Base Model Configulatoins 
// 1: Summer Commercial logistic function: 1: one parameter, 2: two parameter
// 2: Estimate Summer Commercial retention curve -1: No, 1: Yes (phase)
// 3: Estimate Winter Commercial retention curve -1: No, 1: Yes (phase)
// 4: Estimate Anual Molting -1: No, 1: Yes (phase)
// 5: Estimate M  -1: No, 1: yes (phase)  
// 6: M of largest 2 classes: 0: separate, 1: same   
  init_ivector nsel(1,5)  // Number of selectivity functions/Periods
// 1: Trawl survey
// 2: Summer Commercial
// 3: Summer Commercial Retention   
  init_vector likew(1,5) // Likelihood weights 
// 1: Tag recovery 
// 2: Spring Survey  
// 3: hist: discards   
// 4:
  init_number pwh   // Proportion of Winter Commerical Catch
// Working variales         
  init_int nst   // Trawl survey 
  init_int nsc   // Commercial Periods 
  init_int nscr  // Retain 
  init_int ssc
  !! nst  = nsel(1);
  !! nsc  = nsel(2);
  !! nscr = nsel(3);
  
  !! cout << "Data Section Completed" << endl;
  !! cout << "lg " << lg << endl;
  !! cout << "twc " << endl << twc << endl;
  !! cout << "nsel " << endl << nsel << endl;  
  !! cout << "nsc " << endl << nsc << endl;  
  !! cout << "nscr " << endl << nscr << endl; 
  !! cout << "base " << endl << base << endl;    
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
// 14) log_ascr 	->	ln Summer com retention selectivity logistid alpha 
// 15) log_bscr		->	ln Summer com retention selectivity logistid beta 
// 14) log_aswr 	->	ln Winter com retention selectivity logistid alpha 
// 15) log_bswr		->	ln Winter com retention selectivity logistid beta 
// 16) advar 		->	Additional variance for st CPUE
// 17) qtno 		->	1976-1991 NOAA Trawl Survey Catchability 
// 18) M 			->	Natural mortality (Estunate)
// 19) ig			->	linear growth parameters 
// 20) ms1			->	Length specific mortality muliplier
//==========================================================================
  init_bounded_vector log_q(1,2,-20.5,10.0,1)
  init_bounded_number log_qw(-10.5,20.0,-1)  
  init_bounded_number log_initpop(2.0,15.0,1)  
  init_bounded_number log_recscale(2.0,12.0,1)
  init_bounded_dev_vector log_relrec(fyear,lyear-1,-50.0,50.0,1)
  init_bounded_vector flnp(1,na-1,0.0,10.0,1)
  init_bounded_vector rlnp(1,ra-1,0.0,10.0,1)
  init_bounded_number log_amol(-5.0,-1.0,1)     
  init_bounded_number log_bmol(1.0,5.5,1)
  !! phz = (int) base(4);  
  init_bounded_dev_vector log_amol_dev(fyear,lyear,-10.0,10.0,phz)      
  init_bounded_dev_vector log_bmol_dev(fyear,lyear,-10.0,10.0,phz)     
  init_bounded_vector log_ast(1,nst,-5,1.0,1) 
  init_bounded_vector log_bst(1,nst,0.0,6.0,1)
  init_bounded_number log_asw(-5.5,-1.0,1)
  init_bounded_number log_bsw(1.5,6.0,1)   
  init_bounded_vector sw3(1,ra,0.,1.0,1)    
  init_bounded_vector log_asc(1,nsc,-5.0,1.0,1)      
  init_bounded_vector log_bsc(1,nsc,0.0,7.0,1)   
  !! phz = (int) base(2);
  init_bounded_vector log_ascr(1,nscr,-5.0,1.0,phz)      
  init_bounded_vector log_bscr(1,nscr,0.0,7.0,phz)  
  !! phz = (int) base(3);  
  init_bounded_number log_aswr(-5.0,1.0,phz)      
  init_bounded_number log_bswr(0.0,6.0,phz)    
  init_bounded_number advar(0.0,6.0,2) 
  init_bounded_number qtno(0.5,1.0,1) 
  !! phz = (int) base(5);    
  init_bounded_number M(0.02,1.0,phz)           
  init_bounded_number sigma(0.,30.,2)    
  init_bounded_vector ig(1,2,0.,20.,1)
  init_bounded_vector ms1(1,na,0.5,5.0,1) 

//--------------------------------------------------------------------------
// Derived Working Parameters: Total ample sizes 
// 1) st		->	Trawl survey
// 2) sw		->	Winter Pot 
// 3) sc		->	Summer Commercial Retain  
// 4) sod		->	Summer Commercial Observer discards  
// 5) so5		->	Summer Commercial Observer total catch  
// 6) scw		->	Winter Commercial Retain catch  
// 7) ssp		->	Spring Pot survey 
// 8) efst		->	Effective sample size Trawl survey
// 9) efsw		->	Effective sample size Winter Pot 
// 10) efsc		->	Effective sample size Summer Commercial  
// 11) efsod	->	Effective sample size Summer Commercial Observer discards  
// 12) efsot	->	Effective sample size Summer Commercial Observer total catch
// 13) efsop	->	Effective sample size Spring Survey 
//--------------------------------------------------------------------------
  vector st(1,nyt)       
  vector sw(1,nyw)        
  vector sc(fyear,lyear)         
  vector sod(1,nyod)
  vector sot(1,nyot) 
  vector scw(1,nycw)  
  vector ssp(1,nysp)  
  vector efst(1,nyt)       
  vector efsw(1,nyw)        
  vector efsc(fyear,lyear)         
  vector efsod(1,nyod)
  vector efsot(1,nyot)  
  vector efsp(1,nysp)
  vector efcw(1,nycw)  
  
//--------------------------------------------------------------------------
// Derived Observed and Modeled Parameters: Length proportion
// 1) ont oot		->	Observed Trawl survey newshell oldshell
// 2) onw oow		->	Observed Winter Pot survey newshell oldshell
// 3) onc ooc		->	Observed Summer Commercial newshell oldshell
// 4) onod oood		->	Observed Summer Com discards newshell oldshell
// 5) onot ooot		->	Observed Summer Com total catch newshell oldshell
// 6) onwc oowc		->	Observed Winter Com retained catch newshell oldshell
// 7) onsp oosp		->	Observed Spring Survery newshell oldshell
// 8) ent eot		->	Modeled Trawl survey newshell oldshell
// 9) enw eow		->	Modeled Winter Pot survey newshell oldshell
// 10) enc eoc		->	Modeled Summer Commercial newshell oldshell
// 11) enod eood	->	Modeled Summer Com discards newshell oldshell
// 12) enot eoot	->	Modeled Summer Com total catch newshell oldshell
// 13) encw eocw	->	Modeled Winter Com retain newshell oldshell
// 16) ensp eosp	->	Modeled Spring Pot survey newshell oldshell
//--------------------------------------------------------------------------
  matrix ont(1,na,1,nyt) 
  matrix oot(1,na,1,nyt) 
  matrix onw(1,na,1,nyw)
  matrix oow(1,na,1,nyw)
  matrix onc(1,na,fyear,lyear)
  matrix ooc(1,na,fyear,lyear)
  matrix onod(1,na,1,nyod)
  matrix oood(1,na,1,nyod)
  matrix onod1(1,na,1,6)
  matrix oood1(1,na,1,6) 
  matrix onot(1,na,1,nyot)
  matrix ooot(1,na,1,nyot)
  matrix oncw(1,na,1,nycw)
  matrix oocw(1,na,1,nycw) 
  matrix onsp(1,na,1,nysp)                        
  matrix oosp(1,na,1,nysp)   
  matrix ent(1,nyt,1,na)                        
  matrix eot(1,nyt,1,na)                         
  matrix enw(1,nyw,1,na)                        
  matrix eow(1,nyw,1,na)
  matrix enc(fyear,lyear,1,na)
  matrix eoc(fyear,lyear,1,na)
  matrix enod(1,nyod,1,na)
  matrix eood(1,nyod,1,na)
  matrix enod1(1,6,1,na)
  matrix eood1(1,6,1,na)
  matrix enot(1,nyot,1,na)
  matrix eoot(1,nyot,1,na)
  matrix encw(1,nycw,1,na)
  matrix eocw(1,nycw,1,na)  
  matrix ensp(1,nysp,1,na)                        
  matrix eosp(1,nysp,1,na)  

 //--------------------------------------------------------------------------
// Modeled Abundance, catch, cpue
// 1) nps ops		->	Summer abundance newshell oldshell
// 2) npw opw		->	Winter abundance newshell oldshell
// 3) ent0 eot0		->	Trawl survey abundance newshell oldshell
// 4) ett			->  Trawl survey abundance
// 4) ettq			->  Survey q adjusted Trawl survey abundance
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
  vector ett(1,nyt)
  vector ettq(1,nyt)
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
  vector bswd(fyear,lyear)   
  vector FSub(fyear,lyear)
  vector FSubd(fyear,lyear) 

//--------------------------------------------------------------------------
// Molting - Selectivity Probability vectors 
// 1) mlen 			->	Mid length for each length class
// 2) molp			->	Annual Molting Probability 
// 3) mp1			->	Mean molting Probability 
// 4) selc			->  Summer commercial selectivity 
// 5) selcy			->  Annual Summer commercial selectivity
// 6) selry			->  Annual Summer commercial retention selectivity
// 7) selt			->  Trawl survey selectivity
// 8) selty			->  Annual Trawl survey selectivity
// 9) selw			->  winter survey.Commercial selectivity
// 10) selwr		->  winter Commercial retention selectivity
// 11) lselw		->  winter Commercial retained selectivity
// 12) dlselw		->  winter Commercial discarded selectivity
//--------------------------------------------------------------------------
  vector mlen(1,na)                             
  matrix molp(fyear,lyear,1,na)                 
  vector mp1(1,na)                              
  matrix selc(1,nsc,1,na)                       
  matrix selr(1,nscr,1,na)
  matrix selcy(fyear,lyear,1,na) 
  matrix selry(fyear,lyear,1,na)   
  matrix selt(1,2,1,na)
  matrix selty(fyear,lyear,1,na)
  vector selw(1,na)
  vector selwr(1,na)
  vector lselw(1,na)
  vector dlselw(1,na)
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
// 10) wmslg	    ->  lb.*proportion of sublegal
// 11) q			->	Commercial Catchability Coeffient (1 and 2)
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
  vector tf(1,20);           
  vector offset(1,20);  
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
  vector last_legalb(1,na)  
  vector last_sublb(1,na)
  vector legaln(fyear,lyear+1)
  vector legalb(fyear,lyear+1)
  vector mmb(fyear,lyear+1)
  vector last_y(1,na)
//  sdreport_number last_legal  
  sdreport_number last_slegalb;
  sdreport_number last_ssubllb;
  sdreport_number last_ofl
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
  log_ascr         -0.8
  log_bscr         4.6 
  log_aswr         -0.8
  log_bswr         4.6   
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
  sod = colsum(nonod)+colsum(noood);
  sot = colsum(nonot)+colsum(nooot);  
  scw = colsum(noncw)+colsum(noocw);
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
     onod(i) = elem_div(nonod(i),sod);
     oood(i) = elem_div(noood(i),sod); 
     onot(i) = elem_div(nonot(i),sot);
     ooot(i) = elem_div(nooot(i),sot); 	 
     oncw(i) = elem_div(noncw(i),scw);
     oocw(i) = elem_div(noocw(i),scw); 	 
     onsp(i) = elem_div(nonsp(i),ssp);
     oosp(i) = elem_div(noosp(i),ssp); 	
// Extract discards: newshell oldshell
     onod1(i) = onod(i)(1,6);
     oood1(i) = oood(i)(1,6);	 
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
   efsod = sod*efn(2);  
  for (i=1;i<=nyod;i++)
   {
    if (efsod(i) > maxs(2)) efsod(i)=maxs(2);
   }
// Calculate maximum effective smple sieze for Observer
   efsot = sot*efn(2);  
  for (i=1;i<=nyot;i++)
   {
    if (efsot(i) > maxs(2)) efsot(i)=maxs(2);
   }

// Calculate maximum effective smple sieze for Observer
   efcw = scw*efn(2);  
  for (i=1;i<=nycw;i++)
   {
    if (efcw(i) > maxs(2)) efcw(i)=maxs(2);
   }
      
   
   
//--------------------------------------------------------------------------
// 4.4  Calculate cv of cpue index      
//--------------------------------------------------------------------------  
  cvcpue = elem_div(secpue+1.e-3,stcpue+1.e-3);

//--------------------------------------------------------------------------
// 4.5  Calculate Winter subsistence discards (twsd)
//      For early years, total number of crab caught is not available. 
//      Estimate discards based on average proportion of discards   
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
// Trawl survey
  offset(1) = sum(elem_prod(efst,colsum(elem_prod(ont,log(ont+1.e-3)))));
  offset(2) = sum(elem_prod(efst,colsum(elem_prod(oot,log(oot+1.e-3)))));
// Winter survey  
  offset(3) = sum(elem_prod(efsw,colsum(elem_prod(onw,log(onw+1.e-3)))));
  offset(4) = sum(elem_prod(efsw,colsum(elem_prod(oow,log(oow+1.e-3)))));
// Summer Commercial Retain  
  offset(5) = sum(elem_prod(efsc,colsum(elem_prod(onc,log(onc+1.e-3)))));
  offset(6) = sum(elem_prod(efsc,colsum(elem_prod(ooc,log(ooc+1.e-3)))));
// Summer Commercial discards    
  offset(7) = sum(elem_prod(efsod,colsum(elem_prod(onod,log(onod+1.e-3)))));
  offset(8) = sum(elem_prod(efsod,colsum(elem_prod(oood,log(oood+1.e-3)))));
// Summer Commercial total  
  offset(9) = sum(elem_prod(efsot,colsum(elem_prod(onot,log(onot+1.e-3)))));
  offset(10) = sum(elem_prod(efsot,colsum(elem_prod(ooot,log(ooot+1.e-3))))); 
// Winter Commercial discards   
  offset(11) = sum(elem_prod(efcw,colsum(elem_prod(oncw,log(oncw+1.e-3)))));
  offset(12) = sum(elem_prod(efcw,colsum(elem_prod(oocw,log(oocw+1.e-3)))));
// Spring Survey    
  offset(13) = sum(elem_prod(efsp,colsum(elem_prod(onsp,log(onsp+1.e-3)))));
  offset(14) = sum(elem_prod(efsp,colsum(elem_prod(oosp,log(oosp+1.e-3)))));
// Summer Commercial discards larg vessel    
  offset(15) = sum(elem_prod(efsod(1,6),colsum(elem_prod(onod1,log(onod1+1.e-3)))));
  offset(16) = sum(elem_prod(efsod(1,6),colsum(elem_prod(oood1,log(oood1+1.e-3)))));    
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
  maximum_function_evaluations 40000
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
   
//   get_reference_points();
  
//==========================================================================
// 6.0  Function:  Convert parameters into rates      
//==========================================================================
FUNCTION convert_parameters_into_rates
  int i, j;
// working variables: molting probability
  dvariable mol,asw,bsw,aswr,bswr,alpha,beta;    
  dvar_vector ast(1,nst),bst(1,nst), asc(1,nsc),bsc(1,nsc),ascr(1,nscr),bscr(1,nscr);
// Initialize length proportions.   
  npp.initialize();				
  rpp.initialize();				
//catchability coefficient
   q = mfexp(log_q);
//   qw = mfexp(log_qw);   
// Summer Commercial fisheris selectivity parameters    
   asc = mfexp(log_asc);
   bsc = mfexp(log_bsc);  
// Summer Commercial fisheris retention parameters   
   ascr = mfexp(log_ascr);
   bscr = mfexp(log_bscr);
// Summer Trawl Survey selectivity parameters 
   ast = mfexp(log_ast);
   bst = mfexp(log_bst); 
// Winter Pot Survey selectivity parameters   
   asw = mfexp(log_asw);
   bsw = mfexp(log_bsw);
// Winter Pot retention parameters      
   aswr = mfexp(log_aswr);
   bswr = mfexp(log_bswr);
   
// Molting Parameters   
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
// Summer Commercial Catch selectivity 
 for (i=1;i<=nsc;i++)	 
   {
// 1 parameter logistic	   
  if(base(1)<=1) selc(i) = 1.0/(1.0+mfexp(asc(i)*(-mlen+mlen(na))+log(1.0/0.999-1.0))); 
// 2 parameters logistic	     
  if(base(1)>=2) selc(i) = 1.0/(1.0+mfexp(-asc(i)*(mlen - bsc(i)))); 
  selc(i,na) = 1.0;
   }   
   
// Summer Commercial Retention selectivity
  if(base(2)>0) 
  {
   for (i=1;i<=nscr;i++)	 
   {   
   selr(i) = 1.0/(1.0+mfexp(-ascr(i)*(mlen - bscr(i)))); 
   selr(i,na) = 1.0;
   }   
  }
  
// Winter Pot Survey Selectivity  	
   selw = 1.0/(1.0+mfexp(asw*(mlen-bsw)));  
// Directly estimate small length class 
   for (j=1;j<=msn(2);j++) selw(j) = sw3(j);
   selwr = 1.0/(1.0+mfexp(-aswr*(mlen - bswr))); 
   selwr(na) = 1.0;   
  
// Winter commercial retained and discarded 
  if(base(3)>0)
	{
		lselw = elem_prod(selw,selwr);   
		dlselw = elem_prod(selw,1-selwr);  
	}
	else
    {
		lselw = elem_prod(selw,lg);   
		dlselw = elem_prod(selw,slg);  
    }  

      
// Summer Trawl survey Selectivity 
// NOAA summer trawl net selectivity 1976-1992	 
//     selt(1) = 1.0/(1.0+mfexp(-ast(1)*(mlen-bst(1))));   
// ADFG summer trawl net selectivity: 1996- 	 
//     selt(2) = 1.0/(1.0+mfexp(-ast(nst)*(mlen-bst(nst))));   
   selt(1) = 1.0/(1.0+mfexp(ast(1)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));   
   selt(2) = 1.0/(1.0+mfexp(ast(nst)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));        
   selt(1,na) = 1.0;
   selt(2,na) = 1.0;
   

// Define natural mortality  
// Set mortality be 1.0 for some length crab
  for (j=1;j<=(na-msn(1));j++) ms1(j) = 1.0;
    if(base(6)>=1) ms1(7)=ms1(8);    
	Mn = M*ms1;
	    
// Assign Trawl survey seletivity function in different time periods. 
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(1))    // period 1976 - 1992 
      selty(i) = selt(1);   //NOAA trawl survey 
      else                          
      selty(i) = selt(2);   //ADFG trawl survey
    }

// Assign Summer commrecial catch seletivity function in different time periods. 
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(2))  		 // period 1976 - 1992
	 {
        selcy(i) = selc(1); 
		selry(i) = selr(1);
	 }
	  else if(i<=scy(3))            // period 1993 - 2005
	  { 
	    selcy(i) = selc(1);
		selry(i) = selr(1);
	    if(nsc>=2) selcy(i) = selc(2);
		if(nscr>=2) selry(i) = selr(2);
	  }
	  else            // period 2006 - present
	  {
        selcy(i) = selc(1);
		selry(i) = selr(1);
       if(nsc>=2) selcy(i) = selc(2);
	   if(nscr>=2) selry(i) = selr(2);
	   if(nsc>=3)  selcy(i) = selc(3);
	   if(nscr>=3) selry(i) = selr(3);		
	  }
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

// Assign first year abundance: 
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
    egr1(i) = elem_prod(mgr1(i),selc(1));
    egr2(i) = elem_prod(mgr2(i),selc(1));
    egr3(i) = elem_prod(mgr3(i),selc(1));
// Selectivity for 1993-present   
    egr12(i) = elem_prod(mgr1(i),selc(1));
    egr22(i) = elem_prod(mgr2(i),selc(1));
    egr32(i) = elem_prod(mgr3(i),selc(1));    
  if(nsc>=2)
  {
    egr12(i) = elem_prod(mgr1(i),selc(2));
    egr22(i) = elem_prod(mgr2(i),selc(2));
    egr32(i) = elem_prod(mgr3(i),selc(2));    	  
  }	  
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
  dvariable pp, wpm, wpim, tt1, ttt0, ttt1, sb, tw, ews,ewsd;
  dvar_vector nscaf(1,na), tsc(1,na),tscd(1,na), tt0(1,na);    // Abundance after summer fishery
  dvar_vector sselw(1,na),dselw(1,na),sselc(1,na),dselc(1,na);  
  dvar_vector nwc(1,na),owc(1,na),nwcd(1,na),owcd(1,na);
  dvar_vector nws(1,na),ows(1,na),nwsd(1,na),owsd(1,na); 
  dvar_vector par1(1,na),par2(1,na), FOFL(1,na);   
  dvariable TotalCatch, CatNum, cofl,fl;    
  enc.initialize();
  eoc.initialize();
  enod.initialize();
  eood.initialize(); 
  enod1.initialize();
  eood1.initialize();   
  enot.initialize();
  eoot.initialize();   
  tt0.initialize();
  last_ofl.initialize();
  last_subofl.initialize();
// 
// Calculate winter subsistence harvest selectivity
// Assume subsistence harvets take all mature crab caught
   sselw = elem_prod(matc,selw); 
// Assume subsistence harvest discards all immature crab caught   
   dselw = elem_prod(imatc,selw); 
   

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
			
// Calculate Total Winter Legal Crab retainable to commercial fisheries 
      tw = sum(elem_prod(npw(i)+opw(i),lselw));	  

// Calculate number of newshell and oldshell retained in commercial catch   	  
      nwc = (twc(i)/tw)*elem_prod(npw(i),lselw);
      owc = (twc(i)/tw)*elem_prod(opw(i),lselw);
// Calculate number of newshell and oldshell discards in commercial  	  
      nwcd = (twc(i)/tw)*elem_prod(npw(i),dlselw);
      owcd = (twc(i)/tw)*elem_prod(opw(i),dlselw);
// Calculate total number of winter commercial discards mortality
	  bcw(i) = sum(nwcd+owcd)*hm(2); 	  
// Calculate total Winter Crab (Leal + Sublegal) Catchable to commercial
//      twp(i) = sum(elem_prod(npw(i)+opw(i),selw)); 
   
// Calculate exploitable abundance at the middle of the season
//      twp(i) = twp(i) - 0.5*(twc(i)+tws(i));  

//==========================================================================
// 8.2 Winter subsistence catch and discards size composition 
//     Assume that Subsistence fishers took all length >= 3 crab  
//	   and discarded length 1 & 2 crab
//==========================================================================

// Calculate total mature population  to winter subsistence  
     ews = sum(elem_prod(npw(i)+opw(i),sselw));
// Calculate total immature population  to winter subsistence  
     ewsd = sum(elem_prod(npw(i)+opw(i),dselw));
	 
// Calculate number of newshell and oldshell in winter subsistence 
	 nws = (tws(i)/ews)*elem_prod(npw(i),sselw);
     ows = (tws(i)/ews)*elem_prod(opw(i),sselw);

// Calculate number of new and old shell discards in winter subsistence  
	 nwsd = (twsd(i)/ewsd)*elem_prod(npw(i),dselw);
     owsd = (twsd(i)/ewsd)*elem_prod(opw(i),dselw); 	
// Calculate total number of winter commercial discards mortality
	 bswd(i) = sum(nwsd+owsd)*hm(2); 	
// Calculate Subsistence harvest rate
    FSub(i) = tws(i)/ews;
    FSubd(i) = twsd(i)/ewsd;	
	 
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
  if(base(2)<0)
     {
		sselc = elem_prod(selcy(i),lg);   
		dselc = elem_prod(selcy(i),slg);  
     }  
    else
	{
		sselc = elem_prod(selcy(i),selry(i));   
		dselc = elem_prod(selcy(i),1-selry(i));  
	}

// Calculate summer crab legal and sublegal crab selectivity
   
 
// Total number of crab retainable to commercial fishery    
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
// Expected discards = (potential discards)*(actual retain)/(potential retain) 	
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
   }

 for (i=fyear; i<=lyear; i++)
   {
//    mmb(i) = sum(elem_prod(elem_prod(npw(i)+opw(i),matc),wm));  // Mature male biomass
   if(base(2)<0)
   {	   
    legaln(i) = sum(elem_prod(npw(i)+opw(i),lg));  // number of legal crab
    legalb(i) = sum(elem_prod(npw(i)+opw(i),wmlg)); // legal biomass   
	}
     else
   {	   
    legaln(i) = sum(elem_prod(npw(i)+opw(i),selry(i)));  // number of legal crab
    legalb(i) = sum(elem_prod(elem_prod(npw(i)+opw(i),selry(i)),wm)); // legal biomass   
	}
   }
   
  last_y = npw(lyear+1)+opw(lyear+1);
  
  if(base(2)<0)
  {	  
  last_legalb = elem_prod(elem_prod(elem_prod(last_y,selcy(lyear)),lg),wm);
  last_sublb = elem_prod(elem_prod(elem_prod(last_y,selcy(lyear)),slg),wm);
  }
    else 
  {	  
  last_legalb = elem_prod(elem_prod(elem_prod(last_y,selcy(lyear)),selry(lyear)),wm);
  last_sublb = elem_prod(elem_prod(elem_prod(last_y,selcy(lyear)),1-selry(lyear)),wm);
  }
		
  last_mmb= mmb(lyear+1);
  bmsy = (sum(mmb) - mmb(1976)-mmb(1977)-mmb(1978)-mmb(1979))/(lyear-fyear-2);
  cofl = last_mmb/bmsy;
  if (cofl < 0.25){
	  fl = 0;  
        }
	   else if (cofl <=1){ 
	   fl = (cofl-0.1)/0.9;
	   }
	   else{
		   fl = 1;
	}
// Calculate FOFL	
  FOFL = fl*M;
  par1 = 1-exp(-FOFL-0.42*Mn);
  par2 = 1-exp(-0.42*Mn);
  last_ofl = sum(elem_prod(last_legalb,par1-elem_prod(par2,elem_div(1-pwh*par1,1-pwh*par2))))/1000; 
  last_subofl = sum(elem_prod(last_sublb,par1-elem_prod(par2,elem_div(1-pwh*par1,1-pwh*par2))))/1000; 
  last_slegalb = sum(last_legalb)/1000;
  last_ssubllb = sum(last_sublb)/1000;
  
//==========================================================================
// 9.0  produce model estimated observation data  
//==========================================================================
FUNCTION get_proportion_and_effort
  int i,j;
  dvariable bf,af,pp,ts;
  dvar_vector ent0(1,na), eot0(1,na);
  ett.initialize();
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
    ent0 = elem_prod(elem_prod((elem_prod(nps(it(i)),exp(-bf*Mn))-pct(i)*tc(it(i))*enc(it(i))),exp(-af*Mn)),selty(it(i)));
    eot0 = elem_prod(elem_prod((elem_prod(ops(it(i)),exp(-bf*Mn))-pct(i)*tc(it(i))*eoc(it(i))),exp(-af*Mn)),selty(it(i)));
// Stop gap measure	 
	for (j=1;j<=na;j++)
     {
	  if (ent0(j) < 0.0) ent0(j) = 0.0;
	  if (eot0(j) < 0.0) eot0(j) = 0.0;	  
     }
// Calculate trawl survey abundance 
	 ett(i) = sum(ent0+eot0);
// Stop gap measuer
    if (ett(i) <= 0.0) ett(i) = 0.00001;

// Calculate proprion of new and old shell crab by trawl survey	
    ent(i) = ent0/ett(i);
	eot(i) = eot0/ett(i);
 
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
// 9.3  Predicted Winter Commercial length proportion 
//==========================================================================
  for (i=1;i<=nycw;i++)
   {
// Calculate total winter pot catchable 
    pp = sum(elem_prod(npw(icw(i))+opw(icw(i)),lselw));
// Calculate proportion of pot survey shell length  
	encw(i) = elem_prod(npw(icw(i)),lselw)/pp;
	eocw(i) = elem_prod(opw(icw(i)),lselw)/pp;
   }
   
//==========================================================================
// 9.2  Predicted Spring Pot survey length proportion 
//==========================================================================
  for (i=1;i<=nysp;i++)
   {   
// Calculate total winter pot catchable 
    pp = sum(elem_prod(nps(isp(i))+nps(isp(i)),selw));
    if (pp <= 0.0) pp = 0.000001;
// Calculate proportion of pot survey shell length  
	ensp(i) = elem_prod(nps(isp(i)),selw)/pp;
	eosp(i) = elem_prod(ops(isp(i)),selw)/pp;
   }
 

//==========================================================================
// 9.3  Predicted Observer survey length proportion  
//==========================================================================
  for (i=1;i<=nyod;i++)
   {
// Calculate total summer com sublegal   
	pp = sum(elem_prod(elem_prod(nps(iod(i))+ops(iod(i)),selcy(iod(i))),slg));
// Calculate proportion of sublegal shell, length 
    enod(i) = elem_prod(elem_prod(nps(iod(i)),selcy(iod(i))),slg)/pp;
    eood(i) = elem_prod(elem_prod(ops(iod(i)),selcy(iod(i))),slg)/pp;
   }
  for (i=1;i<=6;i++)
   {
// Calculate total summer com sublegal   
	pp = sum(elem_prod(elem_prod(nps(iod(i))+ops(iod(i)),selcy(iod(i))),1-selry(iod(i))));
// Calculate proportion of sublegal shell, length 
    enod1(i) = elem_prod(elem_prod(nps(iod(i)),selcy(iod(i))),1-selry(iod(i)))/pp;
    eood1(i) = elem_prod(elem_prod(ops(iod(i)),selcy(iod(i))),1-selry(iod(i)))/pp;
   }

  for (i=1;i<=nyot;i++)
   {
// Calculate total summer com sublegal   
	pp = sum(elem_prod(nps(iot(i))+ops(iot(i)),selcy(iot(i))));
// Calculate proportion of sublegal shell, length 
    enot(i) = elem_prod(nps(iot(i)),selcy(iot(i)))/pp;
    eoot(i) = elem_prod(ops(iot(i)),selcy(iot(i)))/pp;
   }
   
   
   
//==========================================================================
// 9.4  Predicted CPUE
//==========================================================================

  for (i=fyear;i<=lyear;i++)
  {
    ts = tb(i) - 0.5*tc(i);  //exploitable abundance at the middle of the season
    if (i <= scy(1))
	 {
	 ecpue(i) = q(1)*ts;
	 }
    else
	 {	
	 ecpue(i) = q(2)*ts;
	 }
// if stcpue is zero (no fishery)  ecpue = 0;
	if (stcpue(i) <= 0.0) ecpue(i) = 0.0;
  }
  

//==========================================================================
// 10.0  Likelihood Calculation  
//==========================================================================
FUNCTION evaluate_the_objective_function

// Log likelihood for trawl survey 
  tf(1) = 0.5*norm2(elem_div((log(tt+1.e-3)-log(ettq+1.e-3)),sqrt(log(elem_prod(cv,cv)+1.0))));

// Log likelihood standard cpue in winter fishery 
//  tf(2) = norm2(log(wcpue+1.e-3)-log(ewcpue+1.e-3))/(2*SD(2)*SD(2)); 
   
// Log likelihood standard cpue in summer fishery 
  T_var = sqrt(log(elem_prod(cvcpue,cvcpue)+1.0)+advar); 
  tf(3) = (sum(log(T_var))+0.5*(norm2(elem_div((log(stcpue+1.e-3)-log(ecpue+1.e-3)),T_var))));  

//Log likelhihood for trawl survey multinomial proportion     
  tf(4) = -(sum(elem_prod(efst,rowsum(elem_prod(trans(ont),log(ent+1.e-3))))) - offset(1))      
          -(sum(elem_prod(efst,rowsum(elem_prod(trans(oot),log(eot+1.e-3))))) - offset(2));   
   
//Log likelhihood for winter multinomial proportion   
  tf(5) = -(sum(elem_prod(efsw,rowsum(elem_prod(trans(onw),log(enw+1.e-3))))) - offset(3))     
          -(sum(elem_prod(efsw,rowsum(elem_prod(trans(oow),log(eow+1.e-3))))) - offset(4));       
 
//Log likelhihood size proportion for summer fishery survey    
  tf(6) = -(sum(elem_prod(efsc,rowsum(elem_prod(trans(onc),log(enc+1.e-3))))) - offset(5))         
          -(sum(elem_prod(efsc,rowsum(elem_prod(trans(ooc),log(eoc+1.e-3))))) - offset(6));      
		  
  if(base(2) < 0){  
//Log likelhihood size proportion for observer discard    
  tf(7) = -(sum(elem_prod(efsod,rowsum(elem_prod(trans(onod),log(enod+1.e-3))))) - offset(7)) 
          -(sum(elem_prod(efsod,rowsum(elem_prod(trans(oood),log(eood+1.e-3))))) - offset(8));  
  }
  else
  {
//Log likelhihood size proportion for observer total catch   
  tf(7) = -(sum(elem_prod(efsot,rowsum(elem_prod(trans(onot),log(enot+1.e-3))))) - offset(9)) 
          -(sum(elem_prod(efsot,rowsum(elem_prod(trans(ooot),log(eoot+1.e-3))))) - offset(10)) 
         -likew(3)*(sum(elem_prod(efsod(1,6),rowsum(elem_prod(trans(onod1),log(enod1+1.e-3))))) - offset(15)) 
         -likew(3)*(sum(elem_prod(efsod(1,6),rowsum(elem_prod(trans(oood1),log(eood1+1.e-3))))) - offset(16));	  
  }

//Log likelhihood for Winter Commercial Retained multinomial proportion   
  if(base(3)>0){
  tf(8) = -(sum(elem_prod(efcw,rowsum(elem_prod(trans(oncw),log(encw+1.e-3))))) - offset(11));     
          -(sum(elem_prod(efcw,rowsum(elem_prod(trans(oocw),log(eocw+1.e-3))))) - offset(12));       
	}	
  
//Log likelhihood for Spring survey  multinomial proportion   
  if(likew(2)>0){
  tf(9) = -(sum(elem_prod(efsp,rowsum(elem_prod(trans(onsp),log(ensp+1.e-3))))) - offset(13));     
          -(sum(elem_prod(efsp,rowsum(elem_prod(trans(oosp),log(eosp+1.e-3))))) - offset(14));       
	}	
	
//deviation in recruits.  
  tf(10) = norm2(log_relrec)/(2*SD(1)*SD(1));                            
//deviation in molting parameter.  
  tf(11) = norm2(log_amol_dev)/(2*SD(1)*SD(1));     
//deviation in molting parameter.  
  tf(12) = norm2(log_bmol_dev)/(2*SD(1)*SD(1));     

  
// tag recovery likelihood
  
  ef(1) = -(sum(elem_prod(rowsum(tag1),rowsum(elem_prod(ptag1,log(egr1+1.e-3)))))- toffset(1)); 
  ef(2) = -(sum(elem_prod(rowsum(tag2),rowsum(elem_prod(ptag2,log(egr2+1.e-3)))))-toffset(2));  
  ef(3) = -(sum(elem_prod(rowsum(tag3),rowsum(elem_prod(ptag3,log(egr3+1.e-3)))))-toffset(3));  
  ef(4) = -(nsc-1)*(sum(elem_prod(rowsum(tag12),rowsum(elem_prod(ptag12,log(egr12+1.e-3)))))-toffset(4)); 
  ef(5) = -(nsc-1)*(sum(elem_prod(rowsum(tag22),rowsum(elem_prod(ptag22,log(egr22+1.e-3)))))-toffset(5));  
  ef(6) = -(nsc-1)*(sum(elem_prod(rowsum(tag32),rowsum(elem_prod(ptag32,log(egr32+1.e-3)))))-toffset(6));  
  
  tf(13) = likew(1)*sum(ef);
  
// Total Likelihood
  f = sum(tf);   

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
  dvar_vector sselw(1,na),dselw(1,na);
  dvar_vector sselc(1,na),dselc(1,na); 
  dvar_vector nwc(1,na),owc(1,na),nwcd(1,na),owcd(1,na);
  dvar_vector nws(1,na),ows(1,na),nwsd(1,na),owsd(1,na);
  dvar_vector nsc(1,na),osc(1,na),nscd(1,na),oscd(1,na);  
  dvar_vector par1(1,na),par2(1,na),fx(1,na), ref_Hs(1,na),ref_Hw(1,na);   
  dvariable TotalCatch, CatNum, cofl,fl;    
  dvar_matrix ref_ws(1,na,1,4);  
  dvariable ref_FSub, ref_FSubd, ref_rec;
  dvariable ref_Hwc;
  dvariable ref_Hts;
  dvar_vector ref_Htw(1,301);
  dvar_vector ref_mbio(1,301);
  dvar_vector ref_mmb(1,100);
  dvar_vector ref_fcatch_s(1,301);
  dvar_vector ref_fcatch_sd(1,301);
  dvar_vector ref_fcatch_w_c(1,301);
  dvar_vector ref_fcatch_w_s(1,301);
  dvar_vector ref_fcatch_w_cd(1,301);
  dvar_vector ref_fcatch_w_sd(1,301);   
  dvar_vector ref_legaln(1,301);
//  dvar_vector ref_Hrt(1,100);
//  dvar_vector ref_Hr(1,301);
  dvar_vector ref_legal(1,100);
  dvar_vector ref_legal_b(1,100);
  dvar_vector ref_F(1,301);
  dvariable ref_ys;
  dvariable f35,f40,h35,h40,i35,i40;

  ref_Htw.initialize();
  ref_FSub.initialize();
  ref_FSubd.initialize();  
  
  ref_ys = ys(lyear);  // use mid-point of comfish as reference date
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
  
// Start lyear as starting year   
   ref_ys = ys(lyear);
// Calculate winter subsistence harvest selectivity
// Assume subsistence harvest take all mature crab caught
   sselw = elem_prod(matc,selw); 
// Assume subsistence harvest discards all immature crab caught   
   dselw = elem_prod(imatc,selw); 
// Calculate legal crab and sublegal crab caught by Summer Commercial   
  if(base(2)<0)
     {
		sselc = elem_prod(selcy(lyear),lg);   
		dselc = elem_prod(selcy(lyear),slg);  
     }  
    else
	{
		sselc = elem_prod(selcy(lyear),selry(lyear));   
		dselc = elem_prod(selcy(lyear),1-selry(lyear));  
	}
// Use the lyear's retcurit as reference constant; 
   ref_rec = rec(lyear);


//--------------------------------------------------------------------------
// Simulation for 300 Fishery F itelation    
//--------------------------------------------------------------------------
  ref_Hts = 0.0;
  for (kk = 1; kk <= 301; kk++)
   {
// Set Fishery mortality
  ref_F(kk) = 0.01*kk-0.01;

//--------------------------------------------------------------------------
// 8.1 Based on OFL approach Calculate Winter and Summer commercial Harvest rate    
//--------------------------------------------------------------------------	
  par1 = 1-exp(-ref_F(kk)-0.42*Mn);
  par2 = 1-pwh*(1-exp(-0.42*Mn));
//Winter commercial harvest rate   
  ref_Hw = elem_div(pwh*(1-exp(-ref_F(kk)))*exp(-0.42*Mn),par2);
// Summer commercial harvest rate  
  ref_Hs = elem_div((1-pwh)*(1-exp(-ref_F(kk)))*exp(-0.42*Mn),par2); 
// Initialize temp matrix
  ref_nps.initialize();
  ref_ops.initialize();
  ref_npw.initialize();
  ref_opw.initialize();
// Set First year is winter abundance 
  ref_npw(1) = npw(lyear+1);
  ref_opw(1) = opw(lyear+1);

//==========================================================================
// 8.7  Start 100 year Simulation that reaches to equilibrium
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
//     Assume that Subsistence fishers catch does not change 
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
  ref_nps(i) = elem_prod(ref_npw(i)-nwc-nws-(nwcd+nwsd)*hm(2),exp(-0.417*Mn));  
  ref_ops(i) = elem_prod(ref_opw(i)-owc-ows-(owcd+owsd)*hm(2),exp(-0.417*Mn));  
		 
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
//	ref_Hrt(i) = sum(ref_catch_s(i)+ref_catch_w_c(i))/sum(elem_prod(ref_npw(i)+ref_opw(i),lg));  
//	ref_Hrt(i) = sum(ref_catch_w_c(i))/sum(ref_catch_w_c(i)+ref_catch_s(i));  
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
  ref_legal(i) = sum(elem_prod(elem_prod(ref_npw(i)+ref_opw(i),sselc),wm));  // legal biomass crab  
   }	  

//==========================================================================
// Extract Harvest and biomass of the equilibrium
//==========================================================================
    ref_mbio(kk) = ref_mmb(99);                    
//    ref_totc(j) = ref_catch_b(100);
    ref_legaln(kk) = ref_legal(99);
	ref_fcatch_s(kk) = sum(ref_catch_s(99));
    ref_fcatch_sd(kk) = sum(ref_catch_sd(99));
    ref_fcatch_w_c(kk) = sum(ref_catch_w_c(99));
    ref_fcatch_w_s(kk) = sum(ref_catch_w_s(99));
    ref_fcatch_w_cd(kk) = sum(ref_catch_w_cd(99));
    ref_fcatch_w_sd(kk) = sum(ref_catch_w_sd(99));  
//    cout << "ncrab" << ref_opw(99)+ref_npw(99)<<endl; 
  }

//==========================================================================
// 9.0  Find F35%
//==========================================================================
   
   i35 = 0; i40 = 0;
  
  for (j = 1; j<= 301; j++)
   {
    if (i35 < 1.0)
// Find F35 F35 = 35% of reference mmb		
     if (ref_mbio(j) <= 0.35*ref_mbio(1))
      {
       f35 = ref_F(j);
       i35 = ref_mbio(j);
      }
    if (i40 < 1.0)
     if (ref_mbio(j) <= 0.40*ref_mbio(1))
      {
       f40 = ref_F(j);
       i40 = ref_mbio(j);
      }
   }

  ofstream report1("refp.out");
  report1 <<"F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_F<<endl;
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
  report1 <<"F35: "<<f35<<"  B35%: "<<i35<<endl;
  report1 <<"F40: "<<f40<<"  B40%: "<<i40<<endl;
//  report1 << ref_npw<<endl;    
//  report1 << ref_opw<<endl;    
//  report1 <<"ref_Hwc = (mean of 10 years):  "<<ref_Hwc<<endl;
  cout<<"end of reference points"<<endl;

//==========================================================================

  
//==========================================================================
REPORT_SECTION
  cout << "Report Section" << endl;
  get_reference_points();
  
  report << "nps" << endl << trans(nps) << endl;  // Modeled new shell summer 
  report << "ops" << endl << trans(ops) << endl;  // Modeled old shell summer
  report << "npw" << endl << trans(npw) << endl;  // Modeled new shell winter
  report << "opw" << endl << trans(opw) << endl;  // Modeled oldshell shell winter  
  report << "ettq" << endl << ettq << endl;  // Estimated trawl abundance
  report << "ecpue" << endl << ecpue << endl;  // Estimated Summer fishery cpue
//  report << "ewcpue" << endl << ewcpue << endl;  // Estimated Winter survey cpue
  report << "ent" << endl << trans(ent) << endl;  // Estimated trawl newshell size proportion 
  report << "eot" << endl << trans(eot) << endl;  // Estimated trawl oldshell size proportion   
  report << "enw" << endl << trans(enw) << endl;  // Estimated winter survey newshell size proportion
  report << "eow" << endl << trans(eow) << endl;  // Estimated winter survey oldshell size proportion
  report << "enc" << endl << trans(enc) << endl;  // Estimated summer fishery newshell size proportion
  report << "eoc" << endl << trans(eoc) << endl;  // Estimated summer fishery oldshell size proportion
  report << "enod" << endl << trans(enod) << endl;  // Estimated observer newshell size proportion
  report << "eood" << endl << trans(eood) << endl;  // Estimated observer oldshell size proportion
  report << "enod1" << endl << trans(enod1) << endl;  // Estimated observer newshell size proportion
  report << "eood1" << endl << trans(eood1) << endl;  // Estimated observer oldshell size proportion
  report << "enot" << endl << trans(enot) << endl;  // Estimated observer newshell size proportion
  report << "eoot" << endl << trans(eoot) << endl;  // Estimated observer oldshell size proportion  
  report << "encw" << endl << trans(encw) << endl;  // Estimated spring Pot survey newshell size proportion
  report << "eocw" << endl << trans(eocw) << endl;  // Estimated spring Pot survey oldshell size proportion
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
  report << "selc" << endl << selc << endl;
  report << "selr" << endl << selr << endl;  
  report << "selt" << endl << selt << endl; 
  report << "selw" << endl << selw << endl;
  report << "selwr" << endl << selwr << endl;  
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
  report << "selty" << endl << trans(selty) <<endl; 
  report << "selcy" << endl << trans(selcy) <<endl; 
  report << "selry" << endl << trans(selry) <<endl; 
  
 
//==========================================================================
GLOBALS_SECTION
  #include <math.h>
  #include <admodel.h>
  #include <time.h>
  ofstream mcmc1,mcmc2;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
  int phz;    // phase

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