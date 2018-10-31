  #include <math.h>
  #include <admodel.h>
  #include <time.h>
  ofstream mcmc1,mcmc2;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <NS3n2013cpue_rev.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  ny.allocate("ny");
  na.allocate("na");
  nyt.allocate("nyt");
  nyp.allocate("nyp");
  nyw.allocate("nyw");
  nyf.allocate("nyf");
  nyo.allocate("nyo");
  scy1.allocate("scy1");
  scy2.allocate("scy2");
  scp.allocate("scp");
  slm.allocate("slm");
  slt.allocate("slt");
  M.allocate("M");
  ms.allocate("ms");
  slt6.allocate("slt6");
  qta.allocate("qta");
  spcv.allocate("spcv");
  lamt.allocate("lamt");
  lamp.allocate("lamp");
  lamc.allocate("lamc");
  maxss.allocate("maxss");
  efn.allocate("efn");
  lg.allocate(1,na,"lg");
  gr.allocate(1,na,1,na,"gr");
  it.allocate(1,nyt,"it");
  tt.allocate(1,nyt,"tt");
  pct.allocate(1,nyt,"pct");
  yt.allocate(1,nyt,"yt");
  st.allocate(1,nyt,"st");
  cv1.allocate(1,nyt,"cv1");
  cv2.allocate(1,nyt,"cv2");
  cv3.allocate(1,nyt,"cv3");
  ip.allocate(1,nyp,"ip");
  tp.allocate(1,nyp,"tp");
  yp.allocate(1,nyp,"yp");
  sp.allocate(1,nyp,"sp");
  iw.allocate(1,nyw,"iw");
  sw.allocate(1,nyw,"sw");
  ipre.allocate("ipre");
  spre.allocate("spre");
  sc.allocate(1,ny,"sc");
  tc.allocate(1,ny,"tc");
  te.allocate(1,ny,"te");
  stcpue.allocate(1,ny,"stcpue");
  secpue.allocate(1,ny,"secpue");
  ys.allocate(1,ny,"ys");
  io.allocate(1,nyo,"io");
  so.allocate(1,nyo,"so");
  twc.allocate(1,ny-1,"twc");
  tsc.allocate(1,ny-1,"tsc");
  ont.allocate(1,na,1,nyt,"ont");
  oot.allocate(1,na,1,nyt,"oot");
  onp.allocate(1,na,1,nyp,"onp");
  oop.allocate(1,na,1,nyp,"oop");
  onw.allocate(1,na,1,nyw,"onw");
  oow.allocate(1,na,1,nyw,"oow");
  onf.allocate(1,na,"onf");
  oof.allocate(1,na,"oof");
  onc.allocate(1,na,1,ny,"onc");
  ooc.allocate(1,na,1,ny,"ooc");
  ono.allocate(1,na,1,nyo,"ono");
  ooo.allocate(1,na,1,nyo,"ooo");
  wm.allocate(1,na,"wm");
  hm.allocate("hm");
  log_initp.allocate("log_initp");
  initp_phase.allocate("initp_phase");
  wsp.allocate(1,2,"wsp");
  qtno_phase.allocate("qtno_phase");
  qtad_phase.allocate("qtad_phase");
  M_phase.allocate("M_phase");
  ms_phase.allocate("ms_phase");
  tt1.allocate(1,nyt);
  tt2.allocate(1,nyt);
  tp1.allocate(1,nyp);
  tp2.allocate(1,nyp);
 ad_comm::change_datafile_name("proj.ctl");
  ProjType.allocate("ProjType");
  FirstYr.allocate("FirstYr");
  NprojFmsy.allocate("NprojFmsy");
  NprojSim.allocate("NprojSim");
 Nproj = NprojFmsy;
 if (NprojSim > Nproj) Nproj = NprojSim;
  alpha.allocate("alpha");
  beta.allocate("beta");
  Bmsy_yr1.allocate("Bmsy_yr1");
  Bmsy_yr2.allocate("Bmsy_yr2");
  CVEstWithin.allocate("CVEstWithin");
  CVEstExtra.allocate("CVEstExtra");
  Buffer.allocate("Buffer");
  BurnIn.allocate("BurnIn");
  LastSim.allocate("LastSim");
 LastSim = LastSim + BurnIn;
  SR_FORM.allocate("SR_FORM");
  Steep.allocate("Steep");
  SigmaR.allocate("SigmaR");
  Lag.allocate("Lag");
  PrintAll.allocate("PrintAll");
 mcmc_count = 0;
 fn_count = 0;
 MaxSim = 1000;
  RecInd.allocate(1,MaxSim,ny,ny+Nproj);
 cout << "Data Section Completed" << endl;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_q1.allocate(-32.5,8.5,1,"log_q1");
  log_q2.allocate(-32.5,10.0,1,"log_q2");
  log_q3.allocate(-32.5,10.0,-1,"log_q3");
  log_initpop.allocate(2,15,initp_phase,"log_initpop");
  log_recscale.allocate(2,12,2,"log_recscale");
  log_relrec.allocate(1,ny-2,-20.,20.,2,"log_relrec");
  flnp.allocate(1,na-1,-5.,5.,1,"flnp");
  flop.allocate(1,na-2,-5.,5.,-1,"flop");
  nop.allocate(0.5,0.9,1,"nop");
  r1.allocate(0.5,0.9,1,"r1");
  log_mo1.allocate(-5.5,-2.0,1,"log_mo1");
  log_mo2.allocate(0.55,10.0,1,"log_mo2");
  log_st1.allocate(-10.0,-1.0,1,"log_st1");
  log_st2.allocate(0.51,10.0,1,"log_st2");
  log_sp1.allocate(-3.0,-1.00,1,"log_sp1");
  log_sp2.allocate(3.9,5.50,1,"log_sp2");
  log_sw1.allocate(-10.0,10.00,1,"log_sw1");
  log_sw2.allocate(3.90,5.50,1,"log_sw2");
  sw3.allocate(0.1,1.0,1,"sw3");
  log_sc1.allocate(-3.0,-1.00,1,"log_sc1");
  log_sc2.allocate(3.9,5.5,1,"log_sc2");
  log_sc3.allocate(-3.0,-1.00,1,"log_sc3");
  log_sc4.allocate(3.9,5.5,1,"log_sc4");
  log_sc5.allocate(-3.0,-1.00,-1,"log_sc5");
  log_sc6.allocate(3.9,5.50,-1,"log_sc6");
  sos0.allocate(0.50,1.50,-1,"sos0");
  sow.allocate(0.300,1.70,-1,"sow");
  advar.allocate(0,6.0,2,"advar");
  qtno.allocate(0.1,2.5,qtno_phase,"qtno");
  qtad.allocate(0.1,2.5,qtad_phase,"qtad");
  M.allocate(0.1,0.3,M_phase,"M");
  ms.allocate(1.0,5.0,ms_phase,"ms");
  p4.allocate(0.1,1,1,"p4");
  log_rec.allocate(1,ny-2,"log_rec");
  #ifndef NO_AD_INITIALIZE
    log_rec.initialize();
  #endif
  nps.allocate(1,na,1,ny+Nproj,"nps");
  #ifndef NO_AD_INITIALIZE
    nps.initialize();
  #endif
  ops.allocate(1,na,1,ny+Nproj,"ops");
  #ifndef NO_AD_INITIALIZE
    ops.initialize();
  #endif
  npw.allocate(1,na,1,ny+Nproj,"npw");
  #ifndef NO_AD_INITIALIZE
    npw.initialize();
  #endif
  opw.allocate(1,na,1,ny+Nproj,"opw");
  #ifndef NO_AD_INITIALIZE
    opw.initialize();
  #endif
  rec.allocate(1,ny+Nproj,"rec");
  #ifndef NO_AD_INITIALIZE
    rec.initialize();
  #endif
  ent.allocate(1,na,1,nyt,"ent");
  #ifndef NO_AD_INITIALIZE
    ent.initialize();
  #endif
  et0.allocate(1,na,1,nyt,"et0");
  #ifndef NO_AD_INITIALIZE
    et0.initialize();
  #endif
  enp.allocate(1,na,1,nyp,"enp");
  #ifndef NO_AD_INITIALIZE
    enp.initialize();
  #endif
  eop.allocate(1,na,1,nyp,"eop");
  #ifndef NO_AD_INITIALIZE
    eop.initialize();
  #endif
  enw.allocate(1,na,1,nyw,"enw");
  #ifndef NO_AD_INITIALIZE
    enw.initialize();
  #endif
  eow.allocate(1,na,1,nyw,"eow");
  #ifndef NO_AD_INITIALIZE
    eow.initialize();
  #endif
  enf.allocate(1,na,"enf");
  #ifndef NO_AD_INITIALIZE
    enf.initialize();
  #endif
  eof.allocate(1,na,"eof");
  #ifndef NO_AD_INITIALIZE
    eof.initialize();
  #endif
  enc.allocate(1,na,1,ny,"enc");
  #ifndef NO_AD_INITIALIZE
    enc.initialize();
  #endif
  eoc.allocate(1,na,1,ny,"eoc");
  #ifndef NO_AD_INITIALIZE
    eoc.initialize();
  #endif
  eno.allocate(1,na,1,nyo,"eno");
  #ifndef NO_AD_INITIALIZE
    eno.initialize();
  #endif
  eoo.allocate(1,na,1,nyo,"eoo");
  #ifndef NO_AD_INITIALIZE
    eoo.initialize();
  #endif
  en0.allocate(1,na,1,ny-1,"en0");
  #ifndef NO_AD_INITIALIZE
    en0.initialize();
  #endif
  eo0.allocate(1,na,1,ny-1,"eo0");
  #ifndef NO_AD_INITIALIZE
    eo0.initialize();
  #endif
  en1.allocate(1,na,1,ny-1,"en1");
  #ifndef NO_AD_INITIALIZE
    en1.initialize();
  #endif
  eo1.allocate(1,na,1,ny-1,"eo1");
  #ifndef NO_AD_INITIALIZE
    eo1.initialize();
  #endif
  ettq.allocate(1,nyt,"ettq");
  #ifndef NO_AD_INITIALIZE
    ettq.initialize();
  #endif
  ett.allocate(1,nyt,"ett");
  #ifndef NO_AD_INITIALIZE
    ett.initialize();
  #endif
  ett1.allocate(1,nyt,"ett1");
  #ifndef NO_AD_INITIALIZE
    ett1.initialize();
  #endif
  ett2.allocate(1,nyt,"ett2");
  #ifndef NO_AD_INITIALIZE
    ett2.initialize();
  #endif
  etp.allocate(1,nyp,"etp");
  #ifndef NO_AD_INITIALIZE
    etp.initialize();
  #endif
  etp1.allocate(1,nyp,"etp1");
  #ifndef NO_AD_INITIALIZE
    etp1.initialize();
  #endif
  etp2.allocate(1,nyp,"etp2");
  #ifndef NO_AD_INITIALIZE
    etp2.initialize();
  #endif
  ete.allocate(1,ny,"ete");
  #ifndef NO_AD_INITIALIZE
    ete.initialize();
  #endif
  ecpue.allocate(1,ny,"ecpue");
  #ifndef NO_AD_INITIALIZE
    ecpue.initialize();
  #endif
  cvcpue.allocate(1,ny,"cvcpue");
  #ifndef NO_AD_INITIALIZE
    cvcpue.initialize();
  #endif
  tb.allocate(1,ny,"tb");
  #ifndef NO_AD_INITIALIZE
    tb.initialize();
  #endif
  mp1.allocate(1,na,"mp1");
  #ifndef NO_AD_INITIALIZE
    mp1.initialize();
  #endif
  mp2.allocate(1,na,"mp2");
  #ifndef NO_AD_INITIALIZE
    mp2.initialize();
  #endif
  mp3.allocate(1,na,"mp3");
  #ifndef NO_AD_INITIALIZE
    mp3.initialize();
  #endif
  mp0.allocate(1,na,"mp0");
  #ifndef NO_AD_INITIALIZE
    mp0.initialize();
  #endif
  sel92.allocate(1,na,"sel92");
  #ifndef NO_AD_INITIALIZE
    sel92.initialize();
  #endif
  sel93.allocate(1,na,"sel93");
  #ifndef NO_AD_INITIALIZE
    sel93.initialize();
  #endif
  sel08.allocate(1,na,"sel08");
  #ifndef NO_AD_INITIALIZE
    sel08.initialize();
  #endif
  sel.allocate(1,na,1,ny,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  selt.allocate(1,na,"selt");
  #ifndef NO_AD_INITIALIZE
    selt.initialize();
  #endif
  selp.allocate(1,na,"selp");
  #ifndef NO_AD_INITIALIZE
    selp.initialize();
  #endif
  selw.allocate(1,na,"selw");
  #ifndef NO_AD_INITIALIZE
    selw.initialize();
  #endif
  Mn.allocate(1,na,"Mn");
  #ifndef NO_AD_INITIALIZE
    Mn.initialize();
  #endif
  Mo.allocate(1,na,"Mo");
  #ifndef NO_AD_INITIALIZE
    Mo.initialize();
  #endif
  sos.allocate(1,na,"sos");
  #ifndef NO_AD_INITIALIZE
    sos.initialize();
  #endif
  q1.allocate("q1");
  #ifndef NO_AD_INITIALIZE
  q1.initialize();
  #endif
  q2.allocate("q2");
  #ifndef NO_AD_INITIALIZE
  q2.initialize();
  #endif
  q3.allocate("q3");
  #ifndef NO_AD_INITIALIZE
  q3.initialize();
  #endif
  pre.allocate("pre");
  #ifndef NO_AD_INITIALIZE
  pre.initialize();
  #endif
  bc.allocate(1,ny,"bc");
  #ifndef NO_AD_INITIALIZE
    bc.initialize();
  #endif
  FSum.allocate(1,ny+Nproj,"FSum");
  #ifndef NO_AD_INITIALIZE
    FSum.initialize();
  #endif
  FWin.allocate(1,ny+Nproj,"FWin");
  #ifndef NO_AD_INITIALIZE
    FWin.initialize();
  #endif
  FSub.allocate(1,ny+Nproj,"FSub");
  #ifndef NO_AD_INITIALIZE
    FSub.initialize();
  #endif
  FratSub.allocate("FratSub");
  #ifndef NO_AD_INITIALIZE
  FratSub.initialize();
  #endif
  FratWin.allocate("FratWin");
  #ifndef NO_AD_INITIALIZE
  FratWin.initialize();
  #endif
  Fmult.allocate("Fmult");
  #ifndef NO_AD_INITIALIZE
  Fmult.initialize();
  #endif
  CatSum.allocate("CatSum");
  #ifndef NO_AD_INITIALIZE
  CatSum.initialize();
  #endif
  CatSub.allocate("CatSub");
  #ifndef NO_AD_INITIALIZE
  CatSub.initialize();
  #endif
  CatWin.allocate("CatWin");
  #ifndef NO_AD_INITIALIZE
  CatWin.initialize();
  #endif
  FutCatSum.allocate("FutCatSum");
  #ifndef NO_AD_INITIALIZE
  FutCatSum.initialize();
  #endif
  FutCatSub.allocate("FutCatSub");
  #ifndef NO_AD_INITIALIZE
  FutCatSub.initialize();
  #endif
  FutCatWin.allocate("FutCatWin");
  #ifndef NO_AD_INITIALIZE
  FutCatWin.initialize();
  #endif
  npsEst.allocate(1,na,1,ny+Nproj,"npsEst");
  #ifndef NO_AD_INITIALIZE
    npsEst.initialize();
  #endif
  opsEst.allocate(1,na,1,ny+Nproj,"opsEst");
  #ifndef NO_AD_INITIALIZE
    opsEst.initialize();
  #endif
  npwEst.allocate(1,na,1,ny+Nproj,"npwEst");
  #ifndef NO_AD_INITIALIZE
    npwEst.initialize();
  #endif
  opwEst.allocate(1,na,1,ny+Nproj,"opwEst");
  #ifndef NO_AD_INITIALIZE
    opwEst.initialize();
  #endif
  mmbEst.allocate(1,ny+Nproj,"mmbEst");
  #ifndef NO_AD_INITIALIZE
    mmbEst.initialize();
  #endif
  recEst.allocate(1,ny+Nproj,"recEst");
  #ifndef NO_AD_INITIALIZE
    recEst.initialize();
  #endif
  npsSto.allocate(1,na,1,ny+Nproj,"npsSto");
  #ifndef NO_AD_INITIALIZE
    npsSto.initialize();
  #endif
  opsSto.allocate(1,na,1,ny+Nproj,"opsSto");
  #ifndef NO_AD_INITIALIZE
    opsSto.initialize();
  #endif
  npwSto.allocate(1,na,1,ny+Nproj,"npwSto");
  #ifndef NO_AD_INITIALIZE
    npwSto.initialize();
  #endif
  opwSto.allocate(1,na,1,ny+Nproj,"opwSto");
  #ifndef NO_AD_INITIALIZE
    opwSto.initialize();
  #endif
  mmbSto.allocate(1,ny+Nproj,"mmbSto");
  #ifndef NO_AD_INITIALIZE
    mmbSto.initialize();
  #endif
  recSto.allocate(1,ny+Nproj,"recSto");
  #ifndef NO_AD_INITIALIZE
    recSto.initialize();
  #endif
  EstInd.allocate(1,MaxSim,ny,ny+Nproj,"EstInd");
  #ifndef NO_AD_INITIALIZE
    EstInd.initialize();
  #endif
  EstRec.allocate(1,MaxSim,ny,ny+Nproj,"EstRec");
  #ifndef NO_AD_INITIALIZE
    EstRec.initialize();
  #endif
  MMB0.allocate("MMB0");
  #ifndef NO_AD_INITIALIZE
  MMB0.initialize();
  #endif
  R0.allocate("R0");
  #ifndef NO_AD_INITIALIZE
  R0.initialize();
  #endif
  PassCatch.allocate("PassCatch");
  #ifndef NO_AD_INITIALIZE
  PassCatch.initialize();
  #endif
  RecPass.allocate("RecPass");
  #ifndef NO_AD_INITIALIZE
  RecPass.initialize();
  #endif
  BmsyProx.allocate("BmsyProx");
  #ifndef NO_AD_INITIALIZE
  BmsyProx.initialize();
  #endif
  fbest.allocate("fbest");
  #ifndef NO_AD_INITIALIZE
  fbest.initialize();
  #endif
 fbest = 1.0e+20;
  tf1.allocate("tf1");
  #ifndef NO_AD_INITIALIZE
  tf1.initialize();
  #endif
  tf2.allocate("tf2");
  #ifndef NO_AD_INITIALIZE
  tf2.initialize();
  #endif
  tf3.allocate("tf3");
  #ifndef NO_AD_INITIALIZE
  tf3.initialize();
  #endif
  tf4.allocate("tf4");
  #ifndef NO_AD_INITIALIZE
  tf4.initialize();
  #endif
  tf5.allocate("tf5");
  #ifndef NO_AD_INITIALIZE
  tf5.initialize();
  #endif
  tf6.allocate("tf6");
  #ifndef NO_AD_INITIALIZE
  tf6.initialize();
  #endif
  tf7.allocate("tf7");
  #ifndef NO_AD_INITIALIZE
  tf7.initialize();
  #endif
  tf8.allocate("tf8");
  #ifndef NO_AD_INITIALIZE
  tf8.initialize();
  #endif
  tf9.allocate("tf9");
  #ifndef NO_AD_INITIALIZE
  tf9.initialize();
  #endif
  tf10.allocate("tf10");
  #ifndef NO_AD_INITIALIZE
  tf10.initialize();
  #endif
  T_var.allocate(1,ny,"T_var");
  #ifndef NO_AD_INITIALIZE
    T_var.initialize();
  #endif
  expn.allocate(1,na-1,"expn");
  #ifndef NO_AD_INITIALIZE
    expn.initialize();
  #endif
  expo.allocate(1,na-2,"expo");
  #ifndef NO_AD_INITIALIZE
    expo.initialize();
  #endif
  npp.allocate(1,na,"npp");
  #ifndef NO_AD_INITIALIZE
    npp.initialize();
  #endif
  opp.allocate(1,na-1,"opp");
  #ifndef NO_AD_INITIALIZE
    opp.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  last_y.allocate(1,na,"last_y");
  legaln.allocate(1,ny,"legaln");
  legalb.allocate(1,ny,"legalb");
  mmb0.allocate(1,ny,"mmb0");
  legal.allocate(1,ny+Nproj,"legal");
  #ifndef NO_AD_INITIALIZE
    legal.initialize();
  #endif
  mmb.allocate(1,ny+Nproj,"mmb");
  #ifndef NO_AD_INITIALIZE
    mmb.initialize();
  #endif
  last_legal.allocate("last_legal");
  last_subl.allocate("last_subl");
  last_mmb.allocate("last_mmb");
  last_ofl.allocate("last_ofl");
 cout << "Parameter Section Completed" << endl;
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1E-6}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::initializationfunction(void)
{
  log_q1.set_initial_value(-11.5);
  log_q2.set_initial_value(-11.5);
  log_q3.set_initial_value(-11.05);
  r1.set_initial_value(0.6);
  log_mo1.set_initial_value(-2.41358);
  log_mo2.set_initial_value(4.8211);
  sw3.set_initial_value(1.0);
  log_st1.set_initial_value(-2.5);
  log_st2.set_initial_value(1.00);
  log_sc1.set_initial_value(-2.31267);
  log_sc2.set_initial_value(4.485);
  log_sc3.set_initial_value(-2.41358);
  log_sc4.set_initial_value(4.82115);
  log_sc5.set_initial_value(-1.73);
  log_sc6.set_initial_value(4.68);
  log_sw1.set_initial_value(-2.31267);
  log_sw2.set_initial_value(4.485);
  sos0.set_initial_value(1.0);
  sow.set_initial_value(1.00);
  p4.set_initial_value(0.5);
  advar.set_initial_value(0.5);
  qtad.set_initial_value(1.0);
  qtno.set_initial_value(1.0);
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
  int i,j,Isim,iy;
  dvariable ErrComm;
  cout << "Starting preliminary calcs" << endl;
  qtad = qta;
  qtno = qta;
  log_initpop = log_initp;
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
  
  // Calculate maximum effective smple sieze for Tral survey
  for (i=1;i<=nyt;i++)
   {
    st(i) *= 0.5;
    if (st(i) > maxss) st(i)=maxss*efn;             
   }
  // Calculate maximum effective smple sieze for Pot survey 
  for (i=1;i<=nyp;i++)
   {
    sp(i) *= 0.5;
    if (sp(i) > maxss) sp(i)=maxss*efn;             
   }
  // Calculate maximum effective smple sieze for Winter pot survey
  for (i=1;i<=nyw;i++)
   {
    sw(i) *= 0.10;
    if (sw(i) > 0.5*maxss) sw(i)=0.5*maxss*efn;             
   }
   
  for (i=1;i<=ny;i++)
   {
    sc(i) *= 0.10;
    if (sc(i) > 0.5*maxss) sc(i)=0.5*maxss*efn;             //summer fishery
   }
  for (i=1;i<=nyo;i++)
   {
    so(i) *= 0.10;
    if (so(i) > 0.5*maxss) so(i)=0.5*maxss*efn;             //observer's data
	
   }
  for (i=1;i<=nyt;i++)
   {
    for (j=3;j<=na;j++) ont(j,i) += oot(j,i);
	// tt1 non-legal abundance trawl survey
    tt1(i) = tt(i)*(ont(1,i)+ont(2,i)+ont(3,i));
    // tt2 legal abundance trawl survey
	tt2(i) = tt(i) - tt1(i);
   }
   
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
    // Process error
    for (iy=ny;iy<=ny+Nproj;iy++)
     EstRec(Isim,iy) = mfexp(randn(r)*SigmaR-SigmaR*SigmaR/2.0);
    // Biomass estimation (common to all years)
    ErrComm = mfexp(randn(r)*CVEstExtra-square(CVEstExtra)/2.0);
    for (iy=ny;iy<=ny+Nproj;iy++)
     EstInd(Isim,iy) = ErrComm*mfexp(randn(r)*CVEstWithin-square(CVEstWithin)/2.0);
   }
  cvcpue = elem_div(secpue+1.e-3,stcpue+1.e-3);
  
  mcmc2 << "Year MMB MMN/MMB0 True-OFL Est-OFl ABC/Catch Recr" << endl;
  cout << "End preliminary calcs" << endl; 
  cout << log_initpop << endl;  
  
}

void model_parameters::userfunction(void)
{
  f =0.0;
  dvariable SBPR;
  int i;
  fn_count += 1;
  convert_parameters_into_rates();
  get_first_year_abundance();
  get_number_by_size();
  get_proportion_and_effort();
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
}

void model_parameters::convert_parameters_into_rates(void)
{
  int i, j;
  dvariable mo1,mo2;    
  dvariable sc1,sc2,sc3,sc4, sc5, sc6, st1,st2, sp1,sp2,sw1,sw2;  
  dvariable pp;                                       
   q1 = mfexp(log_q1);
   q2 = mfexp(log_q2);
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
   sp1 = mfexp(log_sp1);
   sp2 = mfexp(log_sp2);
   sw1 = mfexp(log_sw1);
   sw2 = mfexp(log_sw2);
 // calculate molitng and selectivity functions   
   for (j=1;j<=na;j++)
    {
     pp = (double(j)-1.0)*slt;
     pp += slm;
     mp1(j) = 1.0-1.0/(1.0+mfexp(-1.0*mo1*(pp-mo2)));  // molting probability : reverse logistic curve//
     sel92(j) = 1.0/(1.0+mfexp(-1.0*sc1*(pp-sc2)));  // commercial pot selectivity 1976-1992: logistic curve//
     sel93(j) = 1.0/(1.0+mfexp(-1.0*sc3*(pp-sc4)));  // commercial pot selectivity 1993-    : logistic curve//
     selw(j) = 1.0/(1.0+mfexp(-1.0*sw1*(pp-sw2)));   // winter pot selectivity 1976-1992: logistic curve//
     selt(j) = 1.0/(1.0+mfexp(-1.0*st1*(pp-st2)));   // summer trawl net selectivity: logistic curve//
     selp(j) = 1.0/(1.0+mfexp(-1.0*sp1*(pp-sp2))); // summer trawl net selectivity: logistic curve//
    }
   mp1 = mp1/mp1(1);
   sel92 = sel92/sel92(na-1);
   sel92(na) = 1.0;
   sel93 = sel93/sel93(na-1);
   sel93(na) = 1.0;  
   selt = selt/selt(na-1);   
   selt(na) = slt6;  //slt6 is 1.0 by default; however can be changed if we believe in dormshaped selectivity. 
   selp = selp/selp(na-1); 
   selp(na) = 1.0;   
   selw = selw/selw(na-1); 
   selw(na) = sw3;   
   for (i = 1; i<=ny; i++)
    {
     if(i<=scy1)    // period 1976 - 1992 
      for (j=1; j<=na; j++) sel(j,i) = sel92(j); 
	  else        
	  for (j=1; j<=na; j++) sel(j,i) = sel93(j);
    }
}

void model_parameters::get_first_year_abundance(void)
{
  int i,j;
  dvariable first_y, np, op ;
  log_rec = log_relrec+log_recscale;
  first_y = mfexp(log_initpop);
  for (i=1;i<ny;i++) rec(i) = mfexp(log_rec(i));
   // recruits in the last year are assumed to be average of the most recent five years
   rec(ny-1) = 0.2*(rec(ny-2)+rec(ny-3)+rec(ny-4)+rec(ny-5)+rec(ny-6));
   //use the length proportion in the first year to estimate the abundance.
   tb(1) = 0.0;
   ops(1) = 0.0;   
   for (j=1;j<=na-1;j++) expn(j) = mfexp(flnp(j));
   npp(1,na-1) = expn/(1+sum(expn));
   npp(na) = 1-sum(npp(1,na-1));
   for (j=1;j<=na;j++) nps(j,1) = npp(j)*first_y;
   //estimated length proportion for summer commercial fishery.
   for (j=1;j<=na;j++)
    {
     enc(j,1) = 0.0;
     eoc(j,1) = 0.0;
    }
}

void model_parameters::get_number_by_size(void)
{
  int i,j,k;
  dvariable pp,tt1, ttt0;
  dvar_vector tt0(1,na);
  for (i=1;i<ny;i++)
  {
      mp0 = mp1;
    ttt0 = 0.0;
    bc(i) = 0.0;	
    for (k=1;k<=na;k++) {ttt0 += (nps(k,i)+ops(k,i))*lg(k);}	
    if (ttt0 < tc(i)) { ttt0 = tc(i)*1.2;}
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
  last_legal = last_y(3)*lg(3)*wm(3)*sel93(3)+last_y(4)*lg(4)*p4*wm(4)*sel93(4)+last_y(5)*lg(5)*wm(5)*sel93(5)+last_y(6)*lg(6)*wm(6)*sel93(6);
  last_subl = last_y(1)*(1-lg(1))*wm(1)*sel93(1)+last_y(2)*(1-lg(2))*wm(2)*sel93(2)+last_y(3)*(1-lg(3))*wm(3)*sel93(3)+last_y(4)*(1-lg(4))*wm(4)*sel93(4)+last_y(5)*(1-lg(5))*wm(5)*sel93(5)+last_y(6)*(1-lg(6))*wm(6)*sel93(6);
  last_ofl = last_legal*(1.0-exp(-M));
  last_mmb   = mmb(ny);
}

void model_parameters::get_proportion_and_effort(void)
{
  int i,j;
  dvariable bf,af,pp;
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
      et0(j,i) = ((nps(j,it(i))+ops(j,it(i)))*exp(-bf*Mn(j))-(enc(j,it(i))+eoc(j,it(i)))*pct(i)*tc(it(i)))*exp(-af*Mn(j))*selt(j);
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
		if(it(i) < scy1) {ettq(i) = qtno*(ett1(i)+ett2(i));}
			else
		ettq(i) = qtad*(ett1(i)+ett2(i));
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
  for (i=1;i<=ny;i++)
  {
    tb(i) = tb(i) - 0.5*tc(i);  //exploitable abundance at the middle of the season
    if (tb(i) < 0.001) tb(i) = 0.001;
    if (i <= scy1)
	 {
	 ecpue(i) = q1*tb(i);
     ete(i) = tc(i)/(q1*tb(i));     //before 1993
	 }
	else
	 {
	 ete(i) = tc(i)/(q2*tb(i));     //after 1993
	 ecpue(i) = q2*tb(i);	
	 }
  }
  for (i=1;i<=ny;i++)
   {
   if (stcpue(i) <= 0.0) {ecpue(i) = 0.0;}
   }
}

void model_parameters::evaluate_the_objective_function(void)
{
  tf1 = 0.5*norm2(elem_div((log(tt+1.e-3)-log(ettq+1.e-3)),sqrt(log(elem_prod(cv3,cv3)+1.0))));
  f = tf1;
  tf2 = 0.5*norm2(log(tp+1.e-3)-log(qtno*(etp1+etp2)+1.e-3))/log(spcv*spcv+1.0);
  f += wsp(1)*tf2;
   T_var = sqrt(log(elem_prod(cvcpue,cvcpue)+1.0)+advar); 
   tf3 = sum(log(T_var))+0.5*(norm2(elem_div((log(stcpue+1.e-3)-log(ecpue+1.e-3)),T_var))); 
   f += lamc*tf3;   
  tf4 = sum(elem_prod(st,colsum(elem_prod(ont,log(ent+1.e-3))))) - sum(elem_prod(st,colsum(elem_prod(ont,log(ont+1.e-3)))));   
  f -=tf4; 
  tf5 = sum(elem_prod(sp,colsum(elem_prod(onp+oop,log(enp+eop+1.e-3))))) - sum(elem_prod(sp,colsum(elem_prod(onp+oop,log(onp+oop+1.e-3)))));     
  f -= wsp(2)*tf5;
  tf6 = sum(elem_prod(sw,colsum(elem_prod(onw+oow,log(enw+eow+1.e-3))))) - sum(elem_prod(sw,colsum(elem_prod(onw+oow,log(onw+oow+1.e-3)))));      
  f -= tf6;
  tf7 = sum(elem_prod(sc,colsum(elem_prod(onc+ooc,log(enc+eoc+1.e-3))))) - sum(elem_prod(sc,colsum(elem_prod(onc+ooc,log(onc+ooc+1.e-3)))));      
  f -= tf7;
  tf8 = 0.01*norm2(log_relrec);                         //deviation in recruits.
  f += tf8;
  tf9 = sum(elem_prod(so,colsum(elem_prod(ono+ooo,log(eno+eoo+1.e-3))))) - sum(elem_prod(so,colsum(elem_prod(ono+ooo,log(ono+ooo+1.e-3)))));      
  f -= tf9;
}

void model_parameters::Get_Steepness(void)
{
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
   }
}

void model_parameters::Find_OFL(void)
{
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
}

void model_parameters::Proj1Yr(void)
{
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
}

void model_parameters::Find_OFL_EST(void)
{
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
}

void model_parameters::Proj1YrEst(void)
{
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
}

void model_parameters::ProjAll(void)
{
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
}

void model_parameters::ProjConstF(void)
{
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
}

void model_parameters::get_reference_points(void)
{
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
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  cout << "Report Section" << endl;
  get_reference_points();
  FratSub = (FSub(ny-5)+FSub(ny-4)+FSub(ny-3)+FSub(ny-2)+FSub(ny-1))/(FSum(ny-5)+FSum(ny-4)+FSum(ny-3)+FSum(ny-2)+FSum(ny-1));
  FratWin = (FWin(ny-5)+FWin(ny-4)+FWin(ny-3)+FWin(ny-2)+FWin(ny-1))/(FSum(ny-5)+FSum(ny-4)+FSum(ny-3)+FSum(ny-2)+FSum(ny-1));
  report << "nps" << endl << nps << endl;  // Modeled new shell summer 
  report << "ops" << endl << ops << endl;  // Modeled old shell summer
  report << "ett" << endl << ett << endl;  // Estimated trawl abundance
  report << "etp" << endl << etp << endl;  // Estimated Pot abundance
  report << "ete" << endl << ete << endl;  // Estimated Summer fishery effort
  report << "ecpue" << endl << ecpue << endl;  // Estimated Summer fishery cpue
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
  report << "tf1" << endl << tf1 << " " << tf2 << " " << tf3 << " " << tf4 << " " << tf5 << " " << tf6 << " " << tf7 << " " << tf8 << " " << tf9 << " " << tf10 << endl;
}

void model_parameters::final_calcs()
{
 // Output summary stuff
 time(&finish);
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 10000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  time(&start);
  mcmc1.open("Pstar.Out");
  mcmc2.open("FiveY.Out");
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
