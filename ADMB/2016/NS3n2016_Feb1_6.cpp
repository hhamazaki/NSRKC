  #include <math.h>
  #include <admodel.h>
  #include <time.h>
  ofstream mcmc1,mcmc2;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <NS3n2016_Feb1_6.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  ny.allocate("ny");
  na.allocate("na");
  nyt.allocate("nyt");
  nyp.allocate("nyp");
  nyw.allocate("nyw");
  nyf.allocate("nyf");
  nyo.allocate("nyo");
  nys.allocate("nys");
  nyff.allocate("nyff");
  scy.allocate(1,2,"scy");
  scp.allocate("scp");
  slm.allocate("slm");
  slt.allocate("slt");
  M2.allocate("M2");
  ms6.allocate("ms6");
  spcv.allocate("spcv");
  maxss.allocate("maxss");
  maxsc.allocate("maxsc");
  efnt.allocate("efnt");
  efnw.allocate("efnw");
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
  wcpue.allocate(1,nyw,"wcpue");
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
  twc.allocate(1,ny,"twc");
  tws.allocate(1,ny,"tws");
  twst.allocate(1,ny,"twst");
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
  is.allocate(1,nys,"is");
  iss.allocate(1,nys,"iss");
  ins.allocate(1,na,1,nys,"ins");
  ios.allocate(1,na,1,nys,"ios");
  iff.allocate(1,nyff,"iff");
  ifs.allocate(1,nyff,"ifs");
  inff.allocate(1,na,1,nyff,"inff");
  ioff.allocate(1,na,1,nyff,"ioff");
  wm.allocate(1,na,"wm");
  hm.allocate(1,2,"hm");
  tag_recov1.allocate(1,21,1,8,"tag_recov1");
  tag_recov2.allocate(1,21,1,8,"tag_recov2");
  SDRec.allocate("SDRec");
  SDW.allocate("SDW");
  log_initp.allocate("log_initp");
  initp_phase.allocate("initp_phase");
  qtno_phase.allocate("qtno_phase");
  M_phase.allocate("M_phase");
  ms_phase.allocate("ms_phase");
  lamc.allocate("lamc");
  lamw.allocate("lamw");
  lawp.allocate("lawp");
  latag.allocate("latag");
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
 cout << "wm " << wm << endl;
 cout << "latag " << latag << endl;  
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_q1.allocate(-20.5,20.0,1,"log_q1");
  log_q3.allocate(-20.5,20.0,1,"log_q3");
  log_qw.allocate(-10.5,20.0,-1,"log_qw");
  log_initpop.allocate(2.0,15.0,1,"log_initpop");
  log_recscale.allocate(2.0,12.0,1,"log_recscale");
  log_relrec.allocate(1,ny-1,-40.0,40.0,1,"log_relrec");
  flnp.allocate(1,na-1,-10.0,10.0,1,"flnp");
  r1.allocate(0.1,0.9,1,"r1");
  log_mo1.allocate(-10.5,-1.0,1,"log_mo1");
  log_st1.allocate(-15.0,1.0,-1,"log_st1");
  log_st3.allocate(-15.0,1.0,-1,"log_st3");
  log_sw1.allocate(-15.0,1.0,1,"log_sw1");
  sw3.allocate(0.1,1.0001,1,"sw3");
  log_sc1.allocate(-15.0,1.0,1,"log_sc1");
  log_sc3.allocate(-3.0,1.0,-1,"log_sc3");
  log_sc5.allocate(-15.0,1.0,-1,"log_sc5");
  advar.allocate(0.0,6.0,1,"advar");
  qtno.allocate(0.1,1.0,1,"qtno");
  M.allocate(0.1,0.5,M_phase,"M");
  ms.allocate(1.0,5.0,ms_phase,"ms");
  sigma.allocate(0.,30.,2,"sigma");
  ig.allocate(1,2,0.,20.,1,"ig");
  log_rec.allocate(1,ny-1,"log_rec");
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
  eot.allocate(1,na,1,nyt,"eot");
  #ifndef NO_AD_INITIALIZE
    eot.initialize();
  #endif
  ent0.allocate(1,na,1,nyt,"ent0");
  #ifndef NO_AD_INITIALIZE
    ent0.initialize();
  #endif
  eot0.allocate(1,na,1,nyt,"eot0");
  #ifndef NO_AD_INITIALIZE
    eot0.initialize();
  #endif
  enw.allocate(1,na,1,nyw,"enw");
  #ifndef NO_AD_INITIALIZE
    enw.initialize();
  #endif
  eow.allocate(1,na,1,nyw,"eow");
  #ifndef NO_AD_INITIALIZE
    eow.initialize();
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
  enwp.allocate(1,na,1,ny,"enwp");
  #ifndef NO_AD_INITIALIZE
    enwp.initialize();
  #endif
  eowp.allocate(1,na,1,ny,"eowp");
  #ifndef NO_AD_INITIALIZE
    eowp.initialize();
  #endif
  enwpd.allocate(1,na,1,ny,"enwpd");
  #ifndef NO_AD_INITIALIZE
    enwpd.initialize();
  #endif
  eowpd.allocate(1,na,1,ny,"eowpd");
  #ifndef NO_AD_INITIALIZE
    eowpd.initialize();
  #endif
  enwc.allocate(1,na,1,ny,"enwc");
  #ifndef NO_AD_INITIALIZE
    enwc.initialize();
  #endif
  eowc.allocate(1,na,1,ny,"eowc");
  #ifndef NO_AD_INITIALIZE
    eowc.initialize();
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
  ecpue.allocate(1,ny,"ecpue");
  #ifndef NO_AD_INITIALIZE
    ecpue.initialize();
  #endif
  cvcpue.allocate(1,ny,"cvcpue");
  #ifndef NO_AD_INITIALIZE
    cvcpue.initialize();
  #endif
  twsd.allocate(1,ny,"twsd");
  #ifndef NO_AD_INITIALIZE
    twsd.initialize();
  #endif
  tb.allocate(1,ny,"tb");
  #ifndef NO_AD_INITIALIZE
    tb.initialize();
  #endif
  tw.allocate(1,ny,"tw");
  #ifndef NO_AD_INITIALIZE
    tw.initialize();
  #endif
  twp.allocate(1,ny,"twp");
  #ifndef NO_AD_INITIALIZE
    twp.initialize();
  #endif
  ewcpue.allocate(1,nyw,"ewcpue");
  #ifndef NO_AD_INITIALIZE
    ewcpue.initialize();
  #endif
  mp1.allocate(1,na,"mp1");
  #ifndef NO_AD_INITIALIZE
    mp1.initialize();
  #endif
  mp0.allocate(1,na,"mp0");
  #ifndef NO_AD_INITIALIZE
    mp0.initialize();
  #endif
  sel92.allocate(1,na,"sel92");
  #ifndef NO_AD_INITIALIZE
    sel92.initialize();
  #endif
  sel08.allocate(1,na,"sel08");
  #ifndef NO_AD_INITIALIZE
    sel08.initialize();
  #endif
  sel.allocate(1,na,1,ny,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  selt.allocate(1,na,1,ny,"selt");
  #ifndef NO_AD_INITIALIZE
    selt.initialize();
  #endif
  selt1.allocate(1,na,"selt1");
  #ifndef NO_AD_INITIALIZE
    selt1.initialize();
  #endif
  selt2.allocate(1,na,"selt2");
  #ifndef NO_AD_INITIALIZE
    selt2.initialize();
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
  q3.allocate("q3");
  #ifndef NO_AD_INITIALIZE
  q3.initialize();
  #endif
  qw.allocate("qw");
  #ifndef NO_AD_INITIALIZE
  qw.initialize();
  #endif
  pre.allocate("pre");
  #ifndef NO_AD_INITIALIZE
  pre.initialize();
  #endif
  bc.allocate(1,ny,"bc");
  #ifndef NO_AD_INITIALIZE
    bc.initialize();
  #endif
  bcw.allocate(1,ny,"bcw");
  #ifndef NO_AD_INITIALIZE
    bcw.initialize();
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
  tf.allocate(1,13,"tf");
  #ifndef NO_AD_INITIALIZE
    tf.initialize();
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
  tt1.allocate(1,nyt,"tt1");
  #ifndef NO_AD_INITIALIZE
    tt1.initialize();
  #endif
  tt2.allocate(1,nyt,"tt2");
  #ifndef NO_AD_INITIALIZE
    tt2.initialize();
  #endif
  tp1.allocate(1,nyp,"tp1");
  #ifndef NO_AD_INITIALIZE
    tp1.initialize();
  #endif
  tp2.allocate(1,nyp,"tp2");
  #ifndef NO_AD_INITIALIZE
    tp2.initialize();
  #endif
  tagrecap1.allocate(1,na,1,na,"tagrecap1");
  #ifndef NO_AD_INITIALIZE
    tagrecap1.initialize();
  #endif
  tagrecap2.allocate(1,na,1,na,"tagrecap2");
  #ifndef NO_AD_INITIALIZE
    tagrecap2.initialize();
  #endif
  tagrecap3.allocate(1,na,1,na,"tagrecap3");
  #ifndef NO_AD_INITIALIZE
    tagrecap3.initialize();
  #endif
  tagrecap12.allocate(1,na,1,na,"tagrecap12");
  #ifndef NO_AD_INITIALIZE
    tagrecap12.initialize();
  #endif
  tagrecap22.allocate(1,na,1,na,"tagrecap22");
  #ifndef NO_AD_INITIALIZE
    tagrecap22.initialize();
  #endif
  tagrecap32.allocate(1,na,1,na,"tagrecap32");
  #ifndef NO_AD_INITIALIZE
    tagrecap32.initialize();
  #endif
  ptagrecap1.allocate(1,na,1,na,"ptagrecap1");
  #ifndef NO_AD_INITIALIZE
    ptagrecap1.initialize();
  #endif
  ptagrecap2.allocate(1,na,1,na,"ptagrecap2");
  #ifndef NO_AD_INITIALIZE
    ptagrecap2.initialize();
  #endif
  ptagrecap3.allocate(1,na,1,na,"ptagrecap3");
  #ifndef NO_AD_INITIALIZE
    ptagrecap3.initialize();
  #endif
  ptagrecap12.allocate(1,na,1,na,"ptagrecap12");
  #ifndef NO_AD_INITIALIZE
    ptagrecap12.initialize();
  #endif
  ptagrecap22.allocate(1,na,1,na,"ptagrecap22");
  #ifndef NO_AD_INITIALIZE
    ptagrecap22.initialize();
  #endif
  ptagrecap32.allocate(1,na,1,na,"ptagrecap32");
  #ifndef NO_AD_INITIALIZE
    ptagrecap32.initialize();
  #endif
  tgr.allocate(1,na,1,na,"tgr");
  #ifndef NO_AD_INITIALIZE
    tgr.initialize();
  #endif
  mgr1.allocate(1,na,1,na,"mgr1");
  #ifndef NO_AD_INITIALIZE
    mgr1.initialize();
  #endif
  mgr2.allocate(1,na,1,na,"mgr2");
  #ifndef NO_AD_INITIALIZE
    mgr2.initialize();
  #endif
  mgr3.allocate(1,na,1,na,"mgr3");
  #ifndef NO_AD_INITIALIZE
    mgr3.initialize();
  #endif
  egr1.allocate(1,na,1,na,"egr1");
  #ifndef NO_AD_INITIALIZE
    egr1.initialize();
  #endif
  egr2.allocate(1,na,1,na,"egr2");
  #ifndef NO_AD_INITIALIZE
    egr2.initialize();
  #endif
  egr3.allocate(1,na,1,na,"egr3");
  #ifndef NO_AD_INITIALIZE
    egr3.initialize();
  #endif
  egr12.allocate(1,na,1,na,"egr12");
  #ifndef NO_AD_INITIALIZE
    egr12.initialize();
  #endif
  egr22.allocate(1,na,1,na,"egr22");
  #ifndef NO_AD_INITIALIZE
    egr22.initialize();
  #endif
  egr32.allocate(1,na,1,na,"egr32");
  #ifndef NO_AD_INITIALIZE
    egr32.initialize();
  #endif
  ef.allocate(1,6,"ef");
  #ifndef NO_AD_INITIALIZE
    ef.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  last_y.allocate(1,na,"last_y");
  legaln.allocate(1,ny,"legaln");
  legalb.allocate(1,ny,"legalb");
  mmb0.allocate(1,ny+1,"mmb0");
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
  #ifndef NO_AD_INITIALIZE
  last_subl.initialize();
  #endif
  bmsy.allocate("bmsy");
  #ifndef NO_AD_INITIALIZE
  bmsy.initialize();
  #endif
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
  log_q1.set_initial_value(-6.5);
  log_q3.set_initial_value(-6.5);
  log_qw.set_initial_value(-6.5);
  r1.set_initial_value(0.6);
  log_mo1.set_initial_value(-2.41358);
  log_st1.set_initial_value(-2.5);
  log_st3.set_initial_value(-2.5);
  log_sc1.set_initial_value(-2.31267);
  log_sc5.set_initial_value(-2.5);
  log_sw1.set_initial_value(-2.31267);
  advar.set_initial_value(0.5);
  qtno.set_initial_value(1.0);
  sigma.set_initial_value(5.0);
  ig.set_initial_value(10.0);
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  int i,j,Isim,iy;
  dvariable tt0,n0; // Calculated working variable does not need defferentiated
  dvariable ErrComm;
  cout << "Starting preliminary calcs" << endl;
  M = M2;
  ms = ms6;
  Mo = M;
  Mn = M;
  Mo(na) = ms*M;
  Mn(na) = ms*M;
  sos = 1.0;
  selw = 1.0;
  tt0.initialize();
  n0.initialize();
  for (i=1;i<=nyt;i++)
   {
    st(i) *= efnt;
    if (st(i) > maxss) st(i)=maxss;
   }
  for (i=1;i<=nyw;i++)
   {
    sw(i) *= efnw;
    if (sw(i) > maxsc) sw(i)=maxsc;
   }
  for (i=1;i<=ny;i++)
   {
    sc(i) *= efnw;
    if (sc(i) > maxsc) sc(i)=maxsc;
   }
  for (i=1;i<=nyo;i++)
   {
    so(i) *= efnw;
    if (so(i) > maxsc) so(i)=maxsc;
	
   }
  for (i=1;i<=nyt;i++)
   {
    tt1(i) = tt(i)*(ont(1,i)+ont(2,i)+ont(3,i)+oot(1,i)+oot(2,i)+oot(3,i)); // tt1 sub-legal abundance trawl survey
    tt2(i) = tt(i) - tt1(i);                     // tt2 legal abundance trawl survey
   }
  cvcpue = elem_div(secpue+1.e-3,stcpue+1.e-3);
  twsd = twst - tws;
  for (i=1;i<=ny;i++)
  {
  if (twsd(i) > 0)
        {
          tt0 += twsd(i)/tws(i);  // Sum proprotion of discards
          n0 += 1;                // Number of sample 
        }
    }
  for (i=1;i<=ny;i++)
    {
    if (twsd(i)< 0)
        {
          twsd(i) = tws(i)*(tt0/n0);
        }   
    }
 cout << twsd << endl;
  
  
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
  
  mcmc2 << "Year MMB MMN/MMB0 True-OFL Est-OFl ABC/Catch Recr" << endl;
 
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
          
  for(i=1;i<=21;i++)
  {
        tagrecap1(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,3);
        tagrecap2(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,4);
        tagrecap3(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,5);
        tagrecap12(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,3);
        tagrecap22(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,4);       
        tagrecap32(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,5);  
  }
  
  
  for(i=1;i<=na;i++)
  {  
  ptagrecap1(i)=tagrecap1(i)/rowsum(tagrecap1)(i);
  ptagrecap2(i)=tagrecap2(i)/rowsum(tagrecap2)(i);
  ptagrecap3(i)=tagrecap3(i)/rowsum(tagrecap3)(i);
  ptagrecap12(i)=tagrecap12(i)/rowsum(tagrecap12)(i);
  ptagrecap22(i)=tagrecap22(i)/rowsum(tagrecap22)(i);  
  ptagrecap32(i)=tagrecap32(i)/rowsum(tagrecap32)(i);
  }
  ptagrecap1(1) = 0.0;
  ptagrecap22(6) = 0.0;
  ptagrecap32(5) = 0.0;  
  ptagrecap32(6) = 0.0;    
 
 cout << "End preliminary calcs" << endl; 
 
}

void model_parameters::userfunction(void)
{
  f =0.0;
  dvariable SBPR;
  int i;
  fn_count += 1;
  convert_parameters_into_rates();
  get_first_year_abundance();
  growth_matrix();
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
           mcmc2 << " " << FirstYr+i-1 << " " << mmb(i) << " " << mmb(i)/MMB0 << " -1 -1 " << tc(i)+twc(i)+tws(i) << " " << rec(i) << endl;
          ProjAll();
         }
       }
     }
   }
}

void model_parameters::convert_parameters_into_rates(void)
{
  int i, j;
  dvariable mo1;    
  dvariable sc1,sc3,sc5, st1,st3,sw1;  
  double pp, L,Lm;                                       
   q1 = mfexp(log_q1);
   q3 = mfexp(log_q3);
   qw = mfexp(log_qw);   
   mo1 = mfexp(log_mo1);
   sc1 = mfexp(log_sc1);
   sc5 = mfexp(log_sc5);
   st1 = mfexp(log_st1);
   st3 = mfexp(log_st3);
   sw1 = mfexp(log_sw1);
    Lm = slm+slt*(na-1);
 // calculate molting and selectivity functions   
   for (j=1;j<=na;j++)
    {
     L = (double(j)-1.0)*slt;
     L += slm;
     mp1(j) = 1.0-1.0/(1.0+mfexp(mo1*(slm-L)+log(1.0/0.001-1.0)));  // molting probability : reverse logistic function//
     sel92(j) = 1.0/(1.0+mfexp(sc1*(Lm-L)+log(1.0/0.999-1.0)));  // commercial pot selectivity 1976-1992: logistic function//
     sel08(j) = 1.0/(1.0+mfexp(sc1*(Lm-L)+log(1.0/0.999-1.0)));  // commercial pot selectivity 1993 - : logistic function//
     selw(j) = 1.0-1.0/(1.0+mfexp(sw1*(slm-L)+log(1.0/0.001-1.0)));  // winter pot selectivity 1981-2011: reverse logistic function//
     selt1(j) = 1.0/(1.0+mfexp(st1*(Lm-L)+log(1.0/0.999-1.0)));   // NOAA summer trawl net selectivity 1976-1992: logistic function//
     selt2(j) = 1.0/(1.0+mfexp(st1*(Lm-L)+log(1.0/0.999-1.0)));   // ADFG summer trawl net selectivity: 1996- logistic function//     
    }
   selw(1) = sw3;   
   for (i = 1; i<=ny; i++)
    {
     if(i<=scy(1))    // period 1976 - 1992 
      for (j=1; j<=na; j++) selt(j,i) = selt1(j);   //NOAA trawl survey 
      else                          
      for (j=1; j<=na; j++) selt(j,i) = selt1(j);   //ADFG trawl survey
    }
   for (i = 1; i<=ny; i++)
    {
     if(i<=scy(1))    // period 1976 - 1992 
      for (j=1; j<=na; j++) sel(j,i) = sel92(j); 
      else                          // period 2008 - present
      for (j=1; j<=na; j++) sel(j,i) = sel08(j);
    }
}

void model_parameters::get_first_year_abundance(void)
{
  int i,j;
  dvariable first_y, np, op ;
  first_y = mfexp(log_initpop);
   for (j=1;j<=(na-1);j++) expn(j) = mfexp(flnp(j));
   npp(1,na-1) = expn/(1+sum(expn));
   npp(na) = 1-sum(npp(1,na-1));
   for (j=1;j<=na;j++) 
	{
		npw(j,1) = npp(j)*first_y;
		opw(j,1) = 0.0;
	}
  log_rec = log_relrec+log_recscale;  
  for (i=1;i<=(ny-1);i++) rec(i) = mfexp(log_rec(i));  
   rec(ny) = 0.2*(rec(ny-1)+rec(ny-2)+rec(ny-3)+rec(ny-4)+rec(ny-5));
}

void model_parameters::growth_matrix(void)
{
   int i, j ; 
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
   for (i=1;i<=(na-1);i++)
   {
    mu = slm+ig(1)+ig(2)*i;
    for (j=i;j<=na;j++)
     {
     fa = slm+slt*(j-1.5);
     fb = slm+slt*(j-0.5);
     tgr(i,j) = cumd_norm((fb-mu)/sigma)-cumd_norm((fa-mu)/sigma);
     }
   }
   tgr(na,na) = 1;
   for (i=1;i<=na;i++)
   {
    tgr(i) = tgr(i)/rowsum(tgr)(i);
   }
   // growth increment of the last length class is 1.0
   for (i=1;i<=na;i++)
   {
    mgr1(i) = tgr(i)*mp1(i);
   }
   for (i=1;i<=na;i++)
   {
    mgr1(i,i) += (1-mp1(i));
   }
 // Calculate expected matrix for year 2-3
   mgr2 = mgr1*mgr1;    //estimated growth matrix Year 2
   mgr3 = mgr2*mgr1;   //estimated growth matrix Year 3
   for (i=1;i<=na;i++)
    {
    egr1(i) = elem_prod(mgr1(i),sel92);
    egr2(i) = elem_prod(mgr2(i),sel92);
    egr3(i) = elem_prod(mgr3(i),sel92);
    egr12(i) = elem_prod(mgr1(i),sel08);
    egr22(i) = elem_prod(mgr2(i),sel08);
    egr32(i) = elem_prod(mgr3(i),sel08);    
    }
   for (i=1;i<=na;i++)
    {
    egr1(i) = egr1(i)/rowsum(egr1)(i);
    egr2(i) = egr2(i)/rowsum(egr2)(i);
    egr3(i) = egr3(i)/rowsum(egr3)(i);
    egr12(i) = egr12(i)/rowsum(egr12)(i);
    egr22(i) = egr22(i)/rowsum(egr22)(i);
    egr32(i) = egr32(i)/rowsum(egr32)(i);   
    }
}

void model_parameters::get_number_by_size(void)
{
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
  for (i=1;i<=ny;i++)
  {
      mp0 = mp1;
	tt1 = 0.0;  //commercial catch by length class 
	tw(i) = 0.0;
	twp(i) = 0.0;
    for (j=1;j<=na;j++)
     {
      tw(i) += (npw(j,i)+opw(j,i))*selw(j)*lg(j); // For commercial Catch 	  
     }
    for (j=1;j<=na;j++)
     {
      twp(i) += (npw(j,i)+opw(j,i))*selw(j); // For commercial Catch 	  
     }
     twp(i) = twp(i) - 0.5*(twc(i)+tws(i));  
    for (j=1;j<=na;j++)
     {
      enwc(j,i) = npw(j,i)*selw(j)*lg(j)/tw(i);
      eowc(j,i) = opw(j,i)*selw(j)*lg(j)/tw(i);
     }
	for (j=1;j<=na;j++)
     {
      tt1 += (npw(j,i)+opw(j,i))*lg(j); 	  
     }
	pp1 =0.0;
    for (j=3;j<=na;j++)
     {
      pp1 += (npw(j,i)+opw(j,i))*selw(j);  // For subsistence catch   
     }
    for (j=3;j<=na;j++)
     {
	  enwp(j,i) = npw(j,i)*selw(j)/pp1;
      eowp(j,i) = opw(j,i)*selw(j)/pp1;
     }
	pp2 = 0.0;
    for (j=1;j<=2;j++)
     {
      pp2 += (npw(j,i)+opw(j,i))*selw(j);  // For subsistence catch   
     }
    for (j=1;j<=2;j++)
     {
	  enwpd(j,i) = npw(j,i)*selw(j)/pp2;
      eowpd(j,i) = opw(j,i)*selw(j)/pp2;	  
	  }
    // Extract exploitation rate
    FSum(i) = tc(i)/tb(i);
    FSub(i) = tws(i)/pp1;
    FWin(i) = twc(i)/tt1;
  bcw(i) = 0.0;
    for (j=1;j<=na;j++)
	{
     nps(j,i) = (npw(j,i)-enwp(j,i)*tws(i)-enwc(j,i)*twc(i)-(twc(i)/tt1)*(npw(j,i))*selw(j)*(1-lg(j))*hm(2)-twsd(i)*enwpd(j,i)*hm(2))*exp(-0.417*Mn(j));  
     ops(j,i) = (opw(j,i)-eowp(j,i)*tws(i)-eowc(j,i)*twc(i)-(twc(i)/tt1)*(opw(j,i))*selw(j)*(1-lg(j))*hm(2)-twsd(i)*eowpd(j,i)*hm(2))*exp(-0.417*Mo(j));	
	 if (nps(j,i) < 0.0) nps(j,i) = 0.001; // summer abundance should not go bellow zero
	 if (ops(j,i) < 0.0) ops(j,i) = 0.001; // summer abundance should not go bellow zero	
     bcw(i) += (twc(i)/tt1)*(npw(j,i)+opw(j,i))*selw(j)*(1-lg(j))*hm(2); 
	}
	tb(i) = 0.0; 
    for (j=1;j<=na;j++)
	{
        tb(i) += (nps(j,i)+ops(j,i)*sos(j))*sel(j,i)*lg(j); //total summer crab abundance available to 
	}		
	// mean exploitable leagal abundance does not go negative
    if (tb(i) < 0.001) tb(i) = 0.001;
	// Calculate proprotion of newshell and oldshell by comm fiish 	
	for (j=1;j<=na;j++)
     {
        enc(j,1) = 0.0;
        eoc(j,1) = 0.0;
	 	enc(j,i) = nps(j,i)*sel(j,i)*lg(j)/tb(i);
		eoc(j,i) = ops(j,i)*sos(j)*sel(j,i)*lg(j)/tb(i);
	}
    ttt0 = 0.0;   // total number of legal crab on July 1st
    bc(i) = 0.0;  // discards biomass
    for (k=1;k<=na;k++) {ttt0 += (nps(k,i)+ops(k,i))*lg(k);}	
    if (ttt0 < tc(i)) {ttt0 = tc(i)*1.2;}
    else if (ttt0 < 0.001) ttt0 = 0.001;
    for (k=1;k<=na;k++)
     {	  
      tt0(k) = (nps(k,i)+ops(k,i))*exp(-ys(i)*Mn(k))-tc(i)*(enc(k,i)+eoc(k,i))-(tc(i)/ttt0)*(nps(k,i)+ops(k,i))*sel(k,i)*(1.0-lg(k))*hm(1);	  
      if (tt0(k) < 0.0) tt0(k)= 0.001;  // Stock gap measure: abundance of each length class should not go below zero
      bc(i) += (tc(i)/ttt0)*(nps(k,i)+ops(k,i))*sel(k,i)*(1.0-lg(k))*hm(1)*wm(k);  // Bycatch biomass 
	  }
    for (j=1;j<=na;j++)
     {
      pp = 0.0; 
      for (k=1;k<=j;k++) pp += tgr(k,j)*tt0(k)*mp0(k); //Each crab molts right after fishery ended      
      npw(j,i+1) = pp*exp(-(0.583-ys(i))*Mn(j));  //New shell crab on Feb 1st 
      if (j==1) npw(j,i+1) += r1*rec(i);
      if (j==2) npw(j,i+1) += (1.0-r1)*rec(i);
      opw(j,i+1) = tt0(j)*(1.0-mp0(j))*exp(-(0.583-ys(i))*Mn(j));   //Old shell crab are unmolted crab  
     }
   }
 for (i=1; i<=ny; i++)
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
  last_y = column(npw,ny+1);
  for (j=1; j<=na; j++) last_y(j) += opw(j,ny+1);
  last_legal = 0.0;
  for (j=1; j<=na; j++) last_legal += last_y(j)*lg(j)*wm(j)*sel08(j);
  last_subl = 0.0;
  for (j=1; j<=na; j++) last_subl += last_y(j)*(1-lg(j))*wm(j)*sel08(j);
  last_mmb = 0.0;
  for (j=3; j<=na; j++) last_mmb += (npw(j,ny+1)+opw(j,ny+1))*wm(j);
  mmb(ny+1) = last_mmb;
  mmb0(ny+1) = last_mmb;
  bmsy = (sum(mmb)- mmb(1)-mmb(2)-mmb(3)-mmb(4))/(ny-3);
  last_ofl = last_legal*(1.0-exp(-M*(last_mmb/bmsy-0.1)/0.9));
}

void model_parameters::get_proportion_and_effort(void)
{
  int i,j;
  dvariable bf,af,pp;
  ett.initialize();
  ett1.initialize();
  eow.initialize();
  eoo.initialize();
  ecpue.initialize();
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
      ent0(j,i) = ((nps(j,it(i)))*exp(-bf*Mn(j))-(enc(j,it(i)))*pct(i)*tc(it(i)))*exp(-af*Mn(j))*selt(j,i);
      eot0(j,i) = ((ops(j,it(i)))*exp(-bf*Mn(j))-(eoc(j,it(i)))*pct(i)*tc(it(i)))*exp(-af*Mn(j))*selt(j,i); 
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
  if(it(i) < scy(1)) {ettq(i) = qtno*ett(i);}
    else
    ettq(i) = (ett(i));
   }
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
      ewcpue(i) = qw*twp(i); 	 
   }
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
  for (i=1;i<=ny;i++)
  {
    tb(i) = tb(i) - 0.5*tc(i);  //exploitable abundance at the middle of the season
    if (tb(i) < 0.001) tb(i) = 0.001;
    if (i <= scy(1))
	 {
	 ecpue(i) = q1*tb(i);
	 }
    else
	 {	
	 ecpue(i) = q3*tb(i);
	 }
  }
  for (i=1;i<=ny;i++)
   {
   if (stcpue(i) <= 0.0) {ecpue(i) = 0.0;}
   }
}

void model_parameters::evaluate_the_objective_function(void)
{
  tf(1) = 0.5*norm2(elem_div((log(tt+1.e-3)-log(ettq+1.e-3)),sqrt(log(elem_prod(cv3,cv3)+1.0))));
   tf(2) = lamw*norm2(log(wcpue+1.e-3)-log(ewcpue+1.e-3))/(2*SDW*SDW); 
   T_var = sqrt(log(elem_prod(cvcpue,cvcpue)+1.0)+advar); 
   tf(3) = lamc*sum(log(T_var))+0.5*(norm2(elem_div((log(stcpue+1.e-3)-log(ecpue+1.e-3)),T_var)));  
  tf(4) = -(sum(elem_prod(st,colsum(elem_prod(ont,log(ent+1.e-3))))) - sum(elem_prod(st,colsum(elem_prod(ont,log(ont+1.e-3))))));   
  tf(5) = -(sum(elem_prod(st,colsum(elem_prod(oot,log(eot+1.e-3))))) - sum(elem_prod(st,colsum(elem_prod(oot,log(oot+1.e-3))))));   
  tf(6) = -lawp*(sum(elem_prod(sw,colsum(elem_prod(onw,log(enw+1.e-3))))) - sum(elem_prod(sw,colsum(elem_prod(onw,log(onw+1.e-3))))));     
  tf(7) = -lawp*(sum(elem_prod(sw,colsum(elem_prod(oow,log(eow+1.e-3))))) - sum(elem_prod(sw,colsum(elem_prod(oow,log(oow+1.e-3))))));       
  tf(8) = -(sum(elem_prod(sc,colsum(elem_prod(onc,log(enc+1.e-3))))) - sum(elem_prod(sc,colsum(elem_prod(onc,log(onc+1.e-3))))));      
  //Log likelhihood size proportion for summer fishery survey    
  tf(9) = -(sum(elem_prod(sc,colsum(elem_prod(ooc,log(eoc+1.e-3))))) - sum(elem_prod(sc,colsum(elem_prod(ooc,log(ooc+1.e-3))))));      
  tf(10) = -(sum(elem_prod(so,colsum(elem_prod(ono,log(eno+1.e-3))))) - sum(elem_prod(so,colsum(elem_prod(ono,log(ono+1.e-3))))));   
  tf(11) = -(sum(elem_prod(so,colsum(elem_prod(ooo,log(eoo+1.e-3))))) - sum(elem_prod(so,colsum(elem_prod(ooo,log(ooo+1.e-3))))));   
  tf(12) = norm2(log_relrec)/(2*SDRec*SDRec);                            
  ef(1) = -(sum(elem_prod(rowsum(tagrecap1),rowsum(elem_prod(ptagrecap1,log(egr1+1.e-3))))-elem_prod(rowsum(tagrecap1),rowsum(elem_prod(ptagrecap1,log(ptagrecap1+1.e-3)))))); 
  ef(2) = -(sum(elem_prod(rowsum(tagrecap2),rowsum(elem_prod(ptagrecap2,log(egr2+1.e-3))))-elem_prod(rowsum(tagrecap2),rowsum(elem_prod(ptagrecap2,log(ptagrecap2+1.e-3))))));  
  ef(3) = -(sum(elem_prod(rowsum(tagrecap3),rowsum(elem_prod(ptagrecap3,log(egr3+1.e-3))))-elem_prod(rowsum(tagrecap3),rowsum(elem_prod(ptagrecap3,log(ptagrecap3+1.e-3))))));  
  ef(4) = -(sum(elem_prod(rowsum(tagrecap12),rowsum(elem_prod(ptagrecap12,log(egr12+1.e-3))))-elem_prod(rowsum(tagrecap12),rowsum(elem_prod(ptagrecap12,log(ptagrecap12+1.e-3)))))); 
  ef(5) = -(sum(elem_prod(rowsum(tagrecap22),rowsum(elem_prod(ptagrecap22,log(egr22+1.e-3))))-elem_prod(rowsum(tagrecap22),rowsum(elem_prod(ptagrecap22,log(ptagrecap22+1.e-3))))));  
  ef(6) = -(sum(elem_prod(rowsum(tagrecap32),rowsum(elem_prod(ptagrecap32,log(egr32+1.e-3))))-elem_prod(rowsum(tagrecap32),rowsum(elem_prod(ptagrecap32,log(ptagrecap32+1.e-3))))));  
  tf(13) = latag*sum(ef);
  f += sum(tf);   
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
    pp += (npw(Len,YrPass)+opw(Len,YrPass))*selw(Len);
    tt1 += (npw(Len,YrPass)+opw(Len,YrPass))*selw(Len)*lg(Len);
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
    eo0f(Len) = opw(Len,YrPass)*selw(Len)/pp;
    CatNum = opw(Len,YrPass)*selw(Len)*Fmult*FratSub;
    CatSub += CatNum;
    TotalCatch += CatNum*wm(Len);
    en1f(Len) = npw(Len,YrPass)*selw(Len)*lg(Len)/tt1;
    CatNum = npw(Len,YrPass)*selw(Len)*lg(Len)*Fmult*FratWin;
    CatWin += CatNum;
    TotalCatch += CatNum*wm(Len);
    eo1f(Len) = opw(Len,YrPass)*selw(Len)*lg(Len)/tt1;
    CatNum = opw(Len,YrPass)*selw(Len)*lg(Len)*Fmult*FratWin;
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
    pp += (npwEst(Len,YrPass)+opwEst(Len,YrPass))*selw(Len);
    tt1 += (npwEst(Len,YrPass)+opwEst(Len,YrPass))*selw(Len)*lg(Len);
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
    eo0f(Len) = opwEst(Len,YrPass)*selw(Len)/pp;
    CatNum = opwEst(Len,YrPass)*selw(Len)*Fmult*FratSub;
    CatSub += CatNum;
    TotalCatch += CatNum*wm(Len);
    en1f(Len) = npwEst(Len,YrPass)*selw(Len)*lg(Len)/tt1;
    CatNum = npwEst(Len,YrPass)*selw(Len)*lg(Len)*Fmult*FratWin;
    CatWin += CatNum;
    TotalCatch += CatNum*wm(Len);
    eo1f(Len) = opwEst(Len,YrPass)*selw(Len)*lg(Len)/tt1;
    CatNum = opwEst(Len,YrPass)*selw(Len)*lg(Len)*Fmult*FratWin;
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
        pp += (npw(Len,FutYr)+opw(Len,FutYr))*selw(Len);
        tt1 += (npw(Len,FutYr)+opw(Len,FutYr))*selw(Len)*lg(Len);
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
        eo0f(Len) = opw(Len,FutYr)*selw(Len)/pp;
        CatNum = opw(Len,FutYr)*selw(Len)*FutMort*FratSub;
        FutCatSub +=  CatNum;
        PredCatchW += CatNum*wm(Len);
        en1f(Len) = npw(Len,FutYr)*selw(Len)*lg(Len)/tt1;
        CatNum = npw(Len,FutYr)*selw(Len)*lg(Len)*FutMort*FratWin;
        FutCatWin += CatNum;
        PredCatchW += CatNum*wm(Len);
        eo1f(Len) = opw(Len,FutYr)*selw(Len)*lg(Len)/tt1;
        CatNum = opw(Len,FutYr)*selw(Len)*lg(Len)*FutMort*FratWin;
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
      pp += (npw(Len,FutYr)+opw(Len,FutYr))*selw(Len);
      tt1 += (npw(Len,FutYr)+opw(Len,FutYr))*selw(Len)*lg(Len);
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
      eo0f(Len) = opw(Len,FutYr)*selw(Len)/pp;
      CatNum = opw(Len,FutYr)*selw(Len)*Fmult*FratSub;
      CatSub +=  CatNum;
      TotalCatch += CatNum*wm(Len);
      en1f(Len) = npw(Len,FutYr)*selw(Len)*lg(Len)/tt1;
      CatNum = npw(Len,FutYr)*selw(Len)*lg(Len)*Fmult*FratWin;
      CatWin += CatNum;
      TotalCatch += CatNum*wm(Len);
      eo1f(Len) = opw(Len,FutYr)*selw(Len)*lg(Len)/tt1;
      CatNum = opw(Len,FutYr)*selw(Len)*lg(Len)*Fmult*FratWin;
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
    for (j=3;j<=na;j++) tt0 += (npw(j,i)+opw(j,i))*selw(j)*lg(j);
    if (tt0 <= 0.0) tt0 = 0.000001;
    ref_Htw(i) = (twc(i)+tws(i))/tt0;
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
        ref_catch_s(k,i) = (ref_nps(k,i)+ref_ops(k,i)*sos(k))*sel(k,ny)*lg(k)*ref_Hts;
        pp = 0;
        for (kk=1;kk<=k;kk++)
        {
          tt0 = (ref_ps(kk,i)*mfexp(-ref_ys*Mn(kk))-ref_catch_s(kk,i)-ref_Hts*ref_ps(kk,i)*sel(kk,ny)*(1.0-lg(kk))*hm(1))*mfexp(-(0.583-ref_ys)*Mn(kk));
          if (tt0<0.0) tt0 = 0.0;
          pp += gr(kk,k)*tt0*mp0(kk);
        }
        ref_npw(k,i) = pp;
        if (k==1) ref_npw(k,i) += r1*1000.0;
        if (k==2) ref_npw(k,i) += (1.0-r1)*1000.0;
        ref_opw(k,i) = (ref_ps(k,i)*mfexp(-ref_ys*Mn(k))-ref_catch_s(k,i)-ref_Hts*ref_ps(k,i)*sel(k,ny)*(1.0-lg(k))*hm(1))*mfexp(-(0.583-ref_ys)*Mn(k))*(1.0-mp0(k));
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
          ref_catch_w_o(k,i) = ref_opw(k,i)*selw(k)*lg(k)*ref_Hwc;
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

void model_parameters::report(const dvector& gradients)
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
  report << "npw" << endl << npw << endl;  // Modeled new shell winter
  report << "opw" << endl << opw << endl;  // Modeled oldshell shell winter  
  report << "ett" << endl << ett << endl;  // Estimated trawl abundance
  report << "ecpue" << endl << ecpue << endl;  // Estimated Summer fishery cpue
  report << "ewcpue" << endl << ewcpue << endl;  // Estimated Winter survey cpue
  report << "ent" << endl << ent << endl;  // Estimated trawl newshell size proportion 
  report << "eot" << endl << eot << endl;  // Estimated trawl oldshell size proportion   
  report << "enw" << endl << enw << endl;  // Estimated winter survey newshell size proportion
  report << "eow" << endl << eow << endl;  // Estimated winter survey newshell size proportion
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
