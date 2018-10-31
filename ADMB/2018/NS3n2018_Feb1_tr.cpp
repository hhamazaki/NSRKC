#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
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
#include <NS3n2018_Feb1_tr.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  fyear.allocate("fyear");
  lyear.allocate("lyear");
  na.allocate("na");
  ra.allocate("ra");
  sla.allocate("sla");
  nyt.allocate("nyt");
  nyw.allocate("nyw");
  nyo.allocate("nyo");
  ntag1.allocate("ntag1");
  ntag2.allocate("ntag2");
  scy.allocate(1,3,"scy");
  qyear.allocate("qyear");
  slm.allocate("slm");
  slt.allocate("slt");
  M2.allocate("M2");
  ms6.allocate("ms6");
  msn.allocate(1,2,"msn");
  maxs.allocate(1,2,"maxs");
  efn.allocate(1,2,"efn");
  hm.allocate(1,2,"hm");
  lg.allocate(1,na,"lg");
  wm.allocate(1,na,"wm");
  it.allocate(1,nyt,"it");
  tt.allocate(1,nyt,"tt");
  pct.allocate(1,nyt,"pct");
  yt.allocate(1,nyt,"yt");
  cv.allocate(1,nyt,"cv");
  nont.allocate(1,na,1,nyt,"nont");
  noot.allocate(1,na,1,nyt,"noot");
  iw.allocate(1,nyw,"iw");
  wcpue.allocate(1,nyw,"wcpue");
  nonw.allocate(1,na,1,nyw,"nonw");
  noow.allocate(1,na,1,nyw,"noow");
  tc.allocate(fyear,lyear,"tc");
  te.allocate(fyear,lyear,"te");
  stcpue.allocate(fyear,lyear,"stcpue");
  secpue.allocate(fyear,lyear,"secpue");
  ys.allocate(fyear,lyear,"ys");
  nonc.allocate(1,na,fyear,lyear,"nonc");
  nooc.allocate(1,na,fyear,lyear,"nooc");
  twc.allocate(fyear,lyear,"twc");
  tws.allocate(fyear,lyear,"tws");
  twst.allocate(fyear,lyear,"twst");
  io.allocate(1,nyo,"io");
  nono.allocate(1,na,1,nyo,"nono");
  nooo.allocate(1,na,1,nyo,"nooo");
  tag_recov1.allocate(1,ntag1,1,5,"tag_recov1");
  tag_recov2.allocate(1,ntag2,1,5,"tag_recov2");
  SDRec.allocate("SDRec");
  SDW.allocate("SDW");
  log_initp.allocate("log_initp");
  initp_phase.allocate("initp_phase");
  qtno_phase.allocate("qtno_phase");
  M_phase.allocate("M_phase");
  ms_phase.allocate("ms_phase");
  rmol_phase.allocate("rmol_phase");
  lamc.allocate("lamc");
  lamw.allocate("lamw");
  lawp.allocate("lawp");
  latag.allocate("latag");
  nst.allocate("nst");
  nsc.allocate("nsc");
  smol.allocate("smol");
  ssc.allocate("ssc");
  sst.allocate("sst");
  ssw.allocate("ssw");
  sig.allocate("sig");
  ssth.allocate("ssth");
  sthlike3.allocate("sthlike3");
  swm.allocate("swm");
  mol2p.allocate("mol2p");
  pwh.allocate("pwh");
 cout << "Data Section Completed" << endl;
 cout << "lg " << lg << endl;
 cout << "twc " << endl << twc << endl;
 cout << "smol " << endl << smol << endl;
 cout << "ssc " << endl << ssc << endl;  
 cout << "nst " << endl << nst << endl;  
 cout << "ssw " << endl << ssw << endl;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_q.allocate(1,2,-20.5,20.0,1,"log_q");
  log_qw.allocate(-10.5,20.0,-1,"log_qw");
  log_initpop.allocate(2.0,15.0,1,"log_initpop");
  log_recscale.allocate(2.0,12.0,1,"log_recscale");
  log_relrec.allocate(fyear,lyear-1,-40.0,40.0,1,"log_relrec");
  flnp.allocate(1,na-1,0.0,10.0,1,"flnp");
  rlnp.allocate(1,ra-1,0.0,10.0,1,"rlnp");
  log_alpha.allocate(-1.0,5.0,-smol,"log_alpha");
  log_beta.allocate(-5.5,-1.0,-smol,"log_beta");
  log_st1.allocate(1,2,-1.0,5.0,-sst,"log_st1");
  log_st2.allocate(1,2,-5.5,1.0,-sst,"log_st2");
  ist.allocate(1,na,0.0,1.0,sst,"ist");
  ist2.allocate(1,na,0.0,1.0,sst,"ist2");
  log_sw1.allocate(-1.0,5.0,-ssw,"log_sw1");
  log_sw2.allocate(-5.5,-1.0,-ssw,"log_sw2");
  log_sw3.allocate(-15.0,1.0,-ssw,"log_sw3");
  sw3.allocate(1,ra,0.,1.0,1,"sw3");
  log_sc1.allocate(1,3,-15.0,1.0,-ssc,"log_sc1");
  log_sc2.allocate(1,3,-15.0,1.0,-ssc,"log_sc2");
  isc.allocate(1,na,0.,1.0,ssc,"isc");
  isc2.allocate(1,na,0.,1.0,ssc,"isc2");
  advar.allocate(0.0,6.0,1,"advar");
  qtno.allocate(0.1,1.0,1,"qtno");
  M.allocate(0.02,1.0,M_phase,"M");
  ms.allocate(1.0,5.0,ms_phase,"ms");
  sigma.allocate(0.,30.,2,"sigma");
  ig.allocate(1,2,0.,20.,-sig,"ig");
  log_m1.allocate(-5.5,-1.0,1,"log_m1");
  log_m2.allocate(0.5,6.0,1,"log_m2");
  ms1.allocate(1,na,1.0,5.0,1,"ms1");
  st.allocate(1,nyt,"st");
  #ifndef NO_AD_INITIALIZE
    st.initialize();
  #endif
  sw.allocate(1,nyw,"sw");
  #ifndef NO_AD_INITIALIZE
    sw.initialize();
  #endif
  sc.allocate(fyear,lyear,"sc");
  #ifndef NO_AD_INITIALIZE
    sc.initialize();
  #endif
  so.allocate(1,nyo,"so");
  #ifndef NO_AD_INITIALIZE
    so.initialize();
  #endif
  ont.allocate(1,na,1,nyt,"ont");
  #ifndef NO_AD_INITIALIZE
    ont.initialize();
  #endif
  oot.allocate(1,na,1,nyt,"oot");
  #ifndef NO_AD_INITIALIZE
    oot.initialize();
  #endif
  onw.allocate(1,na,1,nyw,"onw");
  #ifndef NO_AD_INITIALIZE
    onw.initialize();
  #endif
  oow.allocate(1,na,1,nyw,"oow");
  #ifndef NO_AD_INITIALIZE
    oow.initialize();
  #endif
  onc.allocate(1,na,fyear,lyear,"onc");
  #ifndef NO_AD_INITIALIZE
    onc.initialize();
  #endif
  ooc.allocate(1,na,fyear,lyear,"ooc");
  #ifndef NO_AD_INITIALIZE
    ooc.initialize();
  #endif
  ono.allocate(1,na,1,nyo,"ono");
  #ifndef NO_AD_INITIALIZE
    ono.initialize();
  #endif
  ooo.allocate(1,na,1,nyo,"ooo");
  #ifndef NO_AD_INITIALIZE
    ooo.initialize();
  #endif
  ent.allocate(1,nyt,1,na,"ent");
  #ifndef NO_AD_INITIALIZE
    ent.initialize();
  #endif
  eot.allocate(1,nyt,1,na,"eot");
  #ifndef NO_AD_INITIALIZE
    eot.initialize();
  #endif
  enw.allocate(1,nyw,1,na,"enw");
  #ifndef NO_AD_INITIALIZE
    enw.initialize();
  #endif
  eow.allocate(1,nyw,1,na,"eow");
  #ifndef NO_AD_INITIALIZE
    eow.initialize();
  #endif
  enc.allocate(fyear,lyear,1,na,"enc");
  #ifndef NO_AD_INITIALIZE
    enc.initialize();
  #endif
  eoc.allocate(fyear,lyear,1,na,"eoc");
  #ifndef NO_AD_INITIALIZE
    eoc.initialize();
  #endif
  eno.allocate(1,nyo,1,na,"eno");
  #ifndef NO_AD_INITIALIZE
    eno.initialize();
  #endif
  eoo.allocate(1,nyo,1,na,"eoo");
  #ifndef NO_AD_INITIALIZE
    eoo.initialize();
  #endif
  enwp.allocate(fyear,lyear,1,na,"enwp");
  #ifndef NO_AD_INITIALIZE
    enwp.initialize();
  #endif
  eowp.allocate(fyear,lyear,1,na,"eowp");
  #ifndef NO_AD_INITIALIZE
    eowp.initialize();
  #endif
  enwpd.allocate(fyear,lyear,1,na,"enwpd");
  #ifndef NO_AD_INITIALIZE
    enwpd.initialize();
  #endif
  eowpd.allocate(fyear,lyear,1,na,"eowpd");
  #ifndef NO_AD_INITIALIZE
    eowpd.initialize();
  #endif
  nps.allocate(fyear,lyear+1,1,na,"nps");
  #ifndef NO_AD_INITIALIZE
    nps.initialize();
  #endif
  ops.allocate(fyear,lyear+1,1,na,"ops");
  #ifndef NO_AD_INITIALIZE
    ops.initialize();
  #endif
  npw.allocate(fyear,lyear+1,1,na,"npw");
  #ifndef NO_AD_INITIALIZE
    npw.initialize();
  #endif
  opw.allocate(fyear,lyear+1,1,na,"opw");
  #ifndef NO_AD_INITIALIZE
    opw.initialize();
  #endif
  ent0.allocate(1,nyt,1,na,"ent0");
  #ifndef NO_AD_INITIALIZE
    ent0.initialize();
  #endif
  eot0.allocate(1,nyt,1,na,"eot0");
  #ifndef NO_AD_INITIALIZE
    eot0.initialize();
  #endif
  ett.allocate(1,nyt,"ett");
  #ifndef NO_AD_INITIALIZE
    ett.initialize();
  #endif
  ettq.allocate(1,nyt,"ettq");
  #ifndef NO_AD_INITIALIZE
    ettq.initialize();
  #endif
  enwc.allocate(fyear,lyear,1,na,"enwc");
  #ifndef NO_AD_INITIALIZE
    enwc.initialize();
  #endif
  eowc.allocate(fyear,lyear,1,na,"eowc");
  #ifndef NO_AD_INITIALIZE
    eowc.initialize();
  #endif
  twsd.allocate(fyear,lyear,"twsd");
  #ifndef NO_AD_INITIALIZE
    twsd.initialize();
  #endif
  tb.allocate(fyear,lyear,"tb");
  #ifndef NO_AD_INITIALIZE
    tb.initialize();
  #endif
  twp.allocate(fyear,lyear,"twp");
  #ifndef NO_AD_INITIALIZE
    twp.initialize();
  #endif
  log_rec.allocate(fyear,lyear-1,"log_rec");
  #ifndef NO_AD_INITIALIZE
    log_rec.initialize();
  #endif
  rec.allocate(fyear,lyear,"rec");
  #ifndef NO_AD_INITIALIZE
    rec.initialize();
  #endif
  ecpue.allocate(fyear,lyear,"ecpue");
  #ifndef NO_AD_INITIALIZE
    ecpue.initialize();
  #endif
  cvcpue.allocate(fyear,lyear,"cvcpue");
  #ifndef NO_AD_INITIALIZE
    cvcpue.initialize();
  #endif
  mlen.allocate(1,na,"mlen");
  #ifndef NO_AD_INITIALIZE
    mlen.initialize();
  #endif
  molp.allocate(fyear,lyear,1,na,"molp");
  #ifndef NO_AD_INITIALIZE
    molp.initialize();
  #endif
  mp1.allocate(1,na,"mp1");
  #ifndef NO_AD_INITIALIZE
    mp1.initialize();
  #endif
  selc1.allocate(1,2,1,na,"selc1");
  #ifndef NO_AD_INITIALIZE
    selc1.initialize();
  #endif
  sel.allocate(fyear,lyear,1,na,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  st1.allocate(1,2,"st1");
  #ifndef NO_AD_INITIALIZE
    st1.initialize();
  #endif
  sc1.allocate(1,3,"sc1");
  #ifndef NO_AD_INITIALIZE
    sc1.initialize();
  #endif
  selt1.allocate(1,2,1,na,"selt1");
  #ifndef NO_AD_INITIALIZE
    selt1.initialize();
  #endif
  selt.allocate(fyear,lyear,1,na,"selt");
  #ifndef NO_AD_INITIALIZE
    selt.initialize();
  #endif
  selw.allocate(1,na,"selw");
  #ifndef NO_AD_INITIALIZE
    selw.initialize();
  #endif
  matc.allocate(1,na,"matc");
  #ifndef NO_AD_INITIALIZE
    matc.initialize();
  #endif
  imatc.allocate(1,na,"imatc");
  #ifndef NO_AD_INITIALIZE
    imatc.initialize();
  #endif
  Mn.allocate(1,na,"Mn");
  #ifndef NO_AD_INITIALIZE
    Mn.initialize();
  #endif
  q.allocate(1,2,"q");
  #ifndef NO_AD_INITIALIZE
    q.initialize();
  #endif
  qw.allocate("qw");
  #ifndef NO_AD_INITIALIZE
  qw.initialize();
  #endif
  bc.allocate(fyear,lyear,"bc");
  #ifndef NO_AD_INITIALIZE
    bc.initialize();
  #endif
  bcw.allocate(fyear,lyear,"bcw");
  #ifndef NO_AD_INITIALIZE
    bcw.initialize();
  #endif
  T_var.allocate(fyear,lyear,"T_var");
  #ifndef NO_AD_INITIALIZE
    T_var.initialize();
  #endif
  expn.allocate(1,na-1,"expn");
  #ifndef NO_AD_INITIALIZE
    expn.initialize();
  #endif
  expr.allocate(1,ra-1,"expr");
  #ifndef NO_AD_INITIALIZE
    expr.initialize();
  #endif
  npp.allocate(1,na,"npp");
  #ifndef NO_AD_INITIALIZE
    npp.initialize();
  #endif
  rpp.allocate(1,na,"rpp");
  #ifndef NO_AD_INITIALIZE
    rpp.initialize();
  #endif
  tf.allocate(1,17,"tf");
  #ifndef NO_AD_INITIALIZE
    tf.initialize();
  #endif
  sth2.allocate(1,6,"sth2");
  #ifndef NO_AD_INITIALIZE
    sth2.initialize();
  #endif
  sth3.allocate(1,6,"sth3");
  #ifndef NO_AD_INITIALIZE
    sth3.initialize();
  #endif
  slg.allocate(1,na,"slg");
  #ifndef NO_AD_INITIALIZE
    slg.initialize();
  #endif
  wmlg.allocate(1,na,"wmlg");
  #ifndef NO_AD_INITIALIZE
    wmlg.initialize();
  #endif
  wmslg.allocate(1,na,"wmslg");
  #ifndef NO_AD_INITIALIZE
    wmslg.initialize();
  #endif
  offset.allocate(1,8,"offset");
  #ifndef NO_AD_INITIALIZE
    offset.initialize();
  #endif
  toffset.allocate(1,6,"toffset");
  #ifndef NO_AD_INITIALIZE
    toffset.initialize();
  #endif
  tag1.allocate(1,na,1,na,"tag1");
  #ifndef NO_AD_INITIALIZE
    tag1.initialize();
  #endif
  tag2.allocate(1,na,1,na,"tag2");
  #ifndef NO_AD_INITIALIZE
    tag2.initialize();
  #endif
  tag3.allocate(1,na,1,na,"tag3");
  #ifndef NO_AD_INITIALIZE
    tag3.initialize();
  #endif
  tag12.allocate(1,na,1,na,"tag12");
  #ifndef NO_AD_INITIALIZE
    tag12.initialize();
  #endif
  tag22.allocate(1,na,1,na,"tag22");
  #ifndef NO_AD_INITIALIZE
    tag22.initialize();
  #endif
  tag32.allocate(1,na,1,na,"tag32");
  #ifndef NO_AD_INITIALIZE
    tag32.initialize();
  #endif
  ptag1.allocate(1,na,1,na,"ptag1");
  #ifndef NO_AD_INITIALIZE
    ptag1.initialize();
  #endif
  ptag2.allocate(1,na,1,na,"ptag2");
  #ifndef NO_AD_INITIALIZE
    ptag2.initialize();
  #endif
  ptag3.allocate(1,na,1,na,"ptag3");
  #ifndef NO_AD_INITIALIZE
    ptag3.initialize();
  #endif
  ptag12.allocate(1,na,1,na,"ptag12");
  #ifndef NO_AD_INITIALIZE
    ptag12.initialize();
  #endif
  ptag22.allocate(1,na,1,na,"ptag22");
  #ifndef NO_AD_INITIALIZE
    ptag22.initialize();
  #endif
  ptag32.allocate(1,na,1,na,"ptag32");
  #ifndef NO_AD_INITIALIZE
    ptag32.initialize();
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
  last_legalb.allocate(1,na,"last_legalb");
  #ifndef NO_AD_INITIALIZE
    last_legalb.initialize();
  #endif
  last_sublb.allocate(1,na,"last_sublb");
  #ifndef NO_AD_INITIALIZE
    last_sublb.initialize();
  #endif
  legaln.allocate(fyear,lyear+1,"legaln");
  legalb.allocate(fyear,lyear+1,"legalb");
  mmb.allocate(fyear,lyear+1,"mmb");
  last_ofl.allocate("last_ofl");
  last_subofl.allocate("last_subofl");
  bmsy.allocate("bmsy");
  #ifndef NO_AD_INITIALIZE
  bmsy.initialize();
  #endif
  last_mmb.allocate("last_mmb");
  #ifndef NO_AD_INITIALIZE
  last_mmb.initialize();
  #endif
 cout << "Parameter Section Completed" << endl;
}

void model_parameters::initializationfunction(void)
{
  log_q.set_initial_value(-6.5);
  log_alpha.set_initial_value(-2.4);
  log_beta.set_initial_value(4.5);
  log_st1.set_initial_value(-2.5);
  log_sc1.set_initial_value(-2.31267);
  log_sw1.set_initial_value(-2.31267);
  log_qw.set_initial_value(-6.5);
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
  int i,j;
  dvariable tt0,n0; // Calculated working variables nondefferentiated
  M = M2;
  ms = ms6;
  tt0.initialize();
  n0.initialize();
  matc.initialize();
  for (i=1;i<=na;i++) mlen(i) = slm + (double(i)-1.0)*slt;
  slg = -lg+1.0;                  // proporiton of sub-legal crab per length
  wmlg  = elem_prod(lg,wm);       // lb.proportion of legal
  wmslg = elem_prod(slg,wm);      // lb.proportion of sublegal
  for (i=1;i<=na;i++)
    {
	if (i>ra) matc(i) = 1.0;          // mature crab is 1.0
	}
  imatc = -matc+1.0;                  // winter subistence discards 
 cout << matc << endl;
 cout << imatc << endl;
 
  st = colsum(nont)+colsum(noot);
  sw = colsum(nonw)+colsum(noow);
  sc = colsum(nonc)+colsum(nooc);
  so = colsum(nono)+colsum(nooo);
  
    
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
  for (i=1;i<=nyt;i++)
   {
    st(i) *= efn(1);
    if (st(i) > maxs(1)) st(i)=maxs(1);
   }
  for (i=1;i<=nyw;i++)
   {
    sw(i) *= efn(2);
    if (sw(i) > maxs(2)) sw(i)=maxs(2);
   }
  for (i=fyear;i<=lyear;i++)
   {
    sc(i) *= efn(2);
    if (sc(i) > maxs(2)) sc(i)=maxs(2);
   }
  for (i=1;i<=nyo;i++)
   {
    so(i) *= efn(2);
    if (so(i) > maxs(2)) so(i)=maxs(2);
   }
  offset(1) = sum(elem_prod(st,colsum(elem_prod(ont,log(ont+1.e-3)))));
  offset(2) = sum(elem_prod(st,colsum(elem_prod(oot,log(oot+1.e-3)))));
  offset(3) = sum(elem_prod(sw,colsum(elem_prod(onw,log(onw+1.e-3)))));
  offset(4) = sum(elem_prod(sw,colsum(elem_prod(oow,log(oow+1.e-3)))));
  offset(5) = sum(elem_prod(sc,colsum(elem_prod(onc,log(onc+1.e-3)))));
  offset(6) = sum(elem_prod(sc,colsum(elem_prod(ooc,log(ooc+1.e-3)))));
  offset(7) = sum(elem_prod(so,colsum(elem_prod(ono,log(ono+1.e-3)))));
  offset(8) = sum(elem_prod(so,colsum(elem_prod(ooo,log(ooo+1.e-3)))));
    
  cvcpue = elem_div(secpue+1.e-3,stcpue+1.e-3);
  twsd = twst - tws;
  for (i=fyear;i<=lyear;i++)
  {
  if (twsd(i) > 0)
        {
          tt0 += twsd(i)/tws(i);  // Sum proprotion of discards
          n0 += 1;                // Number of sample 
        }
    }
  for (i=fyear;i<=lyear;i++)
    {
    if (twsd(i)< 0)
        {
          twsd(i) = tws(i)*(tt0/n0);
        }   
    }
   
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
  toffset(1) = sum(elem_prod(rowsum(tag1),rowsum(elem_prod(ptag1,log(ptag1+1.e-3)))));
  toffset(2) = sum(elem_prod(rowsum(tag2),rowsum(elem_prod(ptag2,log(ptag2+1.e-3)))));
  toffset(3) = sum(elem_prod(rowsum(tag3),rowsum(elem_prod(ptag3,log(ptag3+1.e-3)))));
  toffset(4) = sum(elem_prod(rowsum(tag12),rowsum(elem_prod(ptag12,log(ptag12+1.e-3)))));
  toffset(5) = sum(elem_prod(rowsum(tag22),rowsum(elem_prod(ptag22,log(ptag22+1.e-3)))));
  toffset(6) = sum(elem_prod(rowsum(tag32),rowsum(elem_prod(ptag32,log(ptag32+1.e-3)))));
 cout << "End preliminary calcs" << endl; 
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

void model_parameters::userfunction(void)
{
  f =0.0;
  convert_parameters_into_rates();
  get_first_year_abundance();
  growth_matrix();
  get_number_by_size();
  get_proportion_and_effort();
  evaluate_the_objective_function();
}

void model_parameters::convert_parameters_into_rates(void)
{
  int i, j;
  dvariable mol,sw1,alpha,beta,m1,m2;    
  npp.initialize();				
  rpp.initialize();				
   q = mfexp(log_q);
   qw = mfexp(log_qw);   
   sc1 = mfexp(log_sc1);
   st1 = mfexp(log_st1); 
   sw1 = mfexp(log_sw1);
   alpha = mfexp(log_alpha);
   beta = mfexp(log_beta);
   m1 = mfexp(log_m1);
   m2 = mfexp(log_m2);
    if(smol>=1)
	 {
	  mp1 = 1.0/(1.0+mfexp(alpha*(mlen - beta)));		
	 }
     if(nsc==2) selc1(2) = isc2; // Estimate selectivity two periods individul length class	   		 
	 else  // Estimate selectivity as logistic function 
	  {
     selc1(1) = 1.0/(1.0+mfexp(sc1(1) +sc2(1)*mlen));  // commercial pot selectivity 1976-1992: logistic function//
     selc1(2) = 1.0/(1.0+mfexp(sc1(nsc)+sc2(nsc)*mlen));  // commercial pot selectivity 1993 - : logistic function//
	  }
     selw = 1.0/(1.0+mfexp(-sw1+sw2*mlen));  
	 { 
     selt1(1) = 1.0/(1.0+mfexp(st1(1)-st2(1)*mlen));   // NOAA summer trawl net selectivity 1976-1992: logistic function//
     selt1(2) = 1.0/(1.0+mfexp(st1(ns)-st2(ns)*mlen));   // ADFG summer trawl net selectivity: 1996- logistic function//     
	 }
  for (j=1;j<=msn(2);j++) selw(j) = sw3(j);
    ms1(1) = 1.0;
	ms1(3) = 1.0;
	Mn = M*ms1;
 for (i=fyear;i<=lyear;i++)
  {
	 alpha *= exp(log_alpha_dev(i));
	 beta *= exp(log_beta_dev(i));	 
     molp(i) = 1.0/(1.0+mfexp(-alpha + mlen*beta));		
  } 
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(1))    // period 1976 - 1992 
      selt(i) = selt1(1);   //NOAA trawl survey 
      else                          
      selt(i) = selt1(2);   //ADFG trawl survey
    }
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(2))    // period 1976 - 2007
        sel(i) = selc1(1); 
	  else            // period 2008 - present
        sel(i) = selc1(2);
    }
   for (j=1;j<=(na-1);j++) expn(j) = mfexp(flnp(j));
   npp(1,na-1) = expn/(1+sum(expn));
   npp(na) = 1-sum(npp(1,na-1));
   for (j=1;j<=(ra-1);j++) expr(j) = mfexp(rlnp(j));
   rpp(1,ra-1) = expr/(1+sum(expr));
   rpp(ra) = 1-sum(rpp(1,ra-1));
}

void model_parameters::get_first_year_abundance(void)
{
  int i, j;
  dvariable first_y;
  first_y = mfexp(log_initpop);
   for (j=1;j<=na;j++) 
	{
		npw(fyear,j) = npp(j)*first_y;
		opw(fyear,j) = 0.0;
	}
  log_rec = log_relrec+log_recscale;  
  for (i=fyear;i<=(lyear-1);i++) rec(i) = mfexp(log_rec(i));  
   rec(lyear) = 0.2*(rec(lyear-1)+rec(lyear-2)+rec(lyear-3)+rec(lyear-4)+rec(lyear-5));
}

void model_parameters::growth_matrix(void)
{
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
   mu = slm;
   for (i=1;i<=(na-1);i++)
   {
	if(sig==1) mu += iig(i);		   
    else mu = slm+ig(1)+ig(2)*i;
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
    egr1(i) = elem_prod(mgr1(i),selc1(1));
    egr2(i) = elem_prod(mgr2(i),selc1(1));
    egr3(i) = elem_prod(mgr3(i),selc1(1));
    egr12(i) = elem_prod(mgr1(i),selc1(2));
    egr22(i) = elem_prod(mgr2(i),selc1(2));
    egr32(i) = elem_prod(mgr3(i),selc1(2));    
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
  dvariable pp, pp1, pp2, tt1, ttt0, ttt1, sb, tw;
  dvar_vector nscaf(1,na), tsc(1,na),tscd(1,na), tt0(1,na);    // Abundance after summer fishery
  dvar_vector sselw(1,na),dselw(1,na),lselw(1,na),dlselw(1,na);
  dvar_vector sselc(1,na),dselc(1,na); 
  dvar_vector nwc(1,na),owc(1,na),nwcd(1,na),owcd(1,na);
  dvar_vector nws(1,na),ows(1,na),nwsd(1,na),owsd(1,na);  
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
  last_ofl.initialize();
  last_subofl.initialize();
   sselw = elem_prod(matc,selw);   // subsistence harvets take all mature crab caught
   dselw = elem_prod(imatc,selw);  // subsistence harvest discards all immature crab caught
   lselw = elem_prod(lg,selw);     // legal crab caught by winter gear 
   dlselw = elem_prod(slg,selw);     // sublegal crab caught by winter gear 
  for (i=fyear;i<=lyear;i++)
  {
      tw = sum(elem_prod(npw(i)+opw(i),lselw)); // sum of legal by witner commercial 	  
      enwc(i) = elem_prod(npw(i),lselw)/tw;
      eowc(i) = elem_prod(opw(i),lselw)/tw;
      nwc = twc(i)*enwc(i);
      owc = twc(i)*eowc(i);
      nwcd = (twc(i)/tw)*elem_prod(npw(i),dlselw);
      owcd = (twc(i)/tw)*elem_prod(npw(i),dlselw);
	  bcw(i) = sum(nwcd+owcd)*hm(2); 	  
      twp(i) = sum(elem_prod(npw(i)+opw(i),selw)); //   
      twp(i) = twp(i) - 0.5*(twc(i)+tws(i));  
     pp1 = sum(elem_prod(npw(i)+opw(i),sselw));
	 enwp(i) = elem_prod(npw(i),sselw)/pp1;
     eowp(i) = elem_prod(opw(i),sselw)/pp1;
	 nws = tws(i)*enwp(i);
     ows = tws(i)*eowp(i); 
     pp2 = sum(elem_prod(npw(i)+opw(i),dselw));
	 enwpd(i) = elem_prod(npw(i),dselw)/pp2;
     eowpd(i) = elem_prod(opw(i),dselw)/pp2;
	 nwsd = twsd(i)*enwpd(i);
     owsd = twsd(i)*eowpd(i); 	 
    nps(i) = elem_prod(npw(i)-nwc-nws-nwcd*hm(2)-nwsd*hm(2),exp(-0.417*Mn));  
    ops(i) = elem_prod(opw(i)-owc-ows-owcd*hm(2)-owsd*hm(2),exp(-0.417*Mn));  
    for (j=1;j<=na;j++)
	{
	 if (nps(i,j) < 0.0) nps(i,j) = 0.001; // summer abundance should not go bellow zero
	 if (ops(i,j) < 0.0) ops(i,j) = 0.001; // summer abundance should not go bellow zero	
	}
   sselc = elem_prod(sel(i),lg);   
   dselc = elem_prod(sel(i),slg);     
   tb(i) = sum(elem_prod(nps(i)+ops(i),sselc));	
	// Calculate proprotion of newshell and oldshell by comm fiish 		
	enc(i) = elem_prod(nps(i),sselc)/tb(i);
	eoc(i) = elem_prod(ops(i),sselc)/tb(i);
	enc(fyear) = 0.0;
    eoc(fyear) = 0.0;
	enc(1991) = 0.0;
    eoc(1991) = 0.0;
	tsc = tc(i)*(enc(i)+eoc(i));
    tscd = (tc(i)/tb(i))*elem_prod(nps(i)+ops(i),dselc);
    bc(i) = hm(1)*sum(elem_prod(tscd,wm));  // Bycatch biomass      
 //   else if (ttt0 < 0.001) ttt0 = 0.001;
    nscaf = elem_prod(nps(i)+ops(i),exp(-ys(i)*Mn))-tsc-hm(1)*tscd;	  
    for (k=1;k<=na;k++)
     {	  
      if (nscaf(k) < 0.0) nscaf(k)= 0.001;  // Stop gap measure: abundance of each length class should not go below zero
	  }
    for (j=1;j<=na;j++)
     {
      pp = 0.0; 
      for (k=1;k<=j;k++) pp += tgr(k,j)*nscaf(k)*molp(i,k); //Each crab molts right after fishery ended      
      npw(i+1,j) = pp*exp(-(0.583-ys(i))*Mn(j)) + rpp(j)*rec(i);   //New shell crab molted + recruits 
      opw(i+1,j) = nscaf(j)*(1.0-molp(i,j))*exp(-(0.583-ys(i))*Mn(j));  //Old shell crab are unmolted crab  
     }
   }
 for (i=fyear; i<=lyear+1; i++)
   {
    mmb(i) = sum(elem_prod(elem_prod(npw(i)+opw(i),matc),wm));  // Mature male biomass
    legaln(i) = sum(elem_prod(npw(i)+opw(i),lg));  // number of legal crab
    legalb(i) = sum(elem_prod(npw(i)+opw(i),wmlg)); // legal biomass   
	}
  last_y = npw(lyear+1)+opw(lyear+1);
  last_legalb = elem_prod(elem_prod(last_y,wmlg),selc1(2));
  last_sublb = elem_prod(elem_prod(last_y,wmslg),selc1(2));
  last_mmb = sum(elem_prod(elem_prod(last_y,wm),matc));
  mmb(lyear+1) = last_mmb;
  for(k=1;k<=na;k++)
	{
  last_ofl += last_legalb(k)*(1-exp(-M-0.42*Mn(k))-(1-exp(-0.42*Mn(k)))*((1-pwh*(1-exp(-M-0.42*Mn(k))))/(1-pwh*(1-exp(-0.42*Mn(k))))));
  last_subofl += last_sublb(k)*(1-exp(-M-0.42*Mn(k))-(1-exp(-0.42*Mn(k)))*((1-pwh*(1-exp(-M-0.42*Mn(k))))/(1-pwh*(1-exp(-0.42*Mn(k))))));
	} 
}

void model_parameters::get_proportion_and_effort(void)
{
  int i,j;
  dvariable bf,af,pp;
  ett.initialize();
  eow.initialize();
  eoo.initialize();
  ecpue.initialize();
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
    ent0(i) = elem_prod(elem_prod((elem_prod(nps(it(i)),exp(-bf*Mn))-pct(i)*tc(it(i))*enc(it(i))),exp(-af*Mn)),selt(it(i)));
    eot0(i) = elem_prod(elem_prod((elem_prod(ops(it(i)),exp(-bf*Mn))-pct(i)*tc(it(i))*eoc(it(i))),exp(-af*Mn)),selt(it(i)));
	for (j=1;j<=na;j++)
     {
	  if (ent0(i,j) < 0.0) ent0(i,j) = 0.0;
	  if (eot0(i,j) < 0.0) eot0(i,j) = 0.0;	  
     }
	 ett(i) = sum(ent0(i)+eot0(i));
    if (ett(i) <= 0.0) ett(i) = 0.00001;
    ent(i) = ent0(i)/ett(i);
	eot(i) = eot0(i)/ett(i);
  if(it(i) < qyear) {ettq(i) = qtno*ett(i);}
    else
    ettq(i) = (ett(i));
   }
  for (i=1;i<=nyw;i++)
   {
    pp = sum(elem_prod(npw(iw(i))+opw(iw(i)),selw));
    if (pp <= 0.0) pp = 0.000001;
	enw(i) = elem_prod(npw(iw(i)),selw)/pp;
	eow(i) = elem_prod(opw(iw(i)),selw)/pp;
   }
  for (i=1;i<=nyo;i++)
   {
	pp = sum(elem_prod(elem_prod(nps(io(i))+ops(io(i)),sel(io(i))),slg));
    if (pp <= 0.0) pp = 0.000001;
    eno(i) = elem_prod(elem_prod(nps(io(i)),sel(io(i))),slg)/pp;
    eoo(i) = elem_prod(elem_prod(ops(io(i)),sel(io(i))),slg)/pp;
   }
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
	if (stcpue(i) <= 0.0) ecpue(i) = 0.0;
  }
}

void model_parameters::evaluate_the_objective_function(void)
{
  sth2.initialize();
  sth3.initialize();
  int i,j,k;
  tf(1) = 0.5*norm2(elem_div((log(tt+1.e-3)-log(ettq+1.e-3)),sqrt(log(elem_prod(cv,cv)+1.0))));
  T_var = sqrt(log(elem_prod(cvcpue,cvcpue)+1.0)+advar); 
  tf(3) = lamc*(sum(log(T_var))+0.5*(norm2(elem_div((log(stcpue+1.e-3)-log(ecpue+1.e-3)),T_var))));  
  tf(4) = -(sum(elem_prod(st,rowsum(elem_prod(trans(ont),log(ent+1.e-3))))) - offset(1));   
  tf(5) = -(sum(elem_prod(st,rowsum(elem_prod(trans(oot),log(eot+1.e-3))))) - offset(2));   
  tf(6) = -lawp*(sum(elem_prod(sw,rowsum(elem_prod(trans(onw),log(enw+1.e-3))))) - offset(3));     
  tf(7) = -lawp*(sum(elem_prod(sw,rowsum(elem_prod(trans(oow),log(eow+1.e-3))))) - offset(4));       
  tf(8) = -(sum(elem_prod(sc,rowsum(elem_prod(trans(onc),log(enc+1.e-3))))) - offset(5));      
  //Log likelhihood size proportion for summer fishery survey    
  tf(9) = -(sum(elem_prod(sc,rowsum(elem_prod(trans(ooc),log(eoc+1.e-3))))) - offset(6));      
  tf(10) = -(sum(elem_prod(so,rowsum(elem_prod(trans(ono),log(eno+1.e-3))))) - offset(7));   
  tf(11) = -(sum(elem_prod(so,rowsum(elem_prod(trans(ooo),log(eoo+1.e-3))))) - offset(8));   
  tf(12) = norm2(log_relrec)/(2*SDRec*SDRec);                            
  tf(16) = norm2(log_alpha_dev)/(2*SDRec*SDRec);     
  tf(17) = norm2(log_beta_dev)/(2*SDRec*SDRec);     
  ef(1) = -(sum(elem_prod(rowsum(tag1),rowsum(elem_prod(ptag1,log(egr1+1.e-3)))))- toffset(1)); 
  ef(2) = -(sum(elem_prod(rowsum(tag2),rowsum(elem_prod(ptag2,log(egr2+1.e-3)))))-toffset(2));  
  ef(3) = -(sum(elem_prod(rowsum(tag3),rowsum(elem_prod(ptag3,log(egr3+1.e-3)))))-toffset(3));  
  ef(4) = -(nsc-1)*(sum(elem_prod(rowsum(tag12),rowsum(elem_prod(ptag12,log(egr12+1.e-3)))))-toffset(4)); 
  ef(5) = -(nsc-1)*(sum(elem_prod(rowsum(tag22),rowsum(elem_prod(ptag22,log(egr22+1.e-3)))))-toffset(5));  
  ef(6) = -(nsc-1)*(sum(elem_prod(rowsum(tag32),rowsum(elem_prod(ptag32,log(egr32+1.e-3)))))-toffset(6));  
  tf(14) = latag*sum(ef);
  for (i=1;i<=(na-3);i++)
   {
	if(smol>=1) sth3(1) += square(log(imol(i+3)+1.e-3)-3*log(imol(i+2)+1.e-3)+3*log(imol(i+1)+1.e-3)-log(imol(i)+1.e-3));
	if(sst==1)
	{
	sth3(2) += square(log(ist(i+3))-3*log(ist(i+2))+3*log(ist(i+1))-log(ist(i)));
	if(nst==2) sth3(3) += square(log(ist2(i+3))+3*log(ist2(i+2))+3*log(ist2(i+1))-log(ist2(i)));	
	}
	if(ssc>=1)
	{
	sth3(4) += square(log(isc(i+3))-3*log(isc(i+2))+3*log(isc(i+1))-log(isc(i)));
	if(nsc==2) sth3(5) += square(log(isc2(i+3))-3*log(isc2(i+2))+3*log(isc2(i+1))-log(isc2(i)));
	}
	if(ssw>=1) sth3(6) += square(log(isw(i+3))-3*log(isw(i+2))+3*log(isw(i+1))-log(isw(i)));
	}
  for (i=1;i<=(na-2);i++)
   {
	if(smol>=1) sth2(1) += square(log(imol(i+2))-2*log(imol(i+1))+log(imol(i)));
	if(sst==1)
	{
	sth2(2) += square(log(ist(i+2))-2*log(ist(i+1))+log(ist(i)));
	if(nst==2) sth2(3) += square(log(ist2(i+2))-2*log(ist2(i+1))+log(ist2(i)));	
	}
	if(ssc>=1)
	{
	sth2(4) += square(log(isc(i+2))-2*log(isc(i+1))+log(isc(i)));
	if(nsc==2) sth2(5) += square(log(isc2(i+2))-2*log(isc2(i+1))+log(isc2(i)));
	}
	if(ssw>=1) sth2(6) += square(log(isw(i+2))-2*log(isw(i+1))+log(isw(i)));
	}
	tf(15) = ssth*(sthlike3*sum(sth3)+(1-sthlike3)*sum(sth2));
  f += sum(tf);   
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
  report << "nps" << endl << trans(nps) << endl;  // Modeled new shell summer 
  report << "ops" << endl << trans(ops) << endl;  // Modeled old shell summer
  report << "npw" << endl << trans(npw) << endl;  // Modeled new shell winter
  report << "opw" << endl << trans(opw) << endl;  // Modeled oldshell shell winter  
  report << "ett" << endl << ett << endl;  // Estimated trawl abundance
  report << "ecpue" << endl << ecpue << endl;  // Estimated Summer fishery cpue
  report << "ewcpue" << endl << ewcpue << endl;  // Estimated Winter survey cpue
  report << "ent" << endl << trans(ent) << endl;  // Estimated trawl newshell size proportion 
  report << "eot" << endl << trans(eot) << endl;  // Estimated trawl oldshell size proportion   
  report << "enw" << endl << trans(enw) << endl;  // Estimated winter survey newshell size proportion
  report << "eow" << endl << trans(eow) << endl;  // Estimated winter survey oldshell size proportion
  report << "enc" << endl << trans(enc) << endl;  // Estimated summer fishery newshell size proportion
  report << "eoc" << endl << trans(eoc) << endl;  // Estimated summer fishery oldshell size proportion
  report << "eno" << endl << trans(eno) << endl;  // Estimated observer newshell size proportion
  report << "eoo" << endl << trans(eoo) << endl;  // Estimated observer oldshell size proportion
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
  report << "sth2" << endl << sth2 << endl;
  report << "sth3" << endl << sth3 << endl; 
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
  report << "npp" << endl << npp <<endl; 
  report << "rpp" << endl << rpp <<endl; 
}

void model_parameters::final_calcs()
{
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
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
