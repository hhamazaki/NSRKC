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
#include <NS3n2019_Feb1_tr_new.htp>

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
  nysp.allocate("nysp");
  ntag1.allocate("ntag1");
  ntag2.allocate("ntag2");
  scy.allocate(1,3,"scy");
  qyear.allocate("qyear");
  slm.allocate("slm");
  slt.allocate("slt");
  M2.allocate("M2");
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
  isp.allocate(1,nysp,"isp");
  nonsp.allocate(1,na,1,nysp,"nonsp");
  noosp.allocate(1,na,1,nysp,"noosp");
  tag_recov1.allocate(1,ntag1,1,5,"tag_recov1");
  tag_recov2.allocate(1,ntag2,1,5,"tag_recov2");
  SDRec.allocate("SDRec");
  SDW.allocate("SDW");
  M_phase.allocate("M_phase");
  rmol_phase.allocate("rmol_phase");
  lamc.allocate("lamc");
  lawp.allocate("lawp");
  latag.allocate("latag");
  nst.allocate("nst");
  nsc.allocate("nsc");
  ssc.allocate("ssc");
  bse.allocate(1,2,"bse");
  pwh.allocate("pwh");
 cout << "Data Section Completed" << endl;
 cout << "lg " << lg << endl;
 cout << "twc " << endl << twc << endl;
 cout << "ssc " << endl << ssc << endl;  
 cout << "bse " << endl << bse << endl;  
 cout << "pwh " << endl << pwh << endl;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_q.allocate(1,2,-20.5,10.0,1,"log_q");
  log_qw.allocate(-10.5,20.0,-1,"log_qw");
  log_initpop.allocate(2.0,15.0,1,"log_initpop");
  log_recscale.allocate(2.0,12.0,1,"log_recscale");
  log_relrec.allocate(fyear,lyear-1,-40.0,40.0,1,"log_relrec");
  flnp.allocate(1,na-1,0.0,10.0,1,"flnp");
  rlnp.allocate(1,ra-1,0.0,10.0,1,"rlnp");
  log_amol.allocate(-5.0,-1.0,1,"log_amol");
  log_bmol.allocate(1.0,5.5,1,"log_bmol");
  log_amol_dev.allocate(fyear,lyear,-10.0,10.0,rmol_phase,"log_amol_dev");
  log_bmol_dev.allocate(fyear,lyear,-10.0,10.0,rmol_phase,"log_bmol_dev");
  log_ast.allocate(1,nst,-5,1.0,1,"log_ast");
  log_bst.allocate(1,nst,0.0,6.0,1,"log_bst");
  log_asw.allocate(-5.5,-1.0,1,"log_asw");
  log_bsw.allocate(1.5,6.0,1,"log_bsw");
  sw3.allocate(1,ra,0.,1.0,1,"sw3");
  log_asc.allocate(1,nsc,-5.0,1.0,-ssc,"log_asc");
  log_bsc.allocate(1,nsc,0.0,6.0,-ssc,"log_bsc");
  advar.allocate(0.0,6.0,2,"advar");
  qtno.allocate(0.1,1.0,1,"qtno");
  M.allocate(0.02,1.0,M_phase,"M");
  sigma.allocate(0.,30.,2,"sigma");
  ig.allocate(1,2,0.,20.,1,"ig");
  ms1.allocate(1,na,0.5,5.0,1,"ms1");
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
  ssp.allocate(1,nysp,"ssp");
  #ifndef NO_AD_INITIALIZE
    ssp.initialize();
  #endif
  efst.allocate(1,nyt,"efst");
  #ifndef NO_AD_INITIALIZE
    efst.initialize();
  #endif
  efsw.allocate(1,nyw,"efsw");
  #ifndef NO_AD_INITIALIZE
    efsw.initialize();
  #endif
  efsc.allocate(fyear,lyear,"efsc");
  #ifndef NO_AD_INITIALIZE
    efsc.initialize();
  #endif
  efso.allocate(1,nyo,"efso");
  #ifndef NO_AD_INITIALIZE
    efso.initialize();
  #endif
  efsp.allocate(1,nysp,"efsp");
  #ifndef NO_AD_INITIALIZE
    efsp.initialize();
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
  onsp.allocate(1,na,1,nysp,"onsp");
  #ifndef NO_AD_INITIALIZE
    onsp.initialize();
  #endif
  oosp.allocate(1,na,1,nysp,"oosp");
  #ifndef NO_AD_INITIALIZE
    oosp.initialize();
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
  ensp.allocate(1,nysp,1,na,"ensp");
  #ifndef NO_AD_INITIALIZE
    ensp.initialize();
  #endif
  eosp.allocate(1,nysp,1,na,"eosp");
  #ifndef NO_AD_INITIALIZE
    eosp.initialize();
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
  bc.allocate(fyear,lyear,"bc");
  #ifndef NO_AD_INITIALIZE
    bc.initialize();
  #endif
  bcw.allocate(fyear,lyear,"bcw");
  #ifndef NO_AD_INITIALIZE
    bcw.initialize();
  #endif
  FSub.allocate(fyear,lyear,"FSub");
  #ifndef NO_AD_INITIALIZE
    FSub.initialize();
  #endif
  FSubd.allocate(fyear,lyear,"FSubd");
  #ifndef NO_AD_INITIALIZE
    FSubd.initialize();
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
  q.allocate(1,2,"q");
  #ifndef NO_AD_INITIALIZE
    q.initialize();
  #endif
  T_var.allocate(fyear,lyear,"T_var");
  #ifndef NO_AD_INITIALIZE
    T_var.initialize();
  #endif
  tf.allocate(1,17,"tf");
  #ifndef NO_AD_INITIALIZE
    tf.initialize();
  #endif
  offset.allocate(1,10,"offset");
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
  temp.allocate("temp");
  last_subofl.allocate("last_subofl");
  bmsy.allocate("bmsy");
  last_mmb.allocate("last_mmb");
 cout << "Parameter Section Completed" << endl;
}

void model_parameters::initializationfunction(void)
{
  log_q.set_initial_value(-6.5);
  log_amol.set_initial_value(-2.4);
  log_bmol.set_initial_value(4.5);
  log_ast.set_initial_value(-2.5);
  log_bst.set_initial_value(1.0);
  log_asc.set_initial_value(-2.5);
  log_bsc.set_initial_value(5.5);
  log_asw.set_initial_value(-2.5);
  log_bsw.set_initial_value(4.5);
  log_qw.set_initial_value(-6.5);
  log_amol_dev.set_initial_value(0);
  log_bmol_dev.set_initial_value(0);
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
  tt0.initialize();
  n0.initialize();
  matc.initialize();
  for (i=1;i<=na;i++) mlen(i) = slm + (double(i)-1.0)*slt;
  slg = -lg+1.0;
  wmlg  = elem_prod(lg,wm);       
  wmslg = elem_prod(slg,wm);      
  
  for (i=1;i<=na;i++)
    {
	if (i>ra) matc(i) = 1.0;          // mature crab is 1.0
	}
  imatc = -matc+1.0;                  // inmature crab  
  st = colsum(nont)+colsum(noot);
  sw = colsum(nonw)+colsum(noow);
  sc = colsum(nonc)+colsum(nooc);
  so = colsum(nono)+colsum(nooo);
  ssp = colsum(nonsp)+colsum(noosp);
  
    
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
	
  
  efst = st*efn(1);  
  for (i=1;i<=nyt;i++)
   {
    if (efst(i) > maxs(1)) efst(i) = maxs(1);
   }
   
   efsw = sw*efn(2);
  for (i=1;i<=nyw;i++)
   {
     if (efsw(i) > maxs(2)) efsw(i) = maxs(2);
   }
   efsp = ssp*efn(2);
  for (i=1;i<=nysp;i++)
   {
     if (efsp(i) > maxs(2)) efsp(i) = maxs(2);
   }
   
  efsc = sc*efn(2);
  for (i=fyear;i<=lyear;i++)
   {
    if (efsc(i) > maxs(2)) efsc(i) = maxs(2);
   }
   
   efso = so*efn(2);  
  for (i=1;i<=nyo;i++)
   {
    if (efso(i) > maxs(2)) efso(i)=maxs(2);
   }
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
  dvariable mol,asw,bsw,alpha,beta;    
  dvar_vector ast(1,nst),bst(1,nst), asc(1,nsc),bsc(1,nsc);
  npp.initialize();				
  rpp.initialize();				
   q = mfexp(log_q);
   asc = mfexp(log_asc);
   bsc = mfexp(log_bsc);
   ast = mfexp(log_ast);
   bst = mfexp(log_bst); 
   asw = mfexp(log_asw);
   bsw = mfexp(log_bsw);
   alpha = mfexp(log_amol);
   beta = mfexp(log_bmol);
  mp1 = 1.0/(1.0+mfexp(alpha*(mlen - beta)));
    for (i=fyear;i<=lyear;i++)	
    {	  
	   alpha *= exp(log_amol_dev(i));
	   beta *= exp(log_bmol_dev(i));	 	
       molp(i) = 1.0/(1.0+mfexp(alpha*(mlen - beta)));			
    }
   if(bse(1) <= 1){   
   selc1(1) = 1.0/(1.0+mfexp(asc(1)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));  
   selc1(2) = 1.0/(1.0+mfexp(asc(nsc)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));   
   }  
   if(bse(1) >= 2){    
   selc1(1) = 1.0/(1.0+mfexp(-asc(1)*(mlen - bsc(1)))); 
   selc1(2) = 1.0/(1.0+mfexp(-asc(nsc)*(mlen - bsc(nsc))));  
   }  
   selc1(1,na) = 1.0;
   selc1(2,na) = 1.0;
   selw = 1.0/(1.0+mfexp(asw*(mlen-bsw)));  
   selt1(1) = 1.0/(1.0+mfexp(ast(1)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));   
   selt1(2) = 1.0/(1.0+mfexp(ast(nst)*(-mlen+mlen(na))+log(1.0/0.999-1.0)));        
   selt1(1,na) = 1.0;
   selt1(2,na) = 1.0;
  for (j=1;j<=msn(2);j++) selw(j) = sw3(j);
  for (j=1;j<=(na-msn(1));j++) ms1(j) = 1.0;
    if(bse(2)>=1) ms1(7)=ms1(8);    
	Mn = M*ms1;
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
   sselw = elem_prod(matc,selw); 
   dselw = elem_prod(imatc,selw); 
   lselw = elem_prod(lg,selw);   
   dlselw = elem_prod(slg,selw);      
  for (i=fyear;i<=lyear;i++)
  {
      tw = sum(elem_prod(npw(i)+opw(i),lselw));	  
      nwc = (twc(i)/tw)*elem_prod(npw(i),lselw);
      owc = (twc(i)/tw)*elem_prod(opw(i),lselw);
      nwcd = (twc(i)/tw)*elem_prod(npw(i),dlselw);
      owcd = (twc(i)/tw)*elem_prod(opw(i),dlselw);
	  bcw(i) = sum(nwcd+owcd)*hm(2); 	  
     wpm = sum(elem_prod(npw(i)+opw(i),sselw));
     wpim = sum(elem_prod(npw(i)+opw(i),dselw));
	 nws = (tws(i)/wpm)*elem_prod(npw(i),sselw);
     ows = (tws(i)/wpm)*elem_prod(opw(i),sselw); 
	 nwsd = (twsd(i)/wpim)*elem_prod(npw(i),dselw);
     owsd = (twsd(i)/wpim)*elem_prod(opw(i),dselw); 	 
    FSub(i) = tws(i)/wpm;
    FSubd(i) = twsd(i)/wpim;	
    nps(i) = elem_prod(npw(i)-nwc-nws-nwcd*hm(2)-nwsd*hm(2),exp(-0.417*Mn));  
    ops(i) = elem_prod(opw(i)-owc-ows-owcd*hm(2)-owsd*hm(2),exp(-0.417*Mn));  
    for (j=1;j<=na;j++)
	{
	 if (nps(i,j) < 0.0) nps(i,j) = 0.001; 
	 if (ops(i,j) < 0.0) ops(i,j) = 0.001; 	
	}
   sselc = elem_prod(sel(i),lg);   
   dselc = elem_prod(sel(i),slg);     
   tb(i) = sum(elem_prod(nps(i)+ops(i),sselc));	
	enc(i) = elem_prod(nps(i),sselc)/tb(i);
	eoc(i) = elem_prod(ops(i),sselc)/tb(i);
	enc(fyear) = 0.0;
    eoc(fyear) = 0.0;
	enc(1991) = 0.0;
    eoc(1991) = 0.0;
	tsc = tc(i)*(enc(i)+eoc(i));
    tscd = (tc(i)/tb(i))*elem_prod(nps(i)+ops(i),dselc);
    bc(i) = hm(1)*sum(elem_prod(tscd,wm));     
    nscaf = elem_prod(nps(i)+ops(i),exp(-ys(i)*Mn))-tsc-hm(1)*tscd;	  
    for (k=1;k<=na;k++)
     {	  
      if (nscaf(k) < 0.0) nscaf(k)= 0.001;  
	  }
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
 for (i=fyear; i<=lyear+1; i++)
   {
    mmb(i) = sum(elem_prod(elem_prod(npw(i)+opw(i),matc),wm));  // Mature male biomass
    legaln(i) = sum(elem_prod(npw(i)+opw(i),lg));  // number of legal crab
    legalb(i) = sum(elem_prod(npw(i)+opw(i),wmlg)); // legal biomass   
	}
  last_y = npw(lyear+1)+opw(lyear+1);
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
  FOFL = fl*0.18;
  par1 = 1-exp(-FOFL-0.42*Mn);
  par2 = 1-exp(-0.42*Mn);
  temp = sum(elem_prod(last_legalb,par1-elem_prod(par2,elem_div(1-pwh*par1,1-pwh*par2)))); 
  for(k=1;k<=na;k++)
	{		
  last_ofl += last_legalb(k)*(1-exp(-FOFL(k)-0.42*Mn(k))-(1-exp(-0.42*Mn(k)))*((1-pwh*(1-exp(-FOFL(k)-0.42*Mn(k))))/(1-pwh*(1-exp(-0.42*Mn(k))))));
  last_subofl += last_sublb(k)*(1-exp(-M-0.42*Mn(k))-(1-exp(-0.42*Mn(k)))*((1-pwh*(1-exp(-M-0.42*Mn(k))))/(1-pwh*(1-exp(-0.42*Mn(k))))));
	} 
}

void model_parameters::get_proportion_and_effort(void)
{
  int i,j;
  dvariable bf,af,pp,ts;
  dvar_vector ent0(1,na), eot0(1,na);
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
    ent0 = elem_prod(elem_prod((elem_prod(nps(it(i)),exp(-bf*Mn))-pct(i)*tc(it(i))*enc(it(i))),exp(-af*Mn)),selt(it(i)));
    eot0 = elem_prod(elem_prod((elem_prod(ops(it(i)),exp(-bf*Mn))-pct(i)*tc(it(i))*eoc(it(i))),exp(-af*Mn)),selt(it(i)));
	 ett(i) = sum(ent0+eot0);
    if (ett(i) <= 0.0) ett(i) = 0.00001;
    ent(i) = ent0/ett(i);
	eot(i) = eot0/ett(i); 
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
  if(lawp>=1){
  for (i=1;i<=nysp;i++)
   {   
    pp = sum(elem_prod(nps(isp(i))+nps(isp(i)),selw));
    if (pp <= 0.0) pp = 0.000001;
	ensp(i) = elem_prod(nps(isp(i)),selw)/pp;
	eosp(i) = elem_prod(ops(isp(i)),selw)/pp;
   }
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
    ts = tb(i) - 0.5*tc(i);  //exploitable abundance at the middle of the season
    if (ts < 0.001) ts = 0.001;
    if (i <= scy(1))
	 {
	 ecpue(i) = q(1)*ts;
	 }
    else
	 {	
	 ecpue(i) = q(2)*ts;
	 }
	if (stcpue(i) <= 0.0) ecpue(i) = 0.0;
  }
}

void model_parameters::evaluate_the_objective_function(void)
{
  int i,j,k;
  tf(1) = 0.5*norm2(elem_div((log(tt+1.e-3)-log(ettq+1.e-3)),sqrt(log(elem_prod(cv,cv)+1.0))));
  T_var = sqrt(log(elem_prod(cvcpue,cvcpue)+1.0)+advar); 
  tf(3) = lamc*(sum(log(T_var))+0.5*(norm2(elem_div((log(stcpue+1.e-3)-log(ecpue+1.e-3)),T_var))));  
  tf(4) = -(sum(elem_prod(efst,rowsum(elem_prod(trans(ont),log(ent+1.e-3))))) - offset(1))      
          -(sum(elem_prod(efst,rowsum(elem_prod(trans(oot),log(eot+1.e-3))))) - offset(2));   
  tf(5) = -(sum(elem_prod(efsw,rowsum(elem_prod(trans(onw),log(enw+1.e-3))))) - offset(3))     
          -(sum(elem_prod(efsw,rowsum(elem_prod(trans(oow),log(eow+1.e-3))))) - offset(4));       
  tf(6) = -(sum(elem_prod(efsc,rowsum(elem_prod(trans(onc),log(enc+1.e-3))))) - offset(5))         
          -(sum(elem_prod(efsc,rowsum(elem_prod(trans(ooc),log(eoc+1.e-3))))) - offset(6));      
  tf(7) = -(sum(elem_prod(efso,rowsum(elem_prod(trans(ono),log(eno+1.e-3))))) - offset(7)) 
          -(sum(elem_prod(efso,rowsum(elem_prod(trans(ooo),log(eoo+1.e-3))))) - offset(8));   
  if(lawp>=1){
  tf(8) = -(sum(elem_prod(efsp,rowsum(elem_prod(trans(onsp),log(ensp+1.e-3))))) - offset(9));     
          -(sum(elem_prod(efsp,rowsum(elem_prod(trans(oosp),log(eosp+1.e-3))))) - offset(10));       
	}	  
  tf(9) = norm2(log_relrec)/(2*SDRec*SDRec);                            
  tf(10) = norm2(log_amol_dev)/(2*SDRec*SDRec);     
  tf(11) = norm2(log_bmol_dev)/(2*SDRec*SDRec);     
  ef(1) = -(sum(elem_prod(rowsum(tag1),rowsum(elem_prod(ptag1,log(egr1+1.e-3)))))- toffset(1)); 
  ef(2) = -(sum(elem_prod(rowsum(tag2),rowsum(elem_prod(ptag2,log(egr2+1.e-3)))))-toffset(2));  
  ef(3) = -(sum(elem_prod(rowsum(tag3),rowsum(elem_prod(ptag3,log(egr3+1.e-3)))))-toffset(3));  
  ef(4) = -(nsc-1)*(sum(elem_prod(rowsum(tag12),rowsum(elem_prod(ptag12,log(egr12+1.e-3)))))-toffset(4)); 
  ef(5) = -(nsc-1)*(sum(elem_prod(rowsum(tag22),rowsum(elem_prod(ptag22,log(egr22+1.e-3)))))-toffset(5));  
  ef(6) = -(nsc-1)*(sum(elem_prod(rowsum(tag32),rowsum(elem_prod(ptag32,log(egr32+1.e-3)))))-toffset(6));  
  tf(12) = latag*sum(ef);
  f = sum(tf);   
}

void model_parameters::get_reference_points(void)
{
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
  dvar_vector ref_Htw(1,201);
  dvar_vector ref_mbio(1,201);
  dvar_vector ref_mmb(1,100);
  dvar_vector ref_fcatch_s(1,201);
  dvar_vector ref_fcatch_sd(1,201);
  dvar_vector ref_fcatch_w_c(1,201);
  dvar_vector ref_fcatch_w_s(1,201);
  dvar_vector ref_fcatch_w_cd(1,201);
  dvar_vector ref_fcatch_w_sd(1,201);   
  dvar_vector ref_legaln(1,201);
  dvar_vector ref_Hrt(1,100);
  dvar_vector ref_Hr(1,201);
  dvar_vector ref_legal(1,100);
  dvar_vector ref_legal_b(1,100);
  dvariable ref_ys;
  dvariable f35,f40,h35,h40,i35,i40;
  ref_Htw.initialize();
  ref_FSub.initialize();
  ref_FSubd.initialize();  
  ref_ys = ys(lyear);
  ref_Hwc = 0.0;
  for (i = lyear-10; i<= lyear-1; i++)
   {
    ref_FSub += FSub(i);
	ref_FSubd += FSubd(i);
   }
  ref_FSub = ref_FSub/10.0; 
  ref_FSubd = ref_FSubd/10.0;   
   ref_ys = ys(lyear);
   sselw = elem_prod(matc,selw); 
   dselw = elem_prod(imatc,selw); 
   lselw = elem_prod(lg,selw);   
   dlselw = elem_prod(slg,selw);      
   sselc = elem_prod(sel(lyear),lg);   
   dselc = elem_prod(sel(lyear),slg);  
   ref_rec = rec(lyear);
  ref_Hts = 0.0;
  for (kk = 1; kk <= 201; kk++)
   {
  ref_F = 0.01*kk-0.01;
  par1 = 1-exp(-ref_F-0.42*Mn);
  par2 = 1-pwh*(1-exp(-0.42*Mn));
  ref_Hw = elem_div(pwh*(1-exp(-ref_F))*exp(-0.42*Mn),par2);
  ref_Hs = elem_div((1-pwh)*(1-exp(-ref_F))*exp(-0.42*Mn),par2); 
  ref_nps.initialize();
  ref_ops.initialize();
  ref_npw.initialize();
  ref_opw.initialize();
  ref_npw(1) = npw(lyear);
  ref_opw(1) = opw(lyear);
    for (i=1;i<100;i++)
     { 
	nwc = elem_prod(elem_prod(ref_npw(i),lselw),ref_Hw);
	owc = elem_prod(elem_prod(ref_opw(i),lselw),ref_Hw);
	nwcd = elem_prod(elem_prod(ref_npw(i),dlselw),ref_Hw);
	owcd = elem_prod(elem_prod(ref_opw(i),dlselw),ref_Hw);
    ref_catch_w_c(i) = nwc + owc;
    ref_catch_w_cd(i) = nwcd + owcd;	
		nws = elem_prod(ref_npw(i),sselw)*ref_FSub;
		ows = elem_prod(ref_opw(i),sselw)*ref_FSub;
		nwsd = elem_prod(ref_npw(i),dselw)*ref_FSubd;
        owsd = elem_prod(ref_opw(i),dselw)*ref_FSubd;
    ref_catch_w_s(i) = nws + ows;
    ref_catch_w_sd(i) = nwsd + owsd;
  ref_nps(i) = elem_prod(ref_npw(i)-nwc-nws-nwcd*hm(2)-nwsd*hm(2),exp(-0.417*Mn));  
  ref_ops(i) = elem_prod(ref_opw(i)-owc-ows-owcd*hm(2)-owsd*hm(2),exp(-0.417*Mn));  
	nsc = elem_prod(elem_prod(ref_nps(i),sselc),ref_Hs);
	osc = elem_prod(elem_prod(ref_ops(i),sselc),ref_Hs); 
	nscd = elem_prod(elem_prod(ref_nps(i),dselc),ref_Hs);
	oscd = elem_prod(elem_prod(ref_ops(i),dselc),ref_Hs);
	ref_catch_s(i) = nsc+osc;  	
	ref_catch_sd(i) = nscd+oscd;
  nscaf = elem_prod(ref_nps(i)+ref_ops(i),exp(-ref_ys*Mn))-ref_catch_s(i)-ref_catch_sd(i)*hm(1);	  
    for (k=1;k<=na;k++)
     {	  
      if (nscaf(k) < 0.0) nscaf(k)= 0.001;  
	  }
	ref_Hrt(i) = sum(ref_catch_s(i)+ref_catch_w_c(i))/sum(elem_prod(ref_npw(i)+ref_opw(i),lg));  
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
  ref_mmb(i) = sum(elem_prod(elem_prod(ref_npw(i)+ref_opw(i),matc),wm));  // Mature male biomass
  ref_legal(i) = sum(elem_prod(elem_prod(ref_npw(i)+ref_opw(i),sselc),wm));  // number of legal crab  
   }	  
    ref_mbio(kk) = ref_mmb(99);                    //lbs/R
    ref_legaln(kk) = ref_legal(99);
	ref_fcatch_s(kk) = sum(ref_catch_s(99));
    ref_fcatch_sd(kk) = sum(ref_catch_sd(99));
    ref_fcatch_w_c(kk) = sum(ref_catch_w_c(99));
    ref_fcatch_w_s(kk) = sum(ref_catch_w_s(99));
    ref_fcatch_w_cd(kk) = sum(ref_catch_w_cd(99));
    ref_fcatch_w_sd(kk) = sum(ref_catch_w_sd(99));  
    ref_Hr(kk) = ref_Hrt(99);	 
  }
   i35 = 0; i40 = 0;
  for (j = 1; j<= 201; j++)
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
  report1 <<"Harvest rate in term of legal males as F = 0.00, 0.01, ... 1.0"<<endl;
  report1 << ref_Hr<<endl;
  report1 <<"Total MMB (>93mm) Feb 01 as F = 0.00, 0.01, ... 1.0 (lbs/R)"<<endl;
  report1 << ref_mbio<<endl;
  report1 <<"Total legal as F = 0.00, 0.01, ... 1.0 (lbs/R)"<<endl;
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
  report << "nps" << endl << trans(nps) << endl;  // Modeled new shell summer 
  report << "ops" << endl << trans(ops) << endl;  // Modeled old shell summer
  report << "npw" << endl << trans(npw) << endl;  // Modeled new shell winter
  report << "opw" << endl << trans(opw) << endl;  // Modeled oldshell shell winter  
  report << "ett" << endl << ett << endl;  // Estimated trawl abundance
  report << "ecpue" << endl << ecpue << endl;  // Estimated Summer fishery cpue
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
