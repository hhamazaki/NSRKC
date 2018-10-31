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
#include <NS3n2016_Feb1c3.htp>

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
  lamc.allocate("lamc");
  lamw.allocate("lamw");
  lawp.allocate("lawp");
  latag.allocate("latag");
  st_est.allocate("st_est");
  nst.allocate("nst");
  nsc.allocate("nsc");
 cout << "Data Section Completed" << endl;
 cout << "lg " << lg << endl;
 cout << "nono " << endl << nono << endl;
 cout << "nst " << endl << nst << endl;  
 cout << "nsc " << endl << nsc << endl;
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
  flnp.allocate(1,na-1,-10.0,10.0,1,"flnp");
  rlnp.allocate(1,ra-1,-15.0,15.0,1,"rlnp");
  log_mo1.allocate(-10.5,-1.0,1,"log_mo1");
  log_st1.allocate(1,2,-15.0,1.0,st_est,"log_st1");
  log_sw1.allocate(-15.0,1.0,1,"log_sw1");
  sw3.allocate(1,ra,0.,1.0,1,"sw3");
  log_sc1.allocate(1,2,-15.0,1.0,1,"log_sc1");
  advar.allocate(0.0,6.0,1,"advar");
  qtno.allocate(0.1,1.0,1,"qtno");
  M.allocate(0.02,1.0,M_phase,"M");
  ms.allocate(1.0,5.0,ms_phase,"ms");
  sigma.allocate(0.,30.,2,"sigma");
  ig.allocate(1,2,0.,20.,1,"ig");
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
  log_rec.allocate(fyear,lyear-1,"log_rec");
  #ifndef NO_AD_INITIALIZE
    log_rec.initialize();
  #endif
  nps.allocate(1,na,fyear,lyear+1,"nps");
  #ifndef NO_AD_INITIALIZE
    nps.initialize();
  #endif
  ops.allocate(1,na,fyear,lyear+1,"ops");
  #ifndef NO_AD_INITIALIZE
    ops.initialize();
  #endif
  npw.allocate(1,na,fyear,lyear+1,"npw");
  #ifndef NO_AD_INITIALIZE
    npw.initialize();
  #endif
  opw.allocate(1,na,fyear,lyear+1,"opw");
  #ifndef NO_AD_INITIALIZE
    opw.initialize();
  #endif
  enc.allocate(1,na,fyear,lyear,"enc");
  #ifndef NO_AD_INITIALIZE
    enc.initialize();
  #endif
  eoc.allocate(1,na,fyear,lyear,"eoc");
  #ifndef NO_AD_INITIALIZE
    eoc.initialize();
  #endif
  ecpue.allocate(fyear,lyear,"ecpue");
  #ifndef NO_AD_INITIALIZE
    ecpue.initialize();
  #endif
  cvcpue.allocate(fyear,lyear,"cvcpue");
  #ifndef NO_AD_INITIALIZE
    cvcpue.initialize();
  #endif
  enwp.allocate(1,na,fyear,lyear,"enwp");
  #ifndef NO_AD_INITIALIZE
    enwp.initialize();
  #endif
  eowp.allocate(1,na,fyear,lyear,"eowp");
  #ifndef NO_AD_INITIALIZE
    eowp.initialize();
  #endif
  enwpd.allocate(1,na,fyear,lyear,"enwpd");
  #ifndef NO_AD_INITIALIZE
    enwpd.initialize();
  #endif
  eowpd.allocate(1,na,fyear,lyear,"eowpd");
  #ifndef NO_AD_INITIALIZE
    eowpd.initialize();
  #endif
  enwc.allocate(1,na,fyear,lyear,"enwc");
  #ifndef NO_AD_INITIALIZE
    enwc.initialize();
  #endif
  eowc.allocate(1,na,fyear,lyear,"eowc");
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
  tw.allocate(fyear,lyear,"tw");
  #ifndef NO_AD_INITIALIZE
    tw.initialize();
  #endif
  twp.allocate(fyear,lyear,"twp");
  #ifndef NO_AD_INITIALIZE
    twp.initialize();
  #endif
  rec.allocate(fyear,lyear,"rec");
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
  eno.allocate(1,na,1,nyo,"eno");
  #ifndef NO_AD_INITIALIZE
    eno.initialize();
  #endif
  eoo.allocate(1,na,1,nyo,"eoo");
  #ifndef NO_AD_INITIALIZE
    eoo.initialize();
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
  selc1.allocate(1,2,1,na,"selc1");
  #ifndef NO_AD_INITIALIZE
    selc1.initialize();
  #endif
  sel.allocate(1,na,fyear,lyear,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  st1.allocate(1,2,"st1");
  #ifndef NO_AD_INITIALIZE
    st1.initialize();
  #endif
  sc1.allocate(1,2,"sc1");
  #ifndef NO_AD_INITIALIZE
    sc1.initialize();
  #endif
  selt1.allocate(1,2,1,na,"selt1");
  #ifndef NO_AD_INITIALIZE
    selt1.initialize();
  #endif
  selt.allocate(1,na,fyear,lyear,"selt");
  #ifndef NO_AD_INITIALIZE
    selt.initialize();
  #endif
  selw.allocate(1,na,"selw");
  #ifndef NO_AD_INITIALIZE
    selw.initialize();
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
  FSum.allocate(fyear,lyear,"FSum");
  #ifndef NO_AD_INITIALIZE
    FSum.initialize();
  #endif
  FWin.allocate(fyear,lyear,"FWin");
  #ifndef NO_AD_INITIALIZE
    FWin.initialize();
  #endif
  FSub.allocate(fyear,lyear,"FSub");
  #ifndef NO_AD_INITIALIZE
    FSub.initialize();
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
  tt1.allocate(1,nyt,"tt1");
  #ifndef NO_AD_INITIALIZE
    tt1.initialize();
  #endif
  tt2.allocate(1,nyt,"tt2");
  #ifndef NO_AD_INITIALIZE
    tt2.initialize();
  #endif
  tf.allocate(1,13,"tf");
  #ifndef NO_AD_INITIALIZE
    tf.initialize();
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
  legaln.allocate(fyear,lyear,"legaln");
  legalb.allocate(fyear,lyear,"legalb");
  mmb0.allocate(fyear,lyear+1,"mmb0");
  legal.allocate(fyear,lyear+1,"legal");
  #ifndef NO_AD_INITIALIZE
    legal.initialize();
  #endif
  mmb.allocate(fyear,lyear+1,"mmb");
  #ifndef NO_AD_INITIALIZE
    mmb.initialize();
  #endif
  last_legal.allocate("last_legal");
  last_subl.allocate("last_subl");
  bmsy.allocate("bmsy");
  #ifndef NO_AD_INITIALIZE
  bmsy.initialize();
  #endif
  last_mmb.allocate("last_mmb");
  last_ofl.allocate("last_ofl");
 cout << "Parameter Section Completed" << endl;
}

void model_parameters::initializationfunction(void)
{
  log_q.set_initial_value(-6.5);
  log_mo1.set_initial_value(-2.41358);
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
  dvariable tt0,n0; // Calculated working variable does not need defferentiated
  M = M2;
  ms = ms6;
  tt0.initialize();
  n0.initialize();
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
	
  for (i=1;i<=nyt;i++)
   {
    st(i) *= efn(1);
    if (st(i) > maxs(1)) st(i)=maxs(1);
   }
 cout << st << endl;
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
 cout << twsd << endl;
   
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
          
  for(i=1;i<=ntag1;i++)
  {
        tagrecap1(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,3);
        tagrecap2(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,4);
        tagrecap3(tag_recov1(i,1),tag_recov1(i,2))=tag_recov1(i,5); 
  }
 
  for(i=1;i<=ntag2;i++)
  {
        tagrecap12(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,3);
        tagrecap22(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,4);
        tagrecap32(tag_recov2(i,1),tag_recov2(i,2))=tag_recov2(i,5); 
  }
  cout << tagrecap1 << endl;
  cout << tagrecap12 << endl;
  
  if(nsc == 1)
  {
  tagrecap1 = tagrecap1+tagrecap12;
  tagrecap2 = tagrecap2+tagrecap22;
  tagrecap3 = tagrecap3+tagrecap32;
  }
  
  for(i=1;i<=na;i++)
  {  
  ptagrecap1(i)=tagrecap1(i)/(rowsum(tagrecap1)(i)+0.0000001);
  ptagrecap2(i)=tagrecap2(i)/(rowsum(tagrecap2)(i)+0.0000001);
  ptagrecap3(i)=tagrecap3(i)/(rowsum(tagrecap3)(i)+0.0000001);
  ptagrecap12(i)=tagrecap12(i)/(rowsum(tagrecap12)(i)+0.0000001);
  ptagrecap22(i)=tagrecap22(i)/(rowsum(tagrecap22)(i)+0.0000001);  
  ptagrecap32(i)=tagrecap32(i)/(rowsum(tagrecap32)(i)+0.0000001);
  }
     
 
 cout << tagrecap1 << endl;
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
  dvariable mo1,sw1;    
  double pp, L,Lm, Ls; 
  npp.initialize();				
  rpp.initialize();				
   q = mfexp(log_q);
   qw = mfexp(log_qw);   
   mo1 = mfexp(log_mo1);
   sc1 = mfexp(log_sc1);
   st1 = mfexp(log_st1); 
   sw1 = mfexp(log_sw1);
    Lm = slm + slt*(na-1);
 // calculate molting and selectivity functions   
   for (j=1;j<=na;j++)
    {	
     L = (double(j)-1.0)*slt;
     L += slm;
     mp1(j) = 1.0-1.0/(1.0+mfexp(mo1*(slm-L)+log(1.0/0.001-1.0)));  // molting probability : reverse logistic function//
     selc1(1,j) = 1.0/(1.0+mfexp(sc1(1)*(Lm-L)+log(1.0/0.999-1.0)));  // commercial pot selectivity 1976-1992: logistic function//
     selc1(2,j) = 1.0/(1.0+mfexp(sc1(nsc)*(Lm-L)+log(1.0/0.999-1.0)));  // commercial pot selectivity 1993 - : logistic function//
     selw(j) = 1.0-1.0/(1.0+mfexp(sw1*(slm-L)+log(1.0/0.001-1.0)));  // winter pot selectivity 1981-2011: reverse logistic function//
     selt1(1,j) = 1.0/(1.0+mfexp(st1(1)*(Lm-L)+log(1.0/0.999-1.0)));   // NOAA summer trawl net selectivity 1976-1992: logistic function//
     selt1(2,j) = 1.0/(1.0+mfexp(st1(nst)*(Lm-L)+log(1.0/0.999-1.0)));   // ADFG summer trawl net selectivity: 1996- logistic function//     
    }
  for (j=1;j<=msn(2);j++) selw(j) = sw3(j);
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(1))    // period 1976 - 1992 
      for (j=1; j<=na; j++) selt(j,i) = selt1(1,j);   //NOAA trawl survey 
      else                          
      for (j=1; j<=na; j++) selt(j,i) = selt1(1,j);   //ADFG trawl survey
    }
   for (i=fyear;i<=lyear;i++)
    {
     if(i<=scy(1))    // period 1976 - 1992 
      for (j=1; j<=na; j++) sel(j,i) = selc1(1,j); 
      else                          // period 2008 - present
      for (j=1; j<=na; j++) sel(j,i) = selc1(2,j);
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
		npw(j,fyear) = npp(j)*first_y;
		opw(j,fyear) = 0.0;
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
    Mn = M;
  for (i=1;i<=msn(1);i++)
    { 
    Mn(na+1-i) = ms*M;
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
  for (i=fyear;i<=lyear;i++)
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
    for (j=ra+1;j<=na;j++)
     {
      pp1 += (npw(j,i)+opw(j,i))*selw(j);  // For subsistence catch   
     }
    for (j=ra+1;j<=na;j++)
     {
	  enwp(j,i) = npw(j,i)*selw(j)/pp1;
      eowp(j,i) = opw(j,i)*selw(j)/pp1;
     }
	pp2 = 0.0;
    for (j=1;j<=ra;j++)
     {
      pp2 += (npw(j,i)+opw(j,i))*selw(j);  // For subsistence catch   
     }
    for (j=1;j<=ra;j++)
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
     ops(j,i) = (opw(j,i)-eowp(j,i)*tws(i)-eowc(j,i)*twc(i)-(twc(i)/tt1)*(opw(j,i))*selw(j)*(1-lg(j))*hm(2)-twsd(i)*eowpd(j,i)*hm(2))*exp(-0.417*Mn(j));	
	 if (nps(j,i) < 0.0) nps(j,i) = 0.001; // summer abundance should not go bellow zero
	 if (ops(j,i) < 0.0) ops(j,i) = 0.001; // summer abundance should not go bellow zero	
     bcw(i) += (twc(i)/tt1)*(npw(j,i)+opw(j,i))*selw(j)*(1-lg(j))*hm(2); 
	}
	tb(i) = 0.0; 
    for (j=1;j<=na;j++)
	{
        tb(i) += (nps(j,i)+ops(j,i))*sel(j,i)*lg(j); //total summer crab abundance available to 
	}		
	// mean exploitable leagal abundance does not go negative
    if (tb(i) < 0.001) tb(i) = 0.001;
	// Calculate proprotion of newshell and oldshell by comm fiish 	
	for (j=1;j<=na;j++)
     {
        enc(j,fyear) = 0.0;
        eoc(j,fyear) = 0.0;
	 	enc(j,i) = nps(j,i)*sel(j,i)*lg(j)/tb(i);
		eoc(j,i) = ops(j,i)*sel(j,i)*lg(j)/tb(i);
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
      npw(j,i+1) = pp*exp(-(0.583-ys(i))*Mn(j)) + rpp(j)*rec(i);   //New shell crab molted + recruits 
      opw(j,i+1) = tt0(j)*(1.0-mp0(j))*exp(-(0.583-ys(i))*Mn(j));  //Old shell crab are unmolted crab  
     }
   }
 for (i=fyear; i<=lyear; i++)
   {
    legal(i) = 0.0; legalb(i) = 0.0;
    for (j=1; j<=na; j++)
    {
       legal(i) += (npw(j,i)+opw(j,i))*lg(j);
       legalb(i) += (npw(j,i)+opw(j,i))*lg(j)*wm(j);
    }
    mmb(i) = 0.0;
    for (j=ra+1; j<=na; j++) mmb(i) += (npw(j,i)+opw(j,i))*wm(j);
    mmb0(i) = mmb(i);    
    legaln(i) = legal(i);
	}
  last_y = column(npw,lyear+1);
  for (j=1; j<=na; j++) last_y(j) += opw(j,lyear+1);
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
      ent0(j,i) = ((nps(j,it(i)))*exp(-bf*Mn(j))-(enc(j,it(i)))*pct(i)*tc(it(i)))*exp(-af*Mn(j))*selt(j,it(i));
      eot0(j,i) = ((ops(j,it(i)))*exp(-bf*Mn(j))-(eoc(j,it(i)))*pct(i)*tc(it(i)))*exp(-af*Mn(j))*selt(j,it(i)); 
	  if (ent0(j,i) < 0.0) ent0(j,i) = 0.0;
	  if (eot0(j,i) < 0.0) eot0(j,i) = 0.0;	  
      ett(i) += ent0(j,i)+eot0(j,i);
      if (j<=sla) ett1(i) += ent0(j,i)+eot0(j,i);
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
  if(it(i) < qyear) {ettq(i) = qtno*ett(i);}
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
      ewcpue(i) = qw*twp(iw(i)); 	 
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
  for (i=fyear;i<=lyear;i++)
   {
   if (stcpue(i) <= 0.0) {ecpue(i) = 0.0;}
   }
}

void model_parameters::evaluate_the_objective_function(void)
{
  tf(1) = 0.5*norm2(elem_div((log(tt+1.e-3)-log(ettq+1.e-3)),sqrt(log(elem_prod(cv,cv)+1.0))));
  T_var = sqrt(log(elem_prod(cvcpue,cvcpue)+1.0)+advar); 
  tf(3) = lamc*(sum(log(T_var))+0.5*(norm2(elem_div((log(stcpue+1.e-3)-log(ecpue+1.e-3)),T_var))));  
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
  ef(4) = -(nsc-1)*(sum(elem_prod(rowsum(tagrecap12),rowsum(elem_prod(ptagrecap12,log(egr12+1.e-3))))-elem_prod(rowsum(tagrecap12),rowsum(elem_prod(ptagrecap12,log(ptagrecap12+1.e-3)))))); 
  ef(5) = -(nsc-1)*(sum(elem_prod(rowsum(tagrecap22),rowsum(elem_prod(ptagrecap22,log(egr22+1.e-3))))-elem_prod(rowsum(tagrecap22),rowsum(elem_prod(ptagrecap22,log(ptagrecap22+1.e-3))))));  
  ef(6) = -(nsc-1)*(sum(elem_prod(rowsum(tagrecap32),rowsum(elem_prod(ptagrecap32,log(egr32+1.e-3))))-elem_prod(rowsum(tagrecap32),rowsum(elem_prod(ptagrecap32,log(ptagrecap32+1.e-3))))));  
  tf(13) = latag*sum(ef);
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
