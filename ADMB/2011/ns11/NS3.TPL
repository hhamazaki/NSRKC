DATA_SECTION
  init_int ny                 //number of years
  init_int na                 //number of length groups
  init_int nyt                //number of years with trawl survey
  init_int nyp                //number of years with pot survey
  init_int nyw                //number of years with winter project
  init_int nyf                //number of years with pre-fishery survey
  init_int nyo                //number of years with observer's data
  init_number slm             //mid-size of smallest length group (mm)
  init_number slt             //length interval (mm)
  init_number M               //instantaneous natural mortality
  init_number lamt            //weighting factor for trawl survey
  init_number lamp            //weighting factor for pot survey
  init_number lamc            //weighting factor for fishing effort: summer fishery
  init_number maxss           //maximum effective sample size for length proportion
  init_number efn             //% of effective sample size due to multin.distribution
  init_vector lg(1,na)        //proportion of legals by length group
  init_matrix gr(1,na,1,na)   //growth matrix
  init_ivector it(1,nyt)      //year id for trawl survey data
// !! cout <<it<<endl;
  init_vector tt(1,nyt)       //total annual abundance from trawl survey
  init_vector pct(1,nyt)      //% of summer catch occurred before the mid point of survey
  init_vector yt(1,nyt)       //Mid point of survey from July 1
// !! cout <<yt<<endl;
  init_vector st(1,nyt)       //annual sample size from trawl survey
// !! cout <<st<<endl;
  init_vector cv1(1,nyt)       //annual cv of sublegals for trawl survey
  init_vector cv2(1,nyt)       //annual cv of legals for trawl survey
  init_vector cv3(1,nyt)       //annual cv of totals for trawl survey
  init_ivector ip(1,nyp)      //year id for pot survey data
  init_vector tp(1,nyp)       //total annual abundance from pot survey
  init_vector yp(1,nyp)       //Mid point of survey from July 1
// !! cout <<yp<<endl;
  init_vector sp(1,nyp)       //annual sample size from pot survey
  init_ivector iw(1,nyw)      //year id for winter project
  init_vector sw(1,nyw)       //annual sample size from winter project
  init_int ipre               //year id for pre-fishery survey
  init_number spre            //annual sample size from pre-fishery survey
  init_vector sc(1,ny)        //annual sample size from summer fishery
  init_vector tc(1,ny)        //annual catch for summer fishery
  init_vector te(1,ny)        //annual effort for summer fishery
// !! cout <<te<<endl;
  init_vector ys(1,ny)        //Mid point of summer fishery from July 1
  init_ivector io(1,nyo)      //year id for observer's data
  init_vector so(1,nyo)       //annual sample size from observer's data
  init_vector twc(1,ny-1)     //annual catch for winter fishery
  init_vector tsc(1,ny-1)     //annual catch for subsistence fishery
// !! cout <<tsc<<endl;
  init_matrix ont(1,na,1,nyt) //proportion of length group for new-shell, trawl survey
  init_matrix oot(1,na,1,nyt) //proportion of length group for old-shell, trawl survey
  init_matrix onp(1,na,1,nyp) //proportion of length group for new-shell, pot survey
  init_matrix oop(1,na,1,nyp) //proportion of length group for old-shell, pot survey
  init_matrix onw(1,na,1,nyw) //proportion of length group for new-shell, winter survey
  init_matrix oow(1,na,1,nyw) //proportion of length group for old-shell, winter survey
  init_vector onf(1,na)       //proportion of length group for new-shell, pre-fishery survey
  init_vector oof(1,na)       //proportion of length group for old-shell, pre-fishery survey
//!! cout <<oof<<endl;
  init_matrix onc(1,na,1,ny)  //proportion of length group for new-shell, summer fishery
  init_matrix ooc(1,na,1,ny)  //proportion of length group for old-shell, summer fishery
  init_matrix ono(1,na,1,nyo) //proportion of length group for new-shell, observer's data
  init_matrix ooo(1,na,1,nyo) //proportion of length group for old-shell, observer's data
  init_vector wm(1,na)        //Mean weight
//!! cout <<wm<<endl;
//  !!exit(1);

INITIALIZATION_SECTION
  log_q1         -11.5
  log_q2         -11.5
//M                0.30
//MM               1.50
//log_recscale     6.0
  log_initpop      7.0
  r1               0.6
  log_mo1         -2.41358
  log_mo2          4.8211
  log_mo3         -2.41358
  log_mo4          4.82115
  log_mo5         -2.41358
  log_mo6          4.82115
//  st1              0.7
//  st2              0.85
//  sp1              0.55
//  sp2              0.80
//  sw1              0.5
//  sw2              0.80
  sw3              0.8
  log_sc1         -2.31267
  log_sc2          4.485
  sos0             1.0
  sow              1.00

// TOP_OF_MAIN_SECTION
//  gradient_structure::set_NO_DERIVATIVES();

PARAMETER_SECTION
  init_bounded_number log_q1(-12.5,-8.5,1)
  init_bounded_number log_q2(-12.5,-10.0,1)
//init_bounded_number M(0.02,0.50,2)
//init_bounded_number MM(1.0,2.0,2)
  init_number log_initpop(1)
//init_number log_recscale(1)
//init_bounded_dev_vector log_relrec(1,ny-1,-15.,15.,1)
  init_bounded_dev_vector log_relpop(1,ny,-15.,15.,1)
  init_bounded_number r1(0.5,0.9,1)         //proportion of recruits to length group 1
  init_bounded_number log_mo1(-4.15,-2.0,2) //parameter for molting probability
  init_bounded_number log_mo2(3.0,5.00,2)   //length at 50% molting probability
//  init_bounded_number log_mo3(-3.15,-2.0,2) //parameter for molting probability
//  init_bounded_number log_mo4(2.3,4.95,2)   //length at 50% molting probability
  init_bounded_number log_mo5(-2.85,-1.6,-2) //parameter for molting probability
  init_bounded_number log_mo6(4.0,4.95,-2)   //length at 50% molting probability
//  init_bounded_number st1(0.30,0.9,1)      //selectivity for length group 1: trawl survey
//  init_bounded_number st2(0.50,1.0,1)       //selectivity for length group 2: trawl survey
//  init_bounded_number sp1(0.30,0.8,1)      //selectivity for length group 1: pot survey
//  init_bounded_number sp2(0.80,1.0,1)       //selectivity for length group 2: pot survey
//  init_bounded_number sw1(0.25,1.00,1)      //selectivity for length group 1: winter project
//  init_bounded_number sw2(0.60,1.0,1)       //selectivity for length group 2: winter project
  init_bounded_number log_st1(-5.5,-1.0,1)      //parameter for selectivity: trawl survey
  init_bounded_number log_st2(0.51,3.00,1)      //length at 50% selectivity: trawl survey
  init_bounded_number log_sp1(-3.7,-1.0,1)      //parameter for selectivity: pot survey
  init_bounded_number log_sp2(2.00,5.00,1)      //length at 50% selectivity: pot survey
  init_bounded_number log_sw1(-3.0,-0.3,1)      //parameter for selectivity: winter project
  init_bounded_number log_sw2(3.00,5.00,1)      //length at 50% selectivity: winter project
  init_bounded_number sw3(0.10,1.0,1)       //selectivity for length group NA: winter project
  init_bounded_number log_sc1(-4.5,-1.9,1)  //parameter for selectivity: summer fishery before 1993
  init_bounded_number log_sc2(4.35,6.75,1)  //length at 50% selectivity: summer fishery before 1993
  init_bounded_number log_sc3(-3.5,-1.9,1) //parameter for selectivity: summer fishery after 1992
  init_bounded_number log_sc4(4.35,6.75,1)   //length at 50% selectivity: summer fishery after 1992
//  init_bounded_number s931(0.0,0.5,1) //parameter for molting probability
//  init_bounded_number s932(0.3,1.0,1)   //length at 50% molting probability
//  init_bounded_number s933(0.2,1.0,1)   //length at 50% molting probability
  init_bounded_number sos0(0.50,1.5,-1)      //selectivity of old-shell crabs in summer commercial fishery
  init_bounded_number sow(0.300,1.7,-1)       //selectivity of old-shell crabs in winter project

//vector log_rec(1,ny-1)
  vector log_pop(1,ny)
  matrix nps(1,na,1,ny)        //new-shell population abundance in summer
  matrix ops(1,na,1,ny)        //old-shell population abundance in summer
  matrix npw(1,na,1,ny-1)      //new-shell population abundance in winter
  matrix opw(1,na,1,ny-1)      //old-shell population abundance in winter
  vector rec(1,ny-1)           //annual recruits: starting in year 1.
  matrix ent(1,na,1,nyt)       //estimated proportion: new shell,trawl survey
  matrix et0(1,na,1,nyt)       //estimated trawl survey abundance
  matrix eot(1,na,1,nyt)       //estimated proportion: old shell,trawl survey
  matrix enp(1,na,1,nyp)       //estimated proportion: new shell,pot survey
  matrix eop(1,na,1,nyp)       //estimated proportion: old shell,pot survey
  matrix enw(1,na,1,nyw)       //estimated proportion: new shell,winter project
  matrix eow(1,na,1,nyw)       //estimated proportion: old shell,winter project
  vector enf(1,na)             //estimated proportion: new shell,pre-fishery survey
  vector eof(1,na)             //estimated proportion: old shell,pre-fishery survey
  matrix enc(1,na,1,ny)        //estimated proportion: new shell,summer fishery
  matrix eoc(1,na,1,ny)        //estimated proportion: old shell,summer fishery
  matrix eno(1,na,1,nyo)       //estimated proportion: new shell,observer's data
  matrix eoo(1,na,1,nyo)       //estimated proportion: old shell,observer's data
  matrix en0(1,na,1,ny-1)      //estimated proportion: new shell,subsistence fishery
  matrix eo0(1,na,1,ny-1)      //estimated proportion: old shell,subsistence fishery
  matrix en1(1,na,1,ny-1)      //estimated proportion: new shell,winter fishery
  matrix eo1(1,na,1,ny-1)      //estimated proportion: old shell,winter fishery
  vector ett(1,nyt)            //estimated total annual abundance: trawl survey
  vector ett1(1,nyt)           //estimated annual non-legal crab abundance: trawl survey
  vector ett2(1,nyt)           //estimated annual legal crab abundance: trawl survey
  vector tt1(1,nyt)            //observed annual non-legal crab abundance: trawl survey
  vector tt2(1,nyt)            //observed annual legal crab abundance: trawl survey
  vector etp(1,nyp)            //estimated total annual abundance: pot survey
  vector etp1(1,nyp)           //estimated annual non-legal crab abundance: pot survey
  vector etp2(1,nyp)           //estimated annual legal crab abundance: pot survey
  vector tp1(1,nyp)            //observed annual non-legal crab abundance: pot survey
  vector tp2(1,nyp)            //observed annual legal crab abundance: pot survey
  vector ete(1,ny)             //estimated total annual effort: summer fishery
  vector tb(1,ny)              //estimated mean annual summer exploitable abundance
  vector mp1(1,na)             //molting probability  1976-1981
  vector mp2(1,na)             //molting probability  1982-1990
  vector mp3(1,na)             //molting probability: 1991-1997
  vector mp0(1,na)             //working vector: molting probability
  vector sel92(1,na)           //selectivity of summer fishery before 1993
  vector sel93(1,na)           //selectivity of summer fishery after 1992
  matrix sel(1,na,1,ny)
  vector selt(1,na)            //selectivity of trawl survey
  vector selp(1,na)            //selectivity of pot survey
  vector selw(1,na)            //working variable: selectivity of winter project.
  vector Mn(1,na)              //working variable: natural mortality for new-shell crabs.
  vector Mo(3,na)              //working variable: natural mortality for old-shell crabs.
  vector sos(1,na)             //working variable: additional selectivity for old-shell.
  number q1                    //catchability: summer fishery
  number q2                    //catchability: summer fishery
  number pre                   //working variable:
  number first_y               //working variable: first-year abundance
  number mo1                   //working variable: molting probability
  number mo2                   //working variable: molting probability
  number mo3                   //working variable: molting probability
  number mo4                   //working variable: molting probability
  number mo5                   //working variable: molting probability
  number mo6                   //working variable: molting probability
  number sc1                   //working variable: selectivity
  number sc2                   //working variable: selectivity
  number st1                   //working variable: selectivity
  number st2                   //working variable: selectivity
  number sp1                   //working variable: selectivity
  number sp2                   //working variable: selectivity
  number sw1                   //working variable: selectivity
  number sw2                   //working variable: selectivity
  number sc10                  //working variable: selectivity
  number sc20                  //working variable: selectivity
  number pp                    //working variable
  number tt0                   //working variable
  number tf1                   //working variable
  number tf2                   //working variable
  number tf3                   //working variable
  number tf4                   //working variable
  number tf5                   //working variable
  number tf6                   //working variable
  number tf7                   //working variable
  number tf8                   //working variable
  number tf9                   //working variable
  number tf10                  //working variable
  number tf11                  //working variable
  number tf12                  //working variable
  number tf13                  //working variable
  number tf14                  //working variable
  number tf15                  //working variable

  matrix ref_catch_s(1,na,1,100)
  matrix ref_catch_w_n(1,na,1,100)
  matrix ref_catch_w_o(1,na,1,100)
  matrix ref_catch(1,na,1,100)
  matrix ref_ps(1,na,1,100)
  matrix ref_nps(1,na,1,100)
  matrix ref_ops(1,na,1,100)
  matrix ref_npw(1,na,1,100)
  matrix ref_opw(1,na,1,100)
  vector ref_catch_b(1,100)
  vector ref_catch_n(1,100)
  number ref_F
  number ref_Hwc
  number ref_Hts
  vector ref_Htw(1,ny)
  vector ref_mbio(1,101)
  vector ref_mmb(1,100)
  vector ref_totc(1,101)
  vector ref_legaln(1,101)
  vector ref_legalb(1,101)
  vector ref_Hr(1,101)
  vector ref_legal(1,100)
  vector ref_legal_b(1,100)
  number ref_ys
  number f35
  number f40
  number h35
  number h40
  number i35
  number i40

  objective_function_value f
  sdreport_vector last_y(1,na)
  sdreport_vector legal(1,ny)
  sdreport_vector mmb(1,ny)
  likeprof_number last_legal
  likeprof_number last_mmb

PRELIMINARY_CALCS_SECTION
  Mo = M;
//  Mo(na-1) = 1.00*M;
//  Mo(na) = 1.2*M;
  Mn = M;
//  Mn(na-1) = 1.00*M;
//Mn(na) = 1.6*M;
//  Mn(na) = 1.2*M;
  sos = 1.0;
  selw = 1.0;
  selt = 1.0;
  selp = 1.0;
  for (int i=1;i<=ny;i++)
  {
    for (int j=1;j<3;j++)
    {
      en0(j,i) = 0.0;
      eo0(j,i) = 0.0;
      en1(j,i) = 0.0;
      eo1(j,i) = 0.0;
    }
  }
//set the max effective sample size for length composition.
//all data have effective sample size reduced to 50%
  for (i=1;i<=nyt;i++)
  {
    st(i) = st(i)*0.5;
    if (st(i) > maxss) st(i)=maxss*efn;             //trawl survey
  }
  for (i=1;i<=nyp;i++)
  {
    sp(i) = sp(i)*0.5;
    if (sp(i) > maxss) sp(i)=maxss*efn;             //pot survey
  }
  for (i=1;i<=nyw;i++)
  {
    sw(i) = sw(i)*0.10;
    if (sw(i) > maxss) sw(i)=maxss*efn;             //winter project
  }
//  pre = 0.10*spre;
    pre = 0.0;                                      //do not use pre-fishery data (one year of data only)
  if (pre > maxss) pre=maxss*efn;                   //pre-fishery
  for (i=1;i<=ny;i++)
  {
    sc(i) = sc(i)*0.10;
    if (sc(i) > maxss) sc(i)=maxss*efn;             //summer fishery
  }
  for (i=1;i<=nyo;i++)
  {
    so(i) = so(i)*0.10;
    if (so(i) > maxss) so(i)=maxss*efn;             //observer's data
  }
//combine new- and old-shell data for all trawl survey data & compute observed legal & non-legal abundance
  for (i=1;i<=nyt;i++)
  {
    for (int j=3;j<=na;j++)
    {
       ont(j,i) = ont(j,i) + oot(j,i);
    }
    tt1(i) = tt(i)*(ont(1,i)+ont(2,i)+ont(3,i));
    tt2(i) = tt(i) - tt1(i);
  }
  for (i=1;i<=nyp;i++)
  {
    tp1(i) = tp(i)*(onp(1,i)+onp(2,i)+onp(3,i)+oop(3,i));
    tp2(i) = tp(i) - tp1(i);
  }


PROCEDURE_SECTION
  convert_parameters_into_rates();
//  cout <<"OK for convert_parameters..."<<endl;

  get_first_year_abundance();
//  cout <<"OK for get_first_year..."<<endl;

  get_number_by_size();
//  cout <<"OK for get_number_by_size..."<<endl;

  get_proportion_and_effort();
//  cout <<"OK for get_proportion_and..."<<endl;

  if (sd_phase()) get_reference_points();

  evaluate_the_objective_function();

FUNCTION convert_parameters_into_rates

//Mo(na) = MM*M;
//Mn(na) = MM*M;
//catchability
   q1 = mfexp(log_q1);
   q2 = mfexp(log_q2);
//molting probability and selectivity for summer commercial fisheries
   mo1 = mfexp(log_mo1);
   mo2 = mfexp(log_mo2);
   sc3 = mfexp(log_sc3);
   sc4 = mfexp(log_sc4);
   mo5 = mfexp(log_mo1);
   mo6 = mfexp(log_mo2);
   sc1 = mfexp(log_sc1);
   sc2 = mfexp(log_sc2);
   st1 = mfexp(log_st1);
   st2 = mfexp(log_st2);
   sp1 = mfexp(log_sp1);
   sp2 = mfexp(log_sp2);
   sw1 = mfexp(log_sw1);
   sw2 = mfexp(log_sw2);
   for (int j=1;j<=na;j++)
   {
//   pp = slm + (double(j)-1.0)*slt;
     pp = (double(j)-1.0)*slt;
     pp = pp + slm;
     mp1(j) = 1.0-1.0/(1.0+mfexp(-1.0*mo1*(pp-mo2)));
 //    mp2(j) = 1.0-1.0/(1.0+mfexp(-1.0*mo3*(pp-mo4)));
     mp2(j) = 1.0-1.0/(1.0+mfexp(-1.0*mo1*(pp-mo2)));
     mp3(j) = 1.0-1.0/(1.0+mfexp(-1.0*mo5*(pp-mo6)));
     sel92(j) = 1.0/(1.0+mfexp(-1.0*sc1*(pp-sc2)));
     sel93(j) = 1.0/(1.0+mfexp(-1.0*sc3*(pp-sc4)));
   }
//    sel93(1) = 0.0; sel93(2)=0.0; sel93(3)= s931; sel93(4) = s932; sel93(5)=1.0; sel93(6)=s933;
   sel92 = sel92/sel92(na-1);
   if (sel92(na) > 1.0) sel92(na) = 1.0;
   sel93 = sel93/sel93(na-1);
   if (sel93(na) > 1.0) sel93(na) = 1.0;
   sel93(na) = sel93(na)*0.6;
   for (int i = 1; i<=ny; i++)
   {
     if(i<17)
     {
       for (int j=1; j<=na; j++) sel(j,i) = sel92(j);
     }
     else
     {
       for (int j=1; j<=na; j++) sel(j,i) = sel93(j);
     }
   }

// mp1(na) = 1.0-1.0/(1.0+mfexp(-1.0*mo1*(5.0+pp-mo2)));
// mp2(na) = 1.0-1.0/(1.0+mfexp(-1.0*mo3*(5.0+pp-mo4)));
// mp3(na) = 1.0-1.0/(1.0+mfexp(-1.0*mo5*(5.0+pp-mo6)));
//the first 2 length groups molt 100% and the last 2 length groups are 100% selected.
   mp1 = mp1/mp1(1);
   mp1(1) = 1.0;
   mp2 = mp2/mp2(1);
   mp2(1) = 1.0;
   mp3 = mp3/mp3(1);
   mp3(1) = 1.0;
   sos(na) = sos0;
   selw(1) = 1.0/(1.0+mfexp(-1.0*sw1*(78.5-sw2)));
   selw(2) = 1.0/(1.0+mfexp(-1.0*sw1*(88.5-sw2)));
   selw(na) = sw3;
//   if (last_phase())
//   {
//     if (st2<st1) st2 = st1;
//     if (sp2<sp1) sp2 = sp1;
//   }
   selt(1) = 1.0/(1.0+mfexp(-1.0*st1*(78.5-st2)));
   selt(2) = 1.0/(1.0+mfexp(-1.0*st1*(88.5-st2)));
   selp(1) = 1.0/(1.0+mfexp(-1.0*sp1*(78.5-sp2)));
   selp(2) = 1.0/(1.0+mfexp(-1.0*sp1*(88.5-sp2)));

FUNCTION get_first_year_abundance
// log_rec = log_relrec+log_recscale;
   log_pop = log_relpop+log_initpop;
// first_y = mfexp(log_initpop);
   first_y = mfexp(log_pop(ny));
   for (int i=1;i<ny;i++)
   {
//   rec(i) = mfexp(log_rec(i));
     rec(i) = mfexp(log_pop(i));
   }
// recruits in the last year are assumed to be average of the most recent five years
   rec(ny-1) = 0.2*(rec(ny-2)+rec(ny-3)+rec(ny-4)+rec(ny-5)+rec(ny-6));
//use the length proportion in the first year to estimate the abundance.
   tb(1) = 0.0;
   ops(1) = 0.0;
   ops(2) = 0.0;
   for (int j=1;j<=na;j++)
   {
     if (selt(j) < 0.1) selt(j) = 0.1;
     if (j <= 2)
     {
        nps(j,1) = ont(j,1)*first_y/selt(j);
     }
     else
     {

        nps(j,1) = (ont(j,1)-oot(j,1))*first_y/selt(j);
        ops(j,1) = oot(j,1)*first_y/selt(j);
     }
     tb(1) = tb(1) + (nps(j,1)+ops(j,1)*sos(j))*sel(j,1)*lg(j);
   }
   if (tb(1) < 0.001) tb(1) = 0.001;
//estimated length proportion for summer commercial fishery.
   for (j=1;j<=na;j++)
   {
     enc(j,1) = nps(j,1)*sel(j,1)*lg(j)/tb(1);
     eoc(j,1) = ops(j,1)*sos(j)*sel(j,1)*lg(j)/tb(1);
   }

FUNCTION get_number_by_size
   for (int i=1;i<ny;i++)
   {
//update to abundance in winter
//molting probability for 1993 is about 34% as high as in the other years.
       if (i < 7)
       {
        mp0 = mp1;
       }
       else if (i < 16)
       {
          mp0 = mp2;
       }
       else
       {
          mp0 = mp3;
       }
//   if (i==18)
//   {
//      mp0 = 0.34*mp3;
//      mp0 = 0.34*mp1;
//   }
     for (int j=1;j<=na;j++)
     {
        pp = 0;
        for (int k=1;k<=j;k++)
        {
//        tt0 = (nps(k,i)+ops(k,i)-tc(i)*(enc(k,i)+eoc(k,i)))*exp(-0.5*M);
          tt0 = ((nps(k,i)+ops(k,i))*exp(-ys(i)*Mn(k))-tc(i)*(enc(k,i)+eoc(k,i)))*exp(-(0.583-ys(i))*Mn(k));
          if (tt0<0.0) tt0 = 0.0;
          pp = pp + gr(k,j)*tt0*mp0(k);
        }
        npw(j,i) = pp;
        if (j==1) npw(j,i) = npw(j,i) + r1*rec(i);
        if (j==2) npw(j,i) = npw(j,i) + (1.0-r1)*rec(i);
//      opw(j,i) = (nps(j,i)+ops(j,i)-tc(i)*(enc(j,i)+eoc(j,i)))*exp(-0.5*M)*(1.0-mp0(j));
        opw(j,i) = ((nps(j,i)+ops(j,i))*exp(-ys(i)*Mn(j))-tc(i)*(enc(j,i)+eoc(j,i)))*exp(-(0.583-ys(i))*Mn(j))*(1.0-mp0(j));
        if (opw(j,i)<0.0) opw(j,i) = 0.0001;
        if (i==18)
        {
          npw(1,i) = npw(1,i) + opw(1,i);
          opw(1,i) = 0.0;
          npw(2,i) = npw(2,i) + opw(2,i);
          opw(2,i) = 0.0;
        }
     }
//propotions for winter fishery and subsistence fishery
     pp = 0.0;
     tt0 = 0.0;
     for (j=3;j<=na;j++)
     {
        pp = pp + (npw(j,i)+opw(j,i)*sow)*selw(j);
        tt0 = tt0 + (npw(j,i)+opw(j,i)*sow)*selw(j)*lg(j);
     }
     if (pp <= 0.0) pp = 0.000001;
     if (tt0 <= 0.0) tt0 = 0.000001;
     for (j=3;j<=na;j++)
     {
       en0(j,i) = npw(j,i)*selw(j)/pp;
       eo0(j,i) = opw(j,i)*sow*selw(j)/pp;
       en1(j,i) = npw(j,i)*selw(j)*lg(j)/tt0;
       eo1(j,i) = opw(j,i)*sow*selw(j)*lg(j)/tt0;
     }
//update to abundance in next summer
     for (j=1;j<=na;j++)
     {
//     nps(j,i+1) = (npw(j,i)-en0(j,i)*tsc(i)-en1(j,i)*twc(i))*exp(-0.417*M);
//     ops(j,i+1) = (opw(j,i)-eo0(j,i)*tsc(i)-eo1(j,i)*twc(i))*exp(-0.417*M);
       nps(j,i+1) = (npw(j,i)-en0(j,i)*tsc(i)-en1(j,i)*twc(i))*exp(-0.417*Mn(j));
       ops(j,i+1) = (opw(j,i)-eo0(j,i)*tsc(i)-eo1(j,i)*twc(i))*exp(-0.417*Mo(j));
     }
//proportion for summer fishery. Higher selectivity for new-shell than for old-shell.
     tb(i+1) = 0.0;
     for (j=1;j<=na;j++)
     {
       tb(i+1) = tb(i+1) + (nps(j,i+1)+ops(j,i+1)*sos(j))*sel(j,i+1)*lg(j);
     }
     if (tb(i+1) < 0.001) tb(i+1) = 0.001;
     for (j=1;j<=na;j++)
     {
       enc(j,i+1) = nps(j,i+1)*sel(j,i+1)*lg(j)/tb(i+1);
       eoc(j,i+1) = ops(j,i+1)*sos(j)*sel(j,i+1)*lg(j)/tb(i+1);
     }
   }
// if at the last_phase, calculated abundances in last year & legal abundance for all years.
   if (last_phase())
   {
     last_y = column(nps,ny);
     for (int j=3; j<=na; j++)
     {
        last_y(j) = last_y(j) + ops(j,ny);
     }

     for (int i=1; i<=ny; i++)
     {
        legal(i) = 0.0;
        for (j=1; j<=na; j++)
        {
           legal(i) = legal(i) + (nps(j,i)+ops(j,i))*lg(j);
        }
        mmb(i) = 0.0;
        for (j=3; j<=na; j++)
        {
           mmb(i) = mmb(i) + (nps(j,i)+ops(j,i))*wm(j);
        }
     }
   last_legal = last_y(3)*lg(3)*wm(3)+last_y(4)*lg(4)*wm(4)+last_y(5)*lg(5)*wm(5)+last_y(6)*lg(6)*wm(6);
   last_mmb   = mmb(ny);
   }

FUNCTION get_proportion_and_effort
   ett = 0.0; ett1 = 0.0;
   etp = 0.0; etp1 = 0.0; eop = 0.0;
   eow = 0.0; eoo =0.0; eof = 0.0;

//proprotion and total abundance for trawl survey
   for (int i=1;i<=nyt;i++)
   {
//     ett(i) = 0.0; ett1(i) = 0.0;
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
     for (int j=1;j<=na;j++)
     {
       et0(j,i) = ((nps(j,it(i))+ops(j,it(i)))*exp(-tf8*M)-(enc(j,it(i))+eoc(j,it(i)))*pct(i)*tc(it(i)))*exp(-tf10*M)*selt(j);
       if (et0(j,i) < 0.0) et0(j,i) = 0.0;
       ett(i) = ett(i) + et0(j,i);
       if (j<4) ett1(i) = ett1(i) + et0(j,i);
     }
     if (ett(i) <= 0.0) ett(i) = 0.000001;
     if (ett1(i) <= 0.0) ett1(i) = 0.000001;
     ett2(i) = ett(i) - ett1(i);
     for (j=1;j<=na;j++)
     {
//       ent(j,i) = (nps(j,it(i))+ops(j,it(i))*exp(-tf8*M)-(enc(j,it(i))+eoc(j,it(i)))*pct(i)*tc(it(i)))*exp(-tf10*M)*selt(j)/ett(i);
       ent(j,i) = et0(j,i)/ett(i);
//       if (j > 2)
//       {
//         eot(j,i) = (ops(j,it(i))*exp(-tf8*M)-eoc(j,it(i))*pct(i)*tc(it(i)))*exp(-tf10*M)*selt(j)/ett(i);
//combine new-shell and old-shell together for all survey data.
//         ent(j,i) = ent(j,i) + eot(j,i);
//         if (ent(j,i) < 0.0) ent(j,i) = 0.0;
//       }
       if (ent(j,i) < 0.0) ent(j,i) = 0.0;
     }
   }
    cout <<ett2<<endl;

//proprotion and total abundance for pot survey
   for (i=1;i<=nyp;i++)
   {
//     etp(i) = 0.0; etp1(i) = 0.0;
     for (int j=1;j<=na;j++)
     {
       etp(i) = etp(i) + (nps(j,ip(i))+ops(j,ip(i)))*exp(-yp(i)*M)*selp(j);
       if (j < 4) etp1(i) = etp1(i) + (nps(j,ip(i))+ops(j,ip(i)))*exp(-yp(i)*M)*selp(j);
     }
     if (etp(i) <= 0.0) etp(i) = 0.000001;
     if (etp1(i) <= 0.0) etp1(i) = 0.000001;
     etp2(i) = etp(i) - etp1(i);
     for (j=1;j<=na;j++)
     {
       enp(j,i) = nps(j,ip(i))*exp(-yp(i)*M)*selp(j)/etp(i);
       if (j > 1) eop(j,i) = ops(j,ip(i))*exp(-yp(i)*M)*selp(j)/etp(i);
     }
   }

//proprotion for winter project, higher selectivity for old-shell than for new-shell.
   for (i=1;i<=nyw;i++)
   {
     pp = 0.0;
     for (int j=1;j<=na;j++)
     {
       pp = pp + (npw(j,iw(i))+opw(j,iw(i))*sow)*selw(j);
     }
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
     for (int j=1;j<=na;j++)
     {
       pp = pp + (nps(j,io(i))+ops(j,io(i)))*sel(j,io(i))*(1.0-lg(j));
     }
     if (pp <= 0.0) pp = 0.000001;
     for (j=1;j<=na;j++)
     {
       eno(j,i) = nps(j,io(i))*sel(j,io(i))*(1.0-lg(j))/pp;
       if (j > 1) eoo(j,i) = ops(j,io(i))*sel(j,io(i))*(1.0-lg(j))/pp;
     }
   }

//proprotion for pre-fishery survey
   pp = 0.0;
   for (int j=1;j<=na;j++)
   {
     pp = pp + (nps(j,ipre)+ops(j,ipre))*sel(j,ipre);
   }
   if (pp <= 0.0) pp = 0.000001;
   for (j=1;j<=na;j++)
   {
     enf(j) = nps(j,ipre)*sel(j,ipre)/pp;
     if (j > 1) eof(j) = ops(j,ipre)*sel(j,ipre)/pp;
   }

//fishing effort for summer fishery
   for (i=1;i<=ny;i++)
   {
     tb(i) = tb(i) - 0.5*tc(i);  //exploitable abundance at the middle of the season
     if (tb(i) < 0.001) tb(i) = 0.001;
     if (i <= 17)
     {
       ete(i) = tc(i)/(q1*tb(i));     //before 1993
     }
     else
     {
       ete(i) = tc(i)/(q2*tb(i));     //after 1992
     }
   }


FUNCTION get_reference_points
// cout<<"start reference points"<<endl;

 ref_Htw.initialize();
 ref_ys = ys(ny-1);
 ref_Hwc = 0.0;
 for (int i = ny-10; i<= ny-1; i++)
 {
   tt0 = 0.0;
   for (int j=3;j<=na;j++) tt0 = tt0 + (npw(j,i)+opw(j,i)*sow)*selw(j)*lg(j);
   if (tt0 <= 0.0) tt0 = 0.000001;
   ref_Htw(i) = (twc(i)+tsc(i))/tt0;
   ref_Hwc = ref_Hwc + ref_Htw(i);
 }
 ref_Hwc = ref_Hwc/10.0;        //use average harvest rate of the most recent 10 years as winter and subsistence catch
 ref_Hts = 0.0;
 for (int j = 1; j<= 101; j++)
 {
   ref_F = 0.01*j-0.01;
   if (j > 1) ref_Hts = 1.0-mfexp(-1.0*ref_F);
   ref_catch.initialize();
   ref_nps.initialize();
   ref_ops.initialize();
// cout<<ref_catch<<endl;
//initial year
   for (int k=1; k<=na; k++)
   {
     ref_nps(k,1) = nps(k,ny);
     ref_ops(k,1) = ops(k,ny);
     ref_ps(k,1)  = ref_nps(k,1)  + ref_ops(k,1);
   }

//numbers at length
   for (i=1;i<100;i++)
   {
     for (int k=1;k<=na;k++)
     {
        ref_catch_s(k,i) = (ref_nps(k,i)+ref_ops(k,i)*sos(k))*sel(k,ny-1)*lg(k)*ref_Hts;
        pp = 0;
        for (int kk=1;kk<=k;kk++)
        {
          tt0 = ref_ps(kk,i)*mfexp(-ref_ys*Mn(kk))-ref_catch_s(kk,i)*mfexp(-(0.583-ref_ys)*Mn(kk));
          if (tt0<0.0) tt0 = 0.0;
          pp = pp + gr(kk,k)*tt0*mp0(kk);
        }
        ref_npw(k,i) = pp;
        if (k==1) ref_npw(k,i) = ref_npw(k,i) + r1*1000.0;
        if (k==2) ref_npw(k,i) = ref_npw(k,i) + (1.0-r1)*1000.0;
        ref_opw(k,i) = (ref_ps(k,i)*mfexp(-ref_ys*Mn(k))-ref_catch_s(k,i)*mfexp(-(0.583-ref_ys)*Mn(k)))*(1.0-mp0(k));
        if (ref_opw(k,i)<0.0) ref_opw(k,i) = 0.0001;
     }
//update to abundance in next summer
     for (k=1;k<=na;k++)
     {
       ref_catch_w_n(k,i) = ref_npw(k,i)*selw(k)*lg(k)*ref_Hwc;
       ref_catch_w_o(k,i) = ref_opw(k,i)*sow*selw(k)*lg(k)*ref_Hwc;
       if (j == 1)                            //all fishing is closed.
       {
          ref_catch_w_n(k,i) = 0.0;
          ref_catch_w_o(k,i) = 0.0;
       }
       ref_catch(k,i) = ref_catch_s(k,i) + ref_catch_w_n(k,i) + ref_catch_w_o(k,i);
       ref_nps(k,i+1) = (ref_npw(k,i)-ref_catch_w_n(k,i))*mfexp(-0.417*Mn(k));
       ref_ops(k,i+1) = (ref_opw(k,i)-ref_catch_w_o(k,i))*mfexp(-0.417*Mo(k));
       ref_ps(k,i+1)  = ref_nps(k,i+1)  + ref_ops(k,i+1);
// cout<<ref_nps(k,i+1)<<endl;
     }
     ref_legal(i) = 0.0;
     ref_legal_b(i) = 0.0;
     ref_mmb(i) = 0.0;
     ref_catch_b(i) = 0.0;
     ref_catch_n(i) = 0.0;
     for (k=1; k<=na; k++)
     {
        ref_legal(i) = ref_legal(i) + ref_ps(k,i)*lg(k);
        ref_legal_b(i) = ref_legal_b(i) + ref_ps(k,i)*lg(k)*wm(k);
        ref_catch_b(i) = ref_catch_b(i) + ref_catch(k,i)*wm(k);
        ref_catch_n(i) = ref_catch_n(i) + ref_catch(k,i);
        if (k > 2)
        {
           ref_mmb(i) = ref_mmb(i) + ref_ps(k,i)*wm(k);
        }
     }
   }
   ref_mbio(j) = ref_mmb(99)/1000.0;                    //lbs/R
   ref_totc(j) = ref_catch_b(99)/1000.0;
   ref_legaln(j) = ref_legal(99)/1000.0;
   ref_legalb(j) = ref_legal_b(99)/1000.0;
   ref_Hr(j) = (ref_catch_n(99)/1000.0)/ref_legaln(j);
 }
// cout<<ref_mbio<<endl;
// cout<<ref_totc<<endl;
// cout<<ref_legaln<<endl;
   i35 = 0;
   i40 = 0;
  for (j = 1; j<= 101; j++)
 {
    if (i35 < 1.0)
    {
      if (ref_mbio(j) <= 0.35*ref_mbio(1))
      {
         f35 = 0.01*j-0.01;
         h35 = ref_Hr(j);
         i35 = 2.0;
      }
    }
    if (i40 < 1.0)
    {
      if (ref_mbio(j) <= 0.40*ref_mbio(1))
      {
         f40 = 0.01*j-0.01;
         h40 = ref_Hr(j);
         i40 = 2.0;
      }
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
  report1 <<"Total legal biomass as F = 0.00, 0.01, ... 1.0 (lbs/R)"<<endl;
  report1 << ref_legalb<<endl;
  report1 <<"F35: "<<f35<<"  H35%: "<<h35<<endl;
  report1 <<"F40: "<<f40<<"  H40%: "<<h40<<endl;
  report1 <<"ref_Hwc = (mean of 10 years):  "<<ref_Hwc<<endl;
// cout<<"end of reference points"<<endl;


FUNCTION evaluate_the_objective_function
  // penalty functions to ``regularize '' the solution
// f = 0.01*norm2(log_relrec);                           //deviation in recruits.
   tf13 = 0.01*norm2(log_relpop);                         //deviation in recruits.
   f = tf13;
//   f += lamt*norm2(log(tt+1.e-3)-log(ett+1.e-3));      //total trawl survey
//   tf1 = norm2(log(tt1+1.e-3)-log(ett1+1.e-3));      //total trawl survey
//   f += lamt*tf1;
//weight each year estimate by 1/(2*variance) - use cv of biomass in sqrt(log(cv^2+1)) as sd of log(biomass)
   tf1 = norm2(elem_div((log(tt1+1.e-3)-log(ett1+1.e-3)),
               sqrt(2)*sqrt(log(elem_prod(cv1,cv1)+1.0))));
   f += tf1;
//   tf2 = norm2(log(tt2+1.e-3)-log(ett2+1.e-3));  //total trawl survey
//   f += 1.0*lamt*tf2;                          //total trawl survey
   tf2 = norm2(elem_div((log(tt2+1.e-3)-log(ett2+1.e-3)),
               sqrt(2)*sqrt(log(elem_prod(cv2,cv2)+1.0))));
   f += tf2;
//   f += lamp*norm2(log(tp+1.e-3)-log(etp+1.e-3));      //total pot survey
//   tf3 = norm2(log(tp1+1.e-3)-log(etp1+1.e-3));      //total pot survey
//   f += lamp*tf3;
   tf3 = norm2(log(tp1+1.e-3)-log(etp1+1.e-3))/(2.0*log(0.34*0.34+1.0));
   f += tf3;
//   tf4 = norm2(log(tp2+1.e-3)-log(etp2+1.e-3));  //total pot survey
//   f += 1.0*lamp*tf4;
   tf4 = norm2(log(tp2+1.e-3)-log(etp2+1.e-3))/(2.0*log(0.27*0.27+1.0));
   f += tf4;
   tf5 = norm2(log(te+1.e-3)-log(ete+1.e-3));        //effort in summer fishery
   f += lamc*tf5;
   f -= 0.5*lamc*pow((log(te(17)+1.e-3)-log(ete(17)+1.e-3)),2.);    //downweight 93 effort
   tf6 = sum(elem_prod(st,colsum(elem_prod(ont,log(ent+1.e-3)))));   //trawl survey
   f -=tf6;
   tf7 = sum(elem_prod(sp,colsum(elem_prod(onp+oop,log(enp+eop+1.e-3)))));    //pot survey
   f -= tf7;
 //  tf8 = sum(elem_prod(sp,colsum(elem_prod(oop,log(eop+1.e-3)))));
 //  f -= tf8;
   tf9 = sum(elem_prod(sw,colsum(elem_prod(onw+oow,log(enw+eow+1.e-3)))));      //winter project
   f -= tf9;
//   tf10 = sum(elem_prod(sw,colsum(elem_prod(oow,log(eow+1.e-3)))));
//   f -= tf10;
   tf14 = pre*(sum(elem_prod(onf+oof,log(enf+eof+1.e-3)))+sum(elem_prod(onf+oof,log(enf+eof+1.e-3))));  //pre-fishery
   f -= tf14;
//   f -= pre*sum(elem_prod(oof,log(eof+1.e-3)));
   tf15 = sum(elem_prod(so,colsum(elem_prod(ono+ooo,log(eno+eoo+1.e-3)))));      //observer's data
   f -= tf15;
//   f -= sum(elem_prod(so,colsum(elem_prod(ooo,log(eoo+1.e-3)))));
   tf11 = sum(elem_prod(sc,colsum(elem_prod(onc+ooc,log(enc+eoc+1.e-3)))));      //summer fishery
   f -= tf11;
//   tf12 = sum(elem_prod(sc,colsum(elem_prod(ooc,log(eoc+1.e-3)))));
//   f -= tf12;
  cout <<"Tr_imma Tr-Mat Pot1  Pot2  effort  Tr_pr  Pot_pr  Win_pr  Sum_pr "<<endl;
  cout <<tf1<<" "<<tf2<<" "<<tf3<<" "<<tf4<<" "<<tf5<<" "<<tf6<<" "<<tf7<<" "<<tf9<<" "<<tf11<<endl;
  cout <<"Recruit deviation Pre-fishery  Observer's data total F"<<endl;
  cout <<tf13<<" "<<tf14<<" "<<tf15<<" "<<f<<endl;
//  cout<<tt1<<endl;
//  cout<<ett1<<endl;
//  cout<<tt2<<endl;
//  cout<<ett2<<endl;
//  cout<<tt<<endl;
//  cout<<ett<<endl;
  cout<<"Legal in ny = "<<legal(ny)/1000.0<<" mmb ny = "<<mmb(ny)/1000.0<<endl;
  cout<<"F35% = "<<f35<<" F40% = "<<f40<<endl;
  cout<<"H35% = "<<h35<<" H40% = "<<h40<<endl;


REPORT_SECTION
  report << "Estimated abundance: new shell, old shell "<< endl;
  report << nps << endl;
  report << ops << endl;
  report << ett << endl;
  report << tt<< endl;
  report << etp << endl;
  report << tp<< endl;
  report << ete << endl;
  report << ent << endl;
  report << eot << endl;
  report << enp << endl;
  report << eop << endl;
  report << enw << endl;
  report << eow << endl;
  report << enf << endl;
  report << eof << endl;
  report << enc << endl;
  report << eoc << endl;
  report << eno << endl;
  report << eoo << endl;
  report << rec << endl;
  report << legal << endl;
  report << mmb << endl;
  report << ett1 << endl;
  report << tt1 << endl;
  report << ett2 << endl;
  report << tt2 << endl;
  report << etp1 << endl;
  report << tp1 << endl;
  report << etp2 << endl;
  report << tp2 << endl;
