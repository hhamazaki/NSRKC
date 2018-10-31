DATA_SECTION
// read the data file
  init_int na                 //number of length groups
  init_number slm             //mid-size of smallest length group (mm)
  init_number slt             //length interval (mm)
  init_vector mp1(1,na)
  init_vector sel1(1,na)
  init_vector sel2(1,na)  
  init_imatrix tag_recov1(1,21,1,8)
  init_imatrix tag_recov2(1,21,1,8)  
  
!! cout<<"End of reading data file success and the end of file code is:"<<eof<<endl;

PARAMETER_SECTION
  init_bounded_number sigma(0.,30.,1)                    
  init_bounded_vector ig(1,na-1,0.,30.,1)       //  growth increment  
  matrix tagrecap1(1,na,1,na)   //tagging freq for Year 1 
  matrix tagrecap2(1,na,1,na)   //tagging freq for Year 2
  matrix tagrecap3(1,na,1,na)   //tagging freq for Year 3
  matrix tagrecap12(1,na,1,na)   //tagging freq for Year 4
  matrix tagrecap22(1,na,1,na)   //tagging freq for Year 5
  matrix tagrecap32(1,na,1,na)   //tagging freq for Year 6
  matrix ptagrecap1(1,na,1,na)   //prob tagging for Year 1 
  matrix ptagrecap2(1,na,1,na)   //prob tagging for Year 2
  matrix ptagrecap3(1,na,1,na)   //prob tagging for Year 3
  matrix ptagrecap12(1,na,1,na)   //prob tagging for Year 4
  matrix ptagrecap22(1,na,1,na)   //prob tagging for Year 5
  matrix ptagrecap32(1,na,1,na)   //prob tagging for Year 6
  matrix mgr1(1,na,1,na)          //molting probability adjstete matrix Year 1
  matrix mgr2(1,na,1,na)          //molting probability adjstete matrix Year 2
  matrix mgr3(1,na,1,na)          //molting probability adjstete matrix Year 3
  
  matrix tgr(1,na,1,na)  

  matrix egr1(1,na,1,na)    // estimated growth matrix for sel1 Year 1
  matrix egr2(1,na,1,na)   //estimated growth matrix for sel1 Year 2
  matrix egr3(1,na,1,na)   //estimated growth matrix for sel1 Year 3
  matrix egr12(1,na,1,na)   //estimated growth matrix for sel2 Year 1
  matrix egr22(1,na,1,na)   //estimated growth matrix for sel2 Year 2
  matrix egr32(1,na,1,na)   //estimated growth matrix for sel2 Year 3 
  vector ef(1,6)
//
   objective_function_value f
//
INITIALIZATION_SECTION
  sigma  5.0
  ig 10.0
  
PRELIMINARY_CALCS_SECTION
// 
  int i; 
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
  cout<<tagrecap1<< endl;
  cout<<tagrecap2<< endl;
  cout<<tagrecap3<< endl;
  cout<<tagrecap12<< endl;
  cout<<tagrecap22<< endl;
  cout<<tagrecap32<< endl;
  
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
  ptagrecap1(1) = 0.0;
  ptagrecap22(6) = 0.0;
  ptagrecap32(5) = 0.0;  
  ptagrecap32(6) = 0.0;    
//  mp1 = 1.0;
  
//
PROCEDURE_SECTION
// create growth matrix
  growth_matrix();
  evaluate_the_objective_function();
  
FUNCTION growth_matrix
//estimate  D growth matrix
// use simpson's rule ofor nomal distribution approximation 
// I(a,b)f(x)dx = [(b-a)/6]*[f(a)+4*f((a+b)/2)+f(b)]  where 
// f(a) = (1/sqrt(2*pi*square(sigma)))*exp(-0.5*square((a-(Lm+IG))/sigma))
// In this case 
// a = slm+slt*(j-1) - 0.5*slt (Lower bound of length class j)
// b = slm+slt*(j-1) + 0.5*slt (Upper bound of length class j)
// (a+b)/2 = slm+slt*(j-1)  (mid length of length class j)
// Lm = slm+slt*(i-1)  mid point of length class i)
// Thus, a-(Lm+IG) =  slm+slt*(j-1) - 0.5*slt - slm - slt*(i-1) - IG
// = slt*[(j-1) - 0.5 - (i-1)] - IG = slt*(j-i-0.5) - IG
   int i, j ; 
   dvariable t1, fa, fb, fab;
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
// Assume that crab does not shrink    
    for (j=i;j<=na;j++)
     {
     t1 = 1/sqrt(2*PI*square(sigma));
     fa = exp(-0.5*square((slt*(j-i-0.5)-ig(i))/sigma));
     fb = exp(-0.5*square((slt*(j-i+0.5)-ig(i))/sigma));
     fab = exp(-0.5*square((slt*(j-i)-ig(i))/sigma));
     tgr(i,j) = t1*slt*(fa+4*fab+fb)/6;
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////// 
//    Adjustment for selectivity
////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  
  
// multiply by fishery selectivity sl1
 
   for (i=1;i<=na;i++)
    {
    egr1(i) = elem_prod(mgr1(i),sel1);
    egr2(i) = elem_prod(mgr2(i),sel1);
    egr3(i) = elem_prod(mgr3(i),sel1);
    egr12(i) = elem_prod(mgr1(i),sel2);
    egr22(i) = elem_prod(mgr2(i),sel2);
    egr32(i) = elem_prod(mgr3(i),sel2);    
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
//  cout << egr1 << endl;
//  cout << egr12 << endl;

  
  
FUNCTION evaluate_the_objective_function

  ef(1) = sum(elem_prod(rowsum(tagrecap1),rowsum(elem_prod(ptagrecap1,log(egr1+1.e-3))))); 
  ef(2) = sum(elem_prod(rowsum(tagrecap2),rowsum(elem_prod(ptagrecap2,log(egr2+1.e-3))))); 
  ef(3) = sum(elem_prod(rowsum(tagrecap3),rowsum(elem_prod(ptagrecap3,log(egr3+1.e-3))))); 
  ef(4) = sum(elem_prod(rowsum(tagrecap12),rowsum(elem_prod(ptagrecap12,log(egr12+1.e-3))))); 
  ef(5) = sum(elem_prod(rowsum(tagrecap22),rowsum(elem_prod(ptagrecap22,log(egr22+1.e-3))))); 
  ef(6) = sum(elem_prod(rowsum(tagrecap32),rowsum(elem_prod(ptagrecap32,log(egr32+1.e-3))))); 
  
  f -= sum(ef);

  
REPORT_SECTION
//
  report<<"growth matrix1"<< endl <<tgr<<endl;
  report<<"growth matrix2"<< endl <<mgr1<<endl; 
  report<<"egr"<< endl << egr1<<endl;
  report<<"ptagrecap1"<< endl << ptagrecap1<<endl;  
  report<<"tagrecap1"<< endl << tagrecap1<<endl;    
  report << "ef" << endl << ef <<endl;
//

//
GLOBALS_SECTION
//
RUNTIME_SECTION
//
TOP_OF_MAIN_SECTION
//

//  exit(44);
 //
 

  