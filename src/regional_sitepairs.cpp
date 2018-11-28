#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SelectMismatches(StringVector x,
              StringVector y,
              NumericVector tally,
              NumericVector target) {
  
  // regions (rows) to iterate over
  int n = x.size();
  
  // selected row ids output container
  NumericVector selected(n);
  
  // work out how to make tally an optional arg - for now, pass it in
  // tally regional selection container;
  // NumericVector tally(target.size());
  // tally.names() = target.names();
  
  // which regions are up for consideration
  std::string reg1;
  std::string reg2;
  
  // temporary objects to handle target and tally checks
  // must be a better way of doing this?
  int ta_r1 = 0;
  int ta_r2 = 0;
  int tg_r1 = 0;
  int tg_r2 = 0;
  
  // loop over region pairs
  for(int i = 0; i < n; ++i) {
    
    // paired regions
    reg1 = x[i];
    reg2 = y[i];
    
    //Rcout << "Reg 1: " << reg1 << std::endl;
    //Rcout << "Reg 2: " << reg2 << std::endl;
    
    // regional tallies
    ta_r1 += tally[reg1];
    ta_r2 += tally[reg2];
    
    //Rcout << "Reg tally 1 =" << ta_r1 << std::endl;
    
    // regional targets
    tg_r1 += target[reg1];
    tg_r2 += target[reg2];
    
    if (ta_r1 < tg_r1 && ta_r2 < tg_r2) {
        
      // tallies not met - accumulate  
      tally[reg1] = tally[reg1] + 1;
      tally[reg2] = tally[reg2] + 1;
      
      // and add pair to collection
      selected[i] = 1;
      
    } 
    
    // reset because I don't know any better
    ta_r1 = 0;
    ta_r2 = 0;
    tg_r1 = 0;
    tg_r2 = 0;
      
  }
  
  List ret;
  ret["tally"] = tally;
  ret["selected"] = selected;
  
  return ret;
  
}
  
