#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List FormPairs(int ns){
  
  //ns: n sites
  //np: n pairs
  
  // calculate n possible pairs
  int np = 0;
  np += (pow(ns, 2)-ns)/2;
  
  NumericVector pair1(np);
  NumericVector pair2(np);
  int x = 0;
  
  // loop over sites
  for(int i = 0; i < ns-2; ++i){
    
    for(int j = i+1; j < ns; ++j){
      
      pair1[x] = i+1;
      pair2[x] = j+1;
      
      x+= 1;
      
    }
    
  } 
  
  pair1[x] = ns-1;
  pair2[x] = ns;
  
  List ret;
  ret["pair1"] = pair1;
  ret["pair2"] = pair2;
  
  return ret;
  
}


// [[Rcpp::export]]
List FormPairs2(NumericVector regions, int x){
  
  // regions: must be numeric
  // x: is the number of possible site pairs (e.g. n(unqiue_i) * n(unique_j))
  
  int ll = regions.size();
  NumericVector pair1(x);
  NumericVector pair2(x);
  
  int z = 0;
  
  for(int i = 0; i < ll; ++i){
    
    for(int j = i+1; j < ll; ++j){
    
      if (j > ll) {
        break;
      }
      
      if(regions[i] != regions[j]){
        
        pair1[z] = i+1;
        pair2[z] = j+1;
        
        z+= 1;
        
      }
    
    }
    
  }
  
  List ret;
  ret["pair1"] = pair1;
  ret["pair2"] = pair2;
  
  return ret;

}