#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "sbp.h"

SBP::SBP (arma::mat M0, std::map<int,arma::uvec> node0){
  if(M0.n_cols != M0.n_rows) Rcpp::Rcout << "SBP:: ncols different thant" << std::endl;
  if(M0.n_cols != node0.size()) Rcpp::Rcout << "SBP:: M dimension and elements in node must match" << std::endl;
  if(M0.n_cols < 2) Rcpp::Rcout << "SBP:: matrix dimension lower than two" << std::endl;

  M = M0;
  node = node0;

  L = arma::uvec(M.n_cols);
  L_length = 0;
  R = arma::uvec(M.n_cols);
  R_length = 0;
  sL = 0;
  sR = 0;
  sC = 0;
  variance = 0;

  initialized = false;
}

arma::vec get_two(int n){
  int first = (int)floor(n * arma::randu(1)[0]);
  int delta = 1 + (int)floor((n-1) * arma::randu(1)[0]);
  int second = (first + delta) % n;
  arma::vec res(2);
  res(0) = first;
  res(1) = second;
  return(res);
}

arma::vec arma_sampling_with_replacement(int n, int k){
  return(arma::floor(n*arma::vec(k).randu()));
}
// [[Rcpp::export]]
arma::uvec arma_sampling_without_replacement(int n, int k){
  arma::vec rpick = arma::vec(k).randu();
  arma::uvec random_vec = arma::uvec(k);
  arma::uvec index =  arma::uvec(n);
  index[0] = 0;
  for(int i = 1; i < n; i++) index[i] = index[i-1] + 1;
  for(int i = 0; i < k; i++){
    int N = (n-i);
    int ipick = floor(N * rpick[i]);
    random_vec[i] = index[ipick];
    index[ipick] = index[N-1];
  }
  return(random_vec);
}
void SBP::init(){
  int n = get_n();
  arma::vec v_initial = get_two(n);
  arma::vec v_random = arma_sampling_with_replacement(3, n);
  //Rcpp::IntegerVector v_initial = Rcpp::sample(n, 2, false);
  //arma::uvec v_random = Rcpp::as<arma::uvec>(Rcpp::sample(3, n, true));
  v_random[v_initial[0]] = 1;
  v_random[v_initial[1]] = 2;
  arma::uvec L0 = find(v_random == 1);
  arma::uvec R0 = find(v_random == 2);
  init(L0,R0);
}
void SBP::init(arma::uvec L0, arma::uvec R0){
  L_length = L0.n_elem;
  L.head(L0.n_elem) = L0;
  R_length = R0.n_elem;
  R.head(R0.n_elem) = R0;

  //Rcpp::Rcout << M;
  //Rcpp::Rcout << L0 << R0;
  double nL = (double)get_nL();
  double nR = (double)get_nR();

  sL = (nR/nL) * arma::accu(M(L0,L0));
  sR = (nL/nR) * arma::accu(M(R0,R0));
  sC = - 2*arma::accu(M(R0,L0));

  variance = (sL + sR + sC) / (nL+nR);
  //Rcpp::Rcout << "Variance: "<< variance << std::endl;
  initialized = true;
}
void SBP::best_improve(){
  int n = get_n();

  arma::uvec uL = getL();
  arma::uvec uR = getR();

  arma::uvec O = arma::zeros<arma::uvec>(n);
  O(uL).fill(1);
  O(uR).fill(1);

  arma::uvec uO = find(O == 0);

  double maxVar = variance;
  bool add = false;
  bool left = false;
  int which = -1;

  for(unsigned int i=0; i < uO.n_elem; i++){
    double newVar = v_addL(uO[i]);
    if(newVar > maxVar){
      maxVar = newVar;
      add = true;
      left = true;
      which = uO[i];
    }
    newVar = v_addR(uO[i]);
    if(newVar > maxVar){
      maxVar = newVar;
      add = true;
      left = false;
      which = uO[i];
    }
  }
  if(uL.n_elem > 1){
    for(unsigned int i=0; i<uL.n_elem; i++){
      double newVar = v_removeL(uL[i]);
      if(newVar > maxVar){
        maxVar = newVar;
        add = false;
        left = true;
        which = uL[i];
      }
    }
  }
  if(uR.n_elem > 1){
    for(unsigned int i=0; i<uR.n_elem; i++){
      double newVar = v_removeR(uR[i]);
      if(newVar > maxVar){
        maxVar = newVar;
        add = false;
        left = false;
        which = uR[i];
      }
    }
  }
  if(which > -1){
    if(add){
      if(left){
        addL(which);
      }else{
        addR(which);
      }
    }else{ // remove
      if(left){
        removeL(which);
      }else{
        removeR(which);
      }
    }
  }
}
void SBP::k_best_improve(int k){

  int n = get_n();

  if(k > n) k = n;

  arma::uvec uL = getL();
  arma::uvec uR = getR();

  arma::uvec O = arma::zeros<arma::uvec>(n);
  O(uL).fill(1);
  O(uR).fill(2);

  arma::uvec IN = arma_sampling_without_replacement(n, k);
  arma::uvec uIN = arma::zeros<arma::uvec>(n);
  uIN(IN).fill(1);

  arma::uvec uO = find(O == 0 && uIN == 1);
  arma::uvec uLsub = find(O == 1 && uIN == 1);
  arma::uvec uRsub = find(O == 2 && uIN == 1);

  double maxVar = variance;
  bool add = false;
  bool left = false;
  int which = -1;

  for(unsigned int i=0; i < uO.n_elem; i++){
    double newVar = v_addL(uO[i]);
    if(newVar > maxVar){
      maxVar = newVar;
      add = true;
      left = true;
      which = uO[i];
    }
    newVar = v_addR(uO[i]);
    if(newVar > maxVar){
      maxVar = newVar;
      add = true;
      left = false;
      which = uO[i];
    }
  }
  if(uLsub.n_elem > 1){
    for(unsigned int i=0; i<uLsub.n_elem; i++){
      double newVar = v_removeL(uLsub[i]);
      if(newVar > maxVar){
        maxVar = newVar;
        add = false;
        left = true;
        which = uLsub[i];
      }
    }
  }
  if(uRsub.n_elem > 1){
    for(unsigned int i=0; i<uRsub.n_elem; i++){
      double newVar = v_removeR(uRsub[i]);
      if(newVar > maxVar){
        maxVar = newVar;
        add = false;
        left = false;
        which = uRsub[i];
      }
    }
  }
  if(which > -1){
    if(add){
      if(left){
        addL(which);
      }else{
        addR(which);
      }
    }else{ // remove
      if(left){
        removeL(which);
      }else{
        removeR(which);
      }
    }
  }
}
bool change_state(double prev_variance, double new_variance, double p){
  //Rcpp::Rcout << "Previous: " << prev_variance << ", new: " << new_variance << std::endl;
  double ratio = new_variance/prev_variance;
  if(1 < ratio) return(true);
  if(arma::vec(1).randu()[0] < sqrt(sqrt(ratio))) return(true);
  return(false);
}

void SBP::simulated_annealing(int steps, int optim){
  arma::uvec bestL;
  arma::uvec bestR;
  double bestVar = 0;
  if(node.size() == 2){
    arma::uvec L0(1), R0(1);
    L0(0) = 0;
    R0(0) = 1;
    init(L0,R0);
  }else{
    init();
    for(int i=0; i<steps; i++){
      int n = get_n();

      arma::uvec uL = getL();
      arma::uvec uR = getR();

      arma::uvec O = arma::zeros<arma::uvec>(n);
      O(uL).fill(1);
      O(uR).fill(2);
      if(uL.n_elem == 1){
        O(uL).fill(4);
      }
      if(uR.n_elem == 1){
        O(uR).fill(4);
      }
      double randu = 0;
      unsigned int iact = 0;
      do{
       randu = arma::vec(1).randu()[0];
       iact = (unsigned int)floor(n * randu);
      }while(O[iact] == 4);
      double new_variance = 0;
      if(O[iact] == 0){
        if(randu < 0.5){
          new_variance = v_addL(iact);
          if( change_state(variance, new_variance, (double)i/double(steps)) ) addL(iact);
        }else{
          new_variance = v_addR(iact);
          if( change_state(variance, new_variance, (double)i/double(steps)) ) addR(iact);
        }
      }else{
        if(O[iact] == 1){
          new_variance = v_removeL(iact);
          if( change_state(variance, new_variance, (double)i/double(steps)) ) removeL(iact);
        }else{
          new_variance = v_removeR(iact);
          if( change_state(variance, new_variance, (double)i/double(steps)) ) removeR(iact);
        }
      }
      double prev_variance = 0;
      for(int j=0; j<optim && prev_variance != variance; j++){
        prev_variance = variance;
        best_improve();
      }
      if(bestVar < variance){
        bestVar = variance;
        bestL = getL();
        bestR = getR();
      }
    }
    init(bestL, bestR);
  }
}

void SBP::simulated_annealing2(int steps, int random, int optim, int k){
  arma::uvec bestL;
  arma::uvec bestR;
  double bestVar = 0;
  if(node.size() == 2){
    arma::uvec L0(1), R0(1);
    L0(0) = 0;
    R0(0) = 1;
    init(L0,R0);
  }else{
    init();
    for(int i=0; i<steps; i++){
      if(random == 0){
        // Random start
        init();
      }else{
        // Random walk
        for(int j=0;j<random;j++){
          int n = get_n();

          arma::uvec uL = getL();
          arma::uvec uR = getR();

          arma::uvec O = arma::zeros<arma::uvec>(n);
          O(uL).fill(1);
          O(uR).fill(2);
          if(uL.n_elem == 1){
            O(uL).fill(4);
          }
          if(uR.n_elem == 1){
            O(uR).fill(4);
          }
          double randu = 0;
          unsigned int iact = 0;
          do{
            randu = arma::vec(1).randu()[0];
            iact = (unsigned int)floor(n * randu);
          }while(O[iact] == 4);
          double new_variance = 0;
          if(O[iact] == 0){
            if(randu < 0.5){
              new_variance = v_addL(iact);
              if( change_state(variance, new_variance, (double)i/double(steps)) ) addL(iact);
            }else{
              new_variance = v_addR(iact);
              if( change_state(variance, new_variance, (double)i/double(steps)) ) addR(iact);
            }
          }else{
            if(O[iact] == 1){
              new_variance = v_removeL(iact);
              if( change_state(variance, new_variance, (double)i/double(steps)) ) removeL(iact);
            }else{
              new_variance = v_removeR(iact);
              if( change_state(variance, new_variance, (double)i/double(steps)) ) removeR(iact);
            }
          }
          if(bestVar < variance){
            bestVar = variance;
            bestL = getL();
            bestR = getR();
          }
        }
      }
      // Optimizing walk
      double prev_variance = 0;
      for(int j=0; j<optim && prev_variance != variance; j++){
        prev_variance = variance;
        k_best_improve(k);
      }
      if(bestVar < variance){
        bestVar = variance;
        bestL = getL();
        bestR = getR();
      }
    }
    init(bestL, bestR);
  }
}


void SBP::first_component_approximation(){
  if(node.size() == 2){
    arma::uvec L0(1), R0(1);
    L0(0) = 0;
    R0(0) = 1;
    init(L0,R0);
  }else{
    arma::uvec bestL;
    arma::uvec bestR;
    double bestVar = 0;

    int K = M.n_cols;
    arma::mat G = arma::eye(K,K) - arma::ones<arma::mat>(K,K)/K;
    arma::mat S = G * M * G;
    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, S);
    arma::vec pc1 = eigvec.col(K-1);
    arma::uvec indices = sort_index(abs(pc1), "descent");

    arma::uvec L0(K), R0(K);
    int nL0 = 0, nR0 = 0;
    for(unsigned int i = 0; i < indices.n_elem; i++){
      if(pc1[indices[i]] < 0){ // Left
        L0[nL0] = indices[i];
        nL0++;
      }else{
        R0[nR0] = indices[i];
        nR0++;
      }
      if(nR0 > 0 && nL0 > 0){
        init(L0.head(nL0), R0.head(nR0));
      }
      if(bestVar < variance){
        bestVar = variance;
        bestL = getL();
        bestR = getR();
      }
    }
    init(bestL, bestR);
  }
}

void SBP::first_component_approximation2(arma::mat Mclr){

  if(node.size() == 2){
    arma::uvec L0(1), R0(1);
    L0(0) = 0;
    R0(0) = 1;
    init(L0,R0);
  }else{
    arma::uvec bestL;
    arma::uvec bestR;
    double bestVar = 0;

    int K = M.n_cols;
    int Ki = Mclr.n_cols;

    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, Mclr);
    arma::vec pc1 = eigvec.col(Ki-1);
    //Rcpp::Rcout << pc1 << std::endl;
    arma::vec pc1_restricted = arma::vec(K);
    for(unsigned i = 0; i < pc1_restricted.n_elem; i++){
      pc1_restricted[i] = mean(pc1(node[i]));
    }
    arma::uvec indices = sort_index(abs(pc1_restricted), "descent");
    arma::uvec L0(K), R0(K);
    int nL0 = 0, nR0 = 0;
    for(unsigned int i = 0; i < indices.n_elem; i++){
      if(pc1_restricted[indices[i]] < 0){ // Left
        L0[nL0] = indices[i];
        nL0++;
      }else{
        R0[nR0] = indices[i];
        nR0++;
      }
      if(nR0 > 0 && nL0 > 0){
        init(L0.head(nL0), R0.head(nR0));
      }
      if(bestVar < variance){
        bestVar = variance;
        bestL = getL();
        bestR = getR();
      }
    }
    if(nR0 == 0 || nL0 == 0){
      if(nR0 == 0){
        nL0--;
        R0[0] = indices[indices.n_elem-1];
        nR0 = 1;
      }else{ // nL0 == 0
        nR0--;
        L0[0] = indices[indices.n_elem-1];
        nL0 = 1;
      }
      init(L0.head(nL0), R0.head(nR0));
    }else{
      init(bestL, bestR);
    }
    double prev_variance = -1;
    double new_variance = variance;
    int iter = 0;
    while(prev_variance != new_variance){
      prev_variance = new_variance;
      best_improve();
      new_variance = variance;
      iter++;
      //Rcpp::Rcout << "New variance:" << new_variance << " Old variance: "<< prev_variance << std::endl;
      if (iter % 1000 == 0) Rcpp::checkUserInterrupt();
    }
  }
}

void SBP::first_pc_local_search(arma::vec pc1){
  //Rcpp::Rcout << "Entering first_pc_local_search" << std::endl;
  if(node.size() == 2){
    arma::uvec L0(1), R0(1);
    L0(0) = 0;
    R0(0) = 1;
    init(L0,R0);
  }else{
    arma::uvec bestL;
    arma::uvec bestR;
    //double bestVar = 0;

    int K = M.n_cols;

    //Rcpp::Rcout << pc1 << std::endl;
    arma::vec pc1_restricted = arma::vec(K);
    for(unsigned i = 0; i < pc1_restricted.n_elem; i++){
      pc1_restricted[i] = arma::norm(pc1(node[i]));
    }
    arma::uvec indices = sort_index(abs(pc1_restricted), "descent");
    arma::uvec L0(K), R0(K);
    int nL0 = 0, nR0 = 0;

    for(unsigned int i = 0; i < indices.n_elem; i++){
      // if(nR0 > 0 && nL0 > 0){
      //   init(L0.head(nL0), R0.head(nR0));
      //   if(pc1_restricted[indices[i]] < 0 && v_addL(indices[i]) > variance){ // Left
      //     L0[nL0] = indices[i];
      //     nL0++;
      //   }
      //   if(pc1_restricted[indices[i]] > 0 && v_addR(indices[i]) > variance){ // Right
      //     R0[nR0] = indices[i];
      //     nR0++;
      //   }
      // }else{
        if(pc1_restricted[indices[i]] < 0){ // Left
          L0[nL0] = indices[i];
          nL0++;
        }else{
          R0[nR0] = indices[i];
          nR0++;
        }
      // }
      // if(bestVar < variance){
      //   bestVar = variance;
      //   bestL = getL();
      //   bestR = getR();
      // }
    }
    if(nR0 == 0 || nL0 == 0){
      if(nR0 == 0){
        nL0--;
        R0[0] = indices[indices.n_elem-1];
        nR0 = 1;
      }else{ // nL0 == 0
        nR0--;
        L0[0] = indices[indices.n_elem-1];
        nL0 = 1;
      }

    }//else{
     // init(bestL, bestR);
    //}
    init(L0.head(nL0), R0.head(nR0));
    double prev_variance = -1;
    double new_variance = variance;
    int iter = 0;
    while(prev_variance != new_variance){
      prev_variance = new_variance;
      best_improve();
      new_variance = variance;
      iter++;
      //Rcpp::Rcout << "New variance:" << new_variance << " Old variance: "<< prev_variance << std::endl;
      if (iter % 1000 == 0) Rcpp::checkUserInterrupt();
    }
  }
  //Rcpp::Rcout << "Exiting first_pc_local_search" << std::endl;
}

void SBP::local_search(int rep){
  arma::uvec bestL;
  arma::uvec bestR;
  double bestVar = 0;
  //
  //
  if(node.size() == 2){
    arma::uvec L0(1), R0(1);
    L0(0) = 0;
    R0(0) = 1;
    init(L0,R0);
  }else{
    for(int i=0; i<rep; i++){
      init();

      double prev_variance = -1;
      double new_variance = variance;
      int iter = 0;
      //sbp->print_status();
      while(prev_variance != new_variance){
        prev_variance = new_variance;
        best_improve();
        new_variance = variance;
        iter++;
        //Rcpp::Rcout << "New variance:" << new_variance << " Old variance: "<< prev_variance << std::endl;
        if (iter % 1000 == 0) Rcpp::checkUserInterrupt();
      }
      if(bestVar < variance){
        bestVar = variance;
        bestL = getL();
        bestR = getR();
      }
    }
    init(bestL, bestR);
  }
}

bool SBP::hasNext(){
  int n = get_n();
  arma::uvec uL = arma::uvec(getL());
  arma::uvec uR = arma::uvec(getR());
  if( (uL.n_elem == 1) & (uR.n_elem + 1 == n) & (uL[0] == n-1) ){
    return(false);
  }
  return(true);
}

void SBP::nextSBP(){
  int n = get_n();
  arma::uvec uL = arma::uvec(getL());
  arma::uvec uR = arma::uvec(getR());
  arma::uvec O = arma::zeros<arma::uvec>(n);
  O(uL).fill(1);
  O(uR).fill(2);
  do{
    int pos = 0;
    while(O[pos] == 2){
      O[pos] = 0;
      pos++;
    }
    O[pos]++;
    uL = find(O == 1);
    uR = find(O == 2);
  } while ( (uL.n_elem == 0) | (uR.n_elem == 0) );
  init(uL,uR);
}

void SBP::optimal(){
  arma::uvec bestL;
  arma::uvec bestR;
  double bestVar = 0;
  //
  //
  //Rcpp::Rcout << "Before optimisation!!" << std::endl;
  //print_status(true,true,false);
  arma::uvec L0(1), R0(1);
  L0(0) = 0;
  R0(0) = 1;
  init(L0,R0);
  bestVar = variance;
  bestL = getL();
  bestR = getR();
  if(node.size() > 2){

    double prev_variance = -1;
    double new_variance = variance;
    int iter = 0;
    //sbp->print_status();
    while(hasNext()){

      prev_variance = new_variance;
      nextSBP();
      new_variance = variance;
      iter++;

      //Rcpp::Rcout << "New variance:" << new_variance << " Old variance: "<< prev_variance << std::endl;
      if (iter % 10000 == 0) Rcpp::checkUserInterrupt();
      if(bestVar < variance){
        bestVar = variance;
        bestL = getL();
        bestR = getR();
      }
    }
    init(bestL, bestR);
  }

  //Rcpp::Rcout << "After optimisation!!" << std::endl;
  //print_status(true,true,false);

}


SBP SBP::top(){
  int n = M.n_cols;

  arma::uvec uL = arma::uvec(getL());
  arma::uvec uR = arma::uvec(getR());
  arma::uvec O = arma::zeros<arma::uvec>(n);
  O(uL).fill(1);
  O(uR).fill(1);

  arma::uvec uV = find(O == 1);
  arma::uvec uI = find(O == 0);

  int nL = get_nL();
  int nR = get_nR();

  int nV = nL+nR;
  int nI = n-(uL.n_elem+uR.n_elem);
  // Rcpp::Rcout << " [";
  // Rcpp::Rcout << "n:" << n << ",  nV:" << nV << ", nL:" << nL << ", nR:" << nR;
  // Rcpp::Rcout << ", nI:" << nI;
  // Rcpp:: Rcout << ", uL.n_elem:" << uL.n_elem << ", uR.n_elem:" << uR.n_elem << "]";
  arma::mat new_M = arma::zeros<arma::mat>(1+nI, 1+nI);
  // Diagonal
  new_M(0,0, arma::size(nI, nI)) = M(uI,uI);
  new_M(nI, nI) = arma::accu(M(uV,uV));
  // Off-diagonal
  arma::mat pSum = sum(M(uV,uI),0);
  new_M(nI,0, arma::size(1, nI)) = pSum;
  new_M(0,nI, arma::size(nI, 1)) = pSum.t();

  // New nodes pointers
  std::map<int,arma::uvec> node0;
  for(int i = 0; i<nI; i++){
    node0[i] = node[uI[i]];
  }
  arma::uvec V(nV);
  int k = 0;
  for(unsigned int i=0; i < L_length; i++)
    for(unsigned int j=0; j < node[L[i]].n_elem; j++, k++)
      V(k) = node[L[i]][j];
  for(unsigned int i=0; i < R_length; i++)
    for(unsigned int j=0; j < node[R[i]].n_elem; j++, k++)
      V(k) = node[R[i]][j];

  node0[nI] = V;
  return(SBP(new_M, node0));
}
SBP SBP::left(){
  arma::uvec uL = getL();
  arma::uvec uR = getL();
  arma::mat new_M = M(uL, uR);
  std::map<int,arma::uvec> node0;
  for(unsigned int i = 0; i<uL.n_elem; i++){
    node0[i] = arma::uvec(node[uL[i]]);
  }
  return(SBP(new_M, node0));
}
SBP SBP::right(){
  arma::uvec uL = getR();
  arma::uvec uR = getR();
  arma::mat new_M = M(uL, uR);
  std::map<int,arma::uvec> node0;
  for(unsigned int i = 0; i<uR.n_elem; i++){
    node0[i] = arma::uvec(node[uR[i]]);
  }
  return(SBP(new_M, node0));
}
int SBP::get_nL(){
  double nL = 0;
  for(unsigned int i=0; i<L_length;nL+=node[L[i++]].n_elem);
  return(nL);
}
int SBP::get_nR(){
  double nR = 0;
  for(unsigned int i=0; i<R_length;nR+=node[R[i++]].n_elem);
  return(nR);
}

double SBP::v_addL(int I){
  double nL = (double)get_nL();
  double nR = (double)get_nR();
  double ni = (double)node[I].n_elem;

  arma::uvec uI(1);
  uI[0] = I;

  double sL = this->sL * nL/(nL+ni) + (nR/(nL+ni)) * (2*arma::accu(M(uI,getL())) + M(I,I));
  double sR = this->sR * (nL+ni)/nL;
  double sC = this->sC - 2*arma::accu(M(getR(),uI));

  return( (sL + sR + sC) / (nL+nR+ni) );

}
void SBP::addL(int I){
  double nL = (double)get_nL();
  double nR = (double)get_nR();
  double ni = (double)node[I].n_elem;

  arma::uvec uI(1);
  uI[0] = I;

  sL = sL * nL/(nL+ni) + (nR/(nL+ni)) * (2*arma::accu(M(uI,getL())) + M(I,I));
  sR = sR * (nL+ni)/nL;
  sC = sC - 2*arma::accu(M(getR(),uI));

  variance = (sL + sR + sC) / (nL+nR+ni);

  L[L_length++] = I;
}
double SBP::v_addR(int I){
  double nL = (double)get_nL();
  double nR = (double)get_nR();
  double ni = (double)node[I].n_elem;

  arma::uvec uI(1);
  uI[0] = I;

  double sL = this->sL * (nR+ni)/nR;
  double sR = this->sR * nR/(nR+ni) + (nL/(nR+ni)) * (2*arma::accu(M(uI,getR())) + M(I,I));

  double sC = this->sC - 2*arma::accu(M(getL(),uI));

  return( (sL + sR + sC) / (nL+nR+ni) );

}
void SBP::addR(int I){
  double nL = (double)get_nL();
  double nR = (double)get_nR();
  double ni = (double)node[I].n_elem;

  arma::uvec uI(1);
  uI[0] = I;

  sL = sL * (nR+ni)/nR;
  sR = sR * nR/(nR+ni) + (nL/(nR+ni)) * (2*arma::accu(M(uI,getR())) + M(I,I));

  sC = sC - 2*arma::accu(M(getL(),uI));

  variance = (sL + sR + sC) / (nL+nR+ni);

  R[R_length++] = I;
}
double SBP::v_removeL(int I){
  double ni = (double)node[I].n_elem;
  double nL = (double)get_nL()-ni;
  double nR = (double)get_nR();

  arma::uvec uI(1);
  uI[0] = I;
  arma::uvec uL = find(getL() != I);

  double sL = (this->sL - nR/(nL+ni) * (M(I,I) + 2 * arma::accu(M(uI,L(uL))))) * (nL+ni) / nL;
  double sR = this->sR * nL/(nL+ni);
  double sC = this->sC + 2*arma::accu(M(getR(),uI));

  return( (sL + sR + sC) / (nL+nR) );

}
void SBP::removeL(int I){
  double ni = (double)node[I].n_elem;
  double nL = (double)get_nL()-ni;
  double nR = (double)get_nR();

  arma::uvec uI(1);
  uI[0] = I;
  arma::uvec uL = find(getL() != I);

  sL = (sL - nR/(nL+ni) * (M(I,I) + 2 * arma::accu(M(uI,L(uL))))) * (nL+ni) / nL;
  sR = sR * nL/(nL+ni);
  sC = sC + 2*arma::accu(M(getR(),uI));

  variance = (sL + sR + sC) / (nL+nR);

  L_length--;
  L.head(L_length) = L(uL);

}
double SBP::v_removeR(int I){
  double ni = (double)node[I].n_elem;
  double nL = (double)get_nL();
  double nR = (double)get_nR()-ni;

  arma::uvec uI(1);
  uI[0] = I;
  arma::uvec uR = find(getR() != I);

  double sR = (this->sR - nL/(nR+ni) * (M(I,I) + 2 * arma::accu(M(uI,R(uR))))) * (nR+ni) / nR;
  double sL = this->sL * nR/(nR+ni);
  double sC = this->sC + 2*arma::accu(M(getL(),uI));

  return( (sL + sR + sC) / (nL+nR) );

}
void SBP::removeR(int I){
  double ni = (double)node[I].n_elem;
  double nL = (double)get_nL();
  double nR = (double)get_nR()-ni;

  arma::uvec uI(1);
  uI[0] = I;
  arma::uvec uR = find(getR() != I);

  sR = (sR - nL/(nR+ni) * (M(I,I) + 2 * arma::accu(M(uI,R(uR))))) * (nR+ni) / nR;
  sL = sL * nR/(nR+ni);
  sC = sC + 2*arma::accu(M(getL(),uI));

  variance = (sL + sR + sC) / (nL+nR);

  R_length--;
  R.head(R_length) = R(uR);

}



void SBP::print_nodes(){
  for(unsigned int i=0; i < node.size(); i++){
    Rcpp::Rcout << "Node " << i << ": ";
    for(unsigned int j=0; j < node[i].n_elem; j++){
      Rcpp::Rcout << node[i][j] << " ";
    }
    Rcpp::Rcout << std::endl;
  }
}
void SBP::print_LR(){
  Rcpp::Rcout << "L:";
  for(unsigned int i=0; i < L_length; i++){
    Rcpp::Rcout << " " << L[i];
  }
  Rcpp::Rcout << ", R:";
  for(unsigned int i=0; i < R_length; i++){
    Rcpp::Rcout << " " << R[i];
  }
  Rcpp::Rcout << std::endl;
}
void SBP::print_status(bool show_LR = false, bool show_nodes = false,
                       bool show_M = false){
  if(!initialized) Rcpp::Rcout << "Matrix not initialized!!! ";
  Rcpp::Rcout << "nL:" << get_nL();
  Rcpp::Rcout << " ";
  Rcpp::Rcout << "nR:" << get_nR();
  Rcpp::Rcout << " ";
  Rcpp::Rcout << "sL:" << sL << ", sR:" << sR << ", sC:" << sC << ",var: " << variance;
  Rcpp::Rcout << std::endl;
  if(show_LR){
    print_LR();
  }
  if(show_nodes){
    print_nodes();
  }
  if(show_M){
    print_M();
  }
}
