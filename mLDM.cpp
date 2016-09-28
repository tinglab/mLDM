// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>

using namespace Rcpp;
using namespace Eigen;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]

// Get relative abundance
MatrixXd toRatio(MatrixXd x) {
  int n = x.rows();
  
  VectorXd xRowSums = x.rowwise().sum();
  for(int i = 0; i < n; ++i) {
    x.row(i) /= xRowSums[i];
  }
  
  return x;
}

// Scale matrix with mean=0 and std=1
MatrixXd scale(MatrixXd m) {
  int q = m.cols();
  int n = m.rows();
  
  for(int i = 0; i < q; ++i) {
    ArrayXd perCol = m.col(i);
    perCol -= perCol.mean();
    double perStd = sqrt((perCol * perCol).sum() / (n-1));
    m.col(i) = perCol / perStd;
  }
  
  return m;
}

// compute alphai = exp(u_i + z_i)
VectorXd getAlpha(VectorXd& mui, VectorXd& zi) {
  return (mui + zi).array().exp();
}

// Compute lgamma for vector
VectorXd lgammav(VectorXd a) {
  int length = a.rows();
  for(int i = 0; i < length; ++i)
    a(i) = R::lgammafn(a(i));
  
  return a;
}

MatrixXd lgammam(MatrixXd a) {
  int cols = a.cols();
  for(int i = 0; i < cols; ++i) {
    a.col(i) = lgammav(a.col(i));
  }
  
  return a;
}

// Compute digamma for vector
VectorXd digammav(VectorXd a) {
  int length = a.rows();
  for(int i = 0; i < length; ++i)
    a(i) = R::digamma(a(i));
  
  return a;
}

MatrixXd digammam(MatrixXd a) {
  int cols = a.cols();
  for(int i = 0 ; i < cols; ++i) {
    a.col(i) = digammav(a.col(i));
  }
  
  return a;
}

// compute the lgamms part of objective funciton for zi
double computeAlphaVector(VectorXd& ax, VectorXd& a) {
  double obj = - (lgammav(ax) - lgammav(a)).sum() + (R::lgammafn(ax.sum()) - R::lgammafn(a.sum()));
  return obj;
}

// compute the objective function for the zi
NumericVector objZi(SEXP xv, SEXP env) {
  VectorXd zi(as<Map<VectorXd> >(xv));
  Environment paras = as<Environment>(env);
  VectorXd mui(as<Map<VectorXd> >(paras["mu_i"]));
  VectorXd xi(as<Map<VectorXd> >(paras["x_i"]));
  MatrixXd Theta(as<Map<MatrixXd> >(paras["Theta"]));
  int n = as<int>(paras["n"]);
  VectorXd B0(as<Map<VectorXd> >(paras["B0"]));
  
  VectorXd alphai = getAlpha(mui, zi);
  VectorXd alphaix = alphai + xi;
  double obj = computeAlphaVector(alphaix, alphai);
  obj += ((zi - B0).transpose() * Theta * (zi - B0))(0, 0) / 2;
  
  return NumericVector::create(obj / n);
}

// compute the lgamms part of objective funciton for zi
VectorXd computeAlphaDerVector(VectorXd& ax, VectorXd& a) {
  VectorXd obj = - (digammav(ax) - digammav(a)).array() + (R::digamma(ax.sum()) - R::digamma(a.sum()));
  return obj;
}

// compute the derivate of the objective function for zi
NumericVector derObjZi(SEXP xv, SEXP env) {
  VectorXd zi(as<Map<VectorXd> >(xv));
  Environment paras = as<Environment>(env);
  VectorXd mui(as<Map<VectorXd> >(paras["mu_i"]));
  VectorXd xi(as<Map<VectorXd> >(paras["x_i"]));
  MatrixXd Theta(as<Map<MatrixXd> >(paras["Theta"]));
  int n = as<int>(paras["n"]);
  VectorXd B0(as<Map<VectorXd> >(paras["B0"]));
  
  VectorXd alphai = getAlpha(mui, zi);
  VectorXd alphaix = alphai + xi;
  VectorXd der = computeAlphaDerVector(alphaix, alphai).array() * alphai.array() + (Theta * (zi - B0)).array();
    
  return wrap(der/n);
}

// Compute the edges of Theta
int computeEdgesTheta(MatrixXd& Theta) {
  int edges = 0;
  int p = Theta.rows();
  for(int i = 0; i < p; ++i) {
    for(int j = 0; j < i; ++j) {
      if(Theta(i,j) != 0)
        edges += 1;
    }
  }
  
  return edges;
}

// Compute the edges of B
int computeEdgesB(MatrixXd& B) {
  int edges = 0;
  int q = B.rows();
  int p = B.cols();
  for(int i = 0; i < q; ++i) {
    for(int j = 0; j < p; ++j) {
      if(B(i,j) != 0)
        edges += 1;
    }
  }
  
  return edges;
}

// Compute vary of edges of Theta
int computeEdgesVaryTheta(MatrixXd& Theta, MatrixXd& ThetaOld) {
  int edges = 0;
  int p = Theta.rows();
  bool per1, per2;
  
  for(int i = 0; i < p; ++i) {
    for(int j = 0; j < i; ++j) {
      (Theta(i,j) != 0) ? per1 = true : per1 = false;
      (ThetaOld(i,j) != 0) ? per2 = true : per2 = false;
      if(! (per1&& per2)) {
        edges += 1;
      }
    }
  }
  
  return edges;
}

// Compute vary of edges of B
int computeEdgesVaryB(MatrixXd& B, MatrixXd& BOld) {
  int edges = 0;
  int q = B.rows();
  int p = B.cols();
  bool per1, per2;
  
  for(int i = 0; i < q; ++i) {
    for(int j = 0; j < p; ++j) {
      (B(i,j) != 0) ? per1 = true : per1 = false;
      (BOld(i,j) != 0) ? per2 = true : per2 = false;
      if(! (per1 && per2)) {
        edges += 1;
      }
    }
  }
  
  return edges;
}

double computeAlphaMatrix(MatrixXd& AlphaMX, MatrixXd& AlphaM) {
  double obj = - (lgammam(AlphaMX) - lgammam(AlphaM)).sum() + (lgammav(AlphaMX.colwise().sum()) - lgammav(AlphaM.colwise().sum())).sum();
  return obj;
}

MatrixXd computeAlphaDerMatrix(MatrixXd& AlphaMX, MatrixXd& AlphaM) {
  MatrixXd der = (-(digammam(AlphaMX) - digammam(AlphaM))).rowwise() + (digammav(AlphaMX.colwise().sum()) - digammav(AlphaM.colwise().sum())).transpose();
  return der;
}

double computeLogDet(MatrixXd& x) {
  SelfAdjointEigenSolver<MatrixXd> eigensolver(x);
  return eigensolver.eigenvalues().array().log().sum();
}

// Compute the current objective function for f
double computeObjf(MatrixXd& X, MatrixXd& M, MatrixXd& B, MatrixXd& Theta, MatrixXd& Z, double& lambda1, double& lambda2) {
  int n = X.rows();
  VectorXd B0 = Z.colwise().mean();
  // cout << "B0 rows " << B0.rows() << endl;
  MatrixXd Zcenter = (Z.transpose().colwise() - B0).transpose();
  MatrixXd AlphaM = ((M * B).transpose() + Z.transpose()).array().exp();
  MatrixXd AlphaMX = AlphaM + X.transpose();
  MatrixXd S = Zcenter.transpose() * Zcenter / n;
  /*cout << "max AlphaMX " << AlphaMX.maxCoeff();
  cout << "min AlphaMX " << AlphaMX.minCoeff();
  cout << "max AlphaMX " << AlphaM.maxCoeff();
  cout << "min AlphaMX " << AlphaM.minCoeff();*/
  
  double obj = computeAlphaMatrix(AlphaMX, AlphaM) / n;
/*   cout << "obj in objf " << obj << endl;
   cout << "logdet theta" << computeLogDet(Theta) / 2 << endl;
   cout << "s*theta" << (S * Theta).diagonal().sum() / 2 << endl;*/
  obj = obj - computeLogDet(Theta) / 2 + (S * Theta).diagonal().sum() / 2 + lambda1 * Theta.array().abs().sum() / 2 + lambda2 * B.array().abs().sum();
  return obj;
}

// Compute the EBIC
double computeEBIC(double& objNew, MatrixXd& B, MatrixXd& Theta, int& p, int& q, int& n, double& lambda1, double& lambda2) {
  double g = 0.5;
  int E1 = computeEdgesTheta(Theta);
  int E2 = computeEdgesB(B);
  double EBIC = 2*n*(objNew - lambda2*B.array().abs().sum() - lambda1 * Theta.array().abs().sum() / 2)
              + (E1 + E2) * log(n) + 4 * g * E1 * log(p) + 2 * g * E2 * log(p*q);
  
  return EBIC;               
}

// Compute the objective function for index-th row of the matrix B
double objBp(VectorXd& xv, int& index, MatrixXd& BZ, MatrixXd& X, MatrixXd& M, int& q, int& p, int& n) {
  MatrixXd AlphaM = (BZ + xv * M.col(index).transpose()).array().exp();
  MatrixXd AlphaMX = AlphaM + X.transpose();
  double obj = computeAlphaMatrix(AlphaMX, AlphaM);
  
  // << "obj Bi " << obj / n<< endl;
  
  return obj / n;
}

MatrixXd cBind(MatrixXd a, MatrixXd b) {
  int row = a.rows();
  int col = a.cols() + b.cols();
  MatrixXd c(row, col);
  
  c.block(0,0,row, a.cols()) = a;
  c.block(0,a.cols(),row, b.cols()) = b;
  
  return c;
}

MatrixXd rBind(MatrixXd a, MatrixXd b) {
  int col = a.cols();
  int row = a.rows() + b.rows();
  MatrixXd c(row, col);
  
  c.block(0,0,a.rows(), col) = a;
  c.block(a.rows(), 0,b.rows(), col) = b;
  
  return c;
}

void computeBtQ(int& d, double& gammat, MatrixXd& St, MatrixXd& Yt, int& approx_num, int& approxCount, int& approxIndex, MatrixXd& Q, MatrixXd& Qh) {
  if(approxCount != 0) {
    MatrixXd Stp(d, approxCount);
    MatrixXd Ytp(d, approxCount);

    int basis = 0;
    if(approxCount == approx_num)
      basis = approxIndex + 1;
    for(int i = 0; i < approxCount; ++i) {
      Stp.col(i) = St.col((i + basis) % approx_num);
      Ytp.col(i) = Yt.col((i + basis) % approx_num);
    }
    
    // cout << "Stp " << Stp << endl;
    // cout << "Ytp " << Ytp << endl;
    
    Q = cBind(gammat*Stp, Ytp);
    // cout << "Qin" << Q << endl;
    // cout << "gammat" << gammat << endl;
    MatrixXd SY = Stp.transpose() * Ytp;
    MatrixXd Dt = SY.diagonal().asDiagonal();
    
    MatrixXd Lt = SY.triangularView<StrictlyLower>();
    // cout << "SY" << SY << endl;
    // cout << "Dt" << Dt << endl;
    MatrixXd Rt = rBind(cBind(gammat*Stp.transpose() * Stp, Lt), cBind(Lt.transpose(), - Dt));
    Qh = Rt.inverse() * Q.transpose();
    // cout << "Qhin" << Qh << endl;
  }
  
  return;
}

// Compute the derivate the obj for index-th row of the matrix B
VectorXd derObjBp(VectorXd& xv, int& index, MatrixXd& BZ, MatrixXd& X, MatrixXd& M, int& q, int& p, int& n) {
  MatrixXd AlphaM = (BZ + xv*M.col(index).transpose()).array().exp();
  MatrixXd AlphaMX = AlphaM + X.transpose();
  VectorXd der = (computeAlphaDerMatrix(AlphaMX, AlphaM).array() * AlphaM.array()).matrix() * M.col(index);
  
  // cout  << "der Bi" << der / n << endl;
  
  return der / n;
}

vector<int> getActiveSet(VectorXd& w, VectorXd& gt, double& lambda, int& d) {
  vector<int> index;
  double threshold = 1e-6;
  double a, b, subgt;
  bool addFlag = false;
  
  for(int i = 0; i < d; ++i) {
    a = w[i];
    b = gt[i];
    addFlag = true;
    if(abs(a) < threshold) {
      subgt = max(0.0, max(b - lambda, -(b + lambda)));
      
      if(abs(subgt) < threshold) {
        addFlag = false;
      }
    }
    
    if(addFlag) {
      index.push_back(i);
    }
  }
  
  return index;
}

double softthreshold(double a, double b) {
  return a > 0 ? max(a - b, 0.0) : -1 * max(-a - b, 0.0);
}

// Compute the direction via coordinate descent
VectorXd coordinateDescent(int& d, VectorXd& w, double& gammat, MatrixXd& Q, MatrixXd& Qh, VectorXd& gt, double& lambda,
                          int& max_iteration, double& threshold) {
  VectorXd wt = w;
  VectorXd wtOld;
  vector<int> activeSet = getActiveSet(w, gt, lambda, d);
  int activeSize = activeSet.size();
  
  int round = 0;
  double delta = 1;
  int i;
  double a, b;
  MatrixXd Bt = gammat * MatrixXd::Identity(d,d) - Q * Qh;
  
  while(delta > threshold && round < max_iteration) {
    round += 1;
    wtOld = wt;
    
    for(int index = 0 ; index < activeSize; ++index) {
      i = activeSet[index];
      a = Bt(i,i);
      b = gt[i] + (Bt.row(i) * (wt - w)) - a * wt[i];
      
      wt[i] = softthreshold(-b, lambda) / a;
    }
    
    delta = (wt - wtOld).array().abs().sum() / d;
    
  }
  
  VectorXd direction = wt - w;
  return direction.normalized();
}

void filterDir(VectorXd& w) {
  int d = w.rows();
  for(int i = 0 ; i < d; ++i) {
    if(abs(w[i]) < 1e-10)
      w[i] = 0;
  }
  return;
}

// Find new w via linesearch based on strong wolfe condition
List linesearch(VectorXd& w, VectorXd& direction, double& lambda, int& max_linesearch, double& f0, VectorXd& g0, double& delta1, double& delta2,
                          int& index, MatrixXd& BZ, MatrixXd& X, MatrixXd& M, int& q, int& p, int& n) {
  double beta = 0.5;
  double alpha = 2;
  VectorXd wt = w;
  List ret(4);
  int k = 0;
  
  double f1 = f0;
  double dg1 = 0;
  VectorXd g1 = g0;
  
  if(!direction.allFinite() || !g0.allFinite()) {
    cout << "linesearch direction failed!" << endl;
    ret["exist"] = false;
    ret["wt"] = w;
    ret["value"] = f0;
    ret["grad"] = g0;
    
    return ret;
  }
  
  double d1 = g0.transpose() * direction;
  double dg0 = d1 + lambda * w.array().abs().sum();
  d1 += lambda * ((w + direction).array().abs().sum() - w.array().abs().sum());
  d1 *= delta1;
  
  double part;
  
  while(k <= max_linesearch) {
    alpha *= beta;
    wt = w + alpha * direction;
    f1 = objBp(wt, index, BZ, X, M, q, p, n) + lambda * wt.array().abs().sum();
    part = alpha * d1;
    k += 1;
    
    if(f1 <= f0 + part) {
      g1 = derObjBp(wt , index, BZ, X, M, q, p, n);
      dg1 = g1.transpose() * direction + alpha * lambda * (w + direction).array().abs().sum();
      if(abs(dg1 / dg0) <= delta2) {
        break;
      }
    }
  }
  
  filterDir(wt);
  
  ret["exist"] = true;
  ret["wt"] = wt;
  ret["value"] = f1;
  ret["grad"] = g1;
  return ret;
}

// Estimate the matrix B via proximal qusi newton
List proximalQusiNewton(VectorXd w, double& lambda, int& approx_num, int& max_linesearch, int& max_iteration, double& threshold, double& delta1_threshold, 
                        double& delta2_threshold, double& sy_threshold, int& max_iteration_coor, double& threshold_coor, 
                        MatrixXd& BZ, int&index, MatrixXd& X, MatrixXd& M, int& q, int& p, int& n, bool& verbose) {
  int approxCount = 0;
  int approxIndex = -1;
  double gammat = 1;
  double gammat0 = 1;
  
  int d = w.rows();
  MatrixXd St(d, approx_num);
  St.setZero();
  MatrixXd Yt(d, approx_num);
  Yt.setZero();
  MatrixXd Q(d, 2);
  Q.setZero();
  MatrixXd Q0(d, 2);
  Q0.setZero();
  MatrixXd Qh(2, d);
  Qh.setZero();
  MatrixXd Qh0(2, d);
  Qh0.setZero();
  VectorXd direction(d);
  VectorXd wt = w;
  
  int round = 0;
  // When SYdot < sy_threshold, change to steepest descent
  bool steepestActive = false;
  int steepestActiveHeight = 1;
  int steepestCount = 0;
  
  // Record obj values, the gradient for B
  double objNew = 0;
  double objOld = objBp(w, index, BZ, X, M, q, p, n) + lambda * w.array().abs().sum();
  double delta = 0;
  VectorXd gt0 = derObjBp(w, index, BZ, X, M, q, p, n);
  // cout << "gt0" << gt0 << endl;
  VectorXd gt1(d);
  VectorXd StPer(d);
  VectorXd YtPer(d);
  double SYDot = 0;
  List line(4);
  bool exist;
  
  List ret(3);
  
  while(true) {
    round += 1;
    // cout << "St " << St << endl;
    // cout << "Yt " << Yt << endl;
    // Obtain direction
    if(!steepestActive) {
      computeBtQ(d, gammat, St, Yt, approx_num, approxCount, approxIndex, Q, Qh);
      direction = coordinateDescent(d, w, gammat, Q, Qh, gt0, lambda, max_iteration_coor, threshold_coor);
    } else { // when the hessian nears singular, change to steepest descent
      direction = coordinateDescent(d, w, gammat0, Q0, Qh0, gt0, lambda, max_iteration_coor, threshold_coor);
      steepestCount += 1;
      if(steepestCount >= steepestActiveHeight) {
        steepestActive = false;
      }
    }
    /*cout << "Q " << Q << endl;
    cout << "Qh" << Qh << endl;
    cout << "direction" << direction << endl;*/
    
    line = linesearch(w, direction, lambda, max_linesearch, objOld, gt0, delta1_threshold, delta2_threshold,
                          index, BZ, X, M, q, p, n);
    exist = as<bool>(line["exist"]);
    if(!exist) {
      break;
    }
    
    wt = as<Map<VectorXd> >(line["wt"]);
    gt1 = as<Map<VectorXd> >(line["grad"]);
    objNew = as<double>(line["value"]);
    
    
    // cout << "wt " << wt << endl;
    
    delta = objNew - objOld;
    StPer = wt - w;
    YtPer = gt1 - gt0;
    SYDot = StPer.transpose() * YtPer;
    if(SYDot > sy_threshold) {
      approxCount > approx_num ? approxCount = approx_num: approxCount += 1;
      
      approxIndex = (approxIndex + 1) % approx_num;
      St.col(approxIndex) = StPer;
      Yt.col(approxIndex) = YtPer;
      gammat = SYDot / StPer.squaredNorm();
    } else {
      steepestActive = true;
      steepestCount = 0;
      if(verbose) {
        //cout << "steepest descent active!" << endl;
      }
    }
    
    objOld = objNew;
    gt0 = gt1;
    w = wt;
    
    if(abs(delta) < threshold || round > max_iteration) {
      break;
    }
  }
  ret["par"] = w;
  ret["value"] = objOld;
  ret["gradient"] = gt0;
  
  return ret;
}

MatrixXd derObjB(MatrixXd& B, MatrixXd& ZB0, MatrixXd& X, MatrixXd& M, int& n) {
  MatrixXd AlphaM = ((M * B).transpose() + ZB0).array().exp();
  MatrixXd AlphaMX = AlphaM + X.transpose();
  MatrixXd der = (computeAlphaDerMatrix(AlphaMX, AlphaM).array() * AlphaM.array()).matrix() * M;
  
  return der.transpose() / n;
}

// Optimize the OTU-OTU associations and EF-OTU associations via block coordinate descent
List LognormalDirichletMultinomial(MatrixXd X, MatrixXd M, int& n, int& p, int& q, MatrixXd B, VectorXd B0, MatrixXd Theta, MatrixXd Z, double& lambda1, double& lambda2, 
        int& max_iteration, double& threshold, int& approx_num_Z, int& max_linesearch_Z, int& approx_num_B, 
        int& max_linesearch_B, int& max_iteration_B, double& threshold_B, double& delta1_threshold_B, double& delta2_threshold_B,
        double& sy_threshold_B, int& max_iteration_B_coor, double& threshold_B_coor, bool& verbose, Function& lbfgsc, Function& quicc, bool& loopFlag) {
  // Record values of the objective function between two iterations
  double objOldB = 0.0;
  double objNewB = 0.0;
  double objOldTheta = 0.0;
  double objNewTheta = 0.0;
  double delta = 0.0;
  double deltaB = 0.0;
  double deltaTheta = 0.0;
  double EBIC = 0.0;
  double objNew = 0.0;
  double objOld = 0.0;
  
  // The round of iterations
  int round = 0;
  int roundB = 0;
  int roundTheta = 0;
  
  // Record old values for B and Theta
  MatrixXd BOld = B;
  MatrixXd ThetaOld = Theta;
  // Build environment for Zi
  Environment env = new_env();
  MatrixXd derZ(n, p);
  typedef NumericVector(*funcPtr) (SEXP, SEXP);
  XPtr<funcPtr> obj = XPtr<funcPtr>(new funcPtr(&objZi));
  XPtr<funcPtr> grad = XPtr<funcPtr>(new funcPtr(&derObjZi));
  
  while(true) {
    round += 1;
    if(verbose) {
      //cout << round << endl;
    }
    // Estimate Z
    for(int i = 0; i < n; ++i) {
      VectorXd zi = Z.row(i);
      VectorXd xi = X.row(i);
      VectorXd mui = B.transpose() * M.row(i).transpose();
      
      // Put parameters into env
      env.assign("mu_i", mui);
      env.assign("x_i", xi);
      env.assign("Theta", Theta);
      env.assign("n", n);
      env.assign("B0", B0);
      
      List ziResult = lbfgsc(Named("call_eval", obj), Named("call_grad", grad), Named("vars", zi), Named("environment", env), Named("invisible", 1),
                        Named("m", approx_num_Z), Named("max_linesearch", max_linesearch_Z));
                        
      Z.row(i) = as<Map<VectorXd> >(ziResult["par"]).transpose();
      derZ.row(i) = as<Map<VectorXd> >(derObjZi(wrap(Z.row(i)), env)).transpose();
    }
    
    
    if(verbose) {
      /*cout << "max Z" << Z.maxCoeff() << "min Z " << Z.minCoeff() << endl;
      cout << "Z norm2: " << sqrt(Z.squaredNorm()) << endl;
      cout << "derZ norm2: " << sqrt(derZ.squaredNorm()) << endl;*/
    }
    
    // Estimate B0
    B0 = Z.colwise().mean();
    
    // Estimate B
    if(loopFlag) {
      roundB += 1;
      objOldB = objNewB;
      
      MatrixXd ZB0 = Z.transpose();
      // Optimize the B[i,] respectively
      for(int i = 0 ; i < q; ++i) {
        VectorXd bv = B.row(i).transpose();
        MatrixXd Bp = B;
        Bp.row(i).setZero();
        MatrixXd BZ = (M * Bp).transpose() + ZB0;
        
        // cout << "bv " << bv << endl;
        // cout << "BZ " << BZ << endl;
        
        List biResult = proximalQusiNewton(bv, lambda2, approx_num_B, max_linesearch_B, max_iteration_B, 
                                           threshold_B, delta1_threshold_B, delta2_threshold_B, sy_threshold_B, 
                                           max_iteration_B_coor, threshold_B_coor, BZ, i, X, M, q, p, n, verbose);
        
        B.row(i) = as<Map<VectorXd> >(biResult["par"]).transpose();
      }
      
      objNewB = computeObjf(X, M, B, Theta, Z, lambda1, lambda2);
      // cout << "obj New B " << objNewB << endl;
      deltaB = abs(objNewB - objOldB);
      
      // B is finished
      if(deltaB < threshold || round >= max_iteration) {
        loopFlag = false;
        objNew = objNewB;
        objOld = objNewTheta;
        
        if(verbose) {
          MatrixXd derB = derObjB(B, ZB0, X, M, n);
          cout << "Max B " << B.maxCoeff() << "min B " << B.minCoeff() << endl;
          cout << "Round B is:" << roundB << endl;
          cout << "Norm2 B is:" << sqrt(B.squaredNorm()) << endl;
          cout << "Norm2 derB is:" << sqrt(derB.squaredNorm()) << endl;
        }
        roundB = 0;
      } else {
        continue;
      }
    } else {// Estimate Theta
      roundTheta += 1;
      objOldTheta = objNewTheta;
      
      // Estimate Theta via graphical lasso
      MatrixXd Zcenter = (Z.transpose().colwise() - B0).transpose();
      MatrixXd S = Zcenter.transpose() * Zcenter / n;
      List quicRes = quicc(S, lambda1);
      Theta = as<Map<MatrixXd> >(quicRes["X"]);
      
      objNewTheta = computeObjf(X, M, B, Theta, Z, lambda1, lambda2);
      deltaTheta = abs(objNewTheta - objOldTheta);
      
      if(deltaTheta < threshold || round >= max_iteration) {
        loopFlag = true;
        objNew = objNewTheta;
        objOld = objNewB;
        if(verbose) {
          cout << "Round Theta is:" << roundTheta << endl;
        }
        roundTheta = 0;
      } else {
        continue;
      }
    }
    
    if(verbose) {
      int edgesTheta = computeEdgesTheta(Theta);
      int edgesB = computeEdgesB(B);
      cout << "edges of Theta:" << edgesTheta << endl;
      cout << "edges of B:" << edgesB << endl;
      int varyTheta = computeEdgesVaryTheta(Theta, ThetaOld);
      int varyB = computeEdgesVaryB(B, BOld);
      cout << "vary of edges Theta:" << varyTheta << endl;
      cout << "vary of edges B:" << varyB << endl;
      cout << "objNew " << objNew << endl;
      cout << "objOld " << objOld << endl;
    }
    
    // Judge termination
    delta = abs(objNew - objOld);
    if(verbose) {
      cout << "delta is : " << delta  << endl;
      cout << "lambda1 :" << lambda1 << endl;
      cout << "lambda2 :" << lambda2 << endl;
    }
    
    if(delta < threshold || round >= max_iteration) {
      if(verbose) {
        cout << "Rounds of iterations : " << round << endl;
        cout << "Iteration Success ~" << endl; 
      }
      break;
    }
  }
  
  // Compute the EBIC, return results
  EBIC = computeEBIC(objNew, B, Theta, p, q, n, lambda1, lambda2);
  if(verbose) {
    cout << "EBIC " << EBIC << endl;
  }
  
  VectorXd ThetaDiag = Theta.diagonal();
  MatrixXd OTU_OTU = (-1)*Theta.array() / (ThetaDiag * ThetaDiag.transpose()).array().sqrt();
  MatrixXd EF_OTU = B;
  
  return List::create(B, B0, Theta, Z, lambda1, lambda2, EBIC, EF_OTU, OTU_OTU);
} 

// [[Rcpp::export]]
List mLDM(NumericMatrix X, NumericMatrix M, double Z_mean = 1.0, int max_iteration = 2000, double threshold = 1e-4, int approx_num_Z = 10,
          int max_linesearch_Z = 30, int model_selection_num = 4, int approx_num_B = 10, int max_linesearch_B = 30, 
          int max_iteration_B = 500, double threshold_B = 1e-5, double delta1_threshold_B = 1e-4, double delta2_threshold_B = 0.9, double sy_threshold_B = 1e-6, 
          int max_iteration_B_coor = 20, double threshold_B_coor = 1e-6, double ratio1 = 0.6, double ratio2 = 0.9, bool verbose = false) {
  
  Map<MatrixXd> x(as<Map<MatrixXd> >(X));
  Map<MatrixXd> m(as<Map<MatrixXd> >(M));
  int n = x.rows();
  int p = x.cols();
  int q = m.cols();
  MatrixXd xRatio = toRatio(x);
  MatrixXd mScale = scale(m);
  
  // Set the initial values
  Environment stats("package:stats");
  Function cor = stats["cor"];
  MatrixXd corX(as<Map<MatrixXd> >(cor(Named("x", x), Named("method", "spearman"))));
  MatrixXd corM(as<Map<MatrixXd> >(cor(Named("x", xRatio), Named("y", mScale), Named("method", "spearman"))));
  // For Z
  MatrixXd Z_init = (x.array() + 1).array().log() + Z_mean;
  // For B
  MatrixXd B_init = corM.transpose();
  // For Theta
  MatrixXd Theta_cov = corX;
  double diagValue = 1;
  while(Theta_cov.determinant() < p) {
    Theta_cov.diagonal().array() += diagValue;
  }
  MatrixXd Theta_init = Theta_cov.inverse();
  // For B0
  VectorXd B0_init = Z_init.colwise().mean();
  
  // Set combinations of lambda1 and lambda2
  Environment base("package:base");
  Function unique = base["unique"];
  Function rep = base["rep"];
  Function quantile = stats["quantile"];
  Function seq = base["seq"];
  NumericVector corX_list = unique(rep(corX.array().abs()));
  NumericVector corM_list = unique(rep(corM.array().abs()));
  
  double lambda1_left = as<double>(quantile(corX_list, ratio1));
  double lambda1_right = as<double>(quantile(corX_list, ratio2));
  double lambda2_left = as<double>(quantile(corM_list, ratio1));
  double lambda2_right = as<double>(quantile(corM_list, ratio2));
  
  if(abs(ratio1 - ratio2) < 1e-3) {
    model_selection_num = 1;
  }
  
  NumericVector lambda1_list = seq(lambda1_left, lambda1_right, Named("length", model_selection_num));
  NumericVector lambda2_list = seq(lambda2_left, lambda2_right, Named("length", model_selection_num));
  
  int length1 = lambda1_list.size();
  int length2 = lambda2_list.size();
  
  // Record the minimum of EBIC
  double EBIC_min = 1e+20;
  // Record the optimal solution
  List LDM_result(0);
  // Record all solutions
  List LDM_result_all(length1*length2);
  int count = 0;
  // If the optimal lambda1 is found
  bool exist_best = false;
  // If is the optimal lambda1
  LogicalVector exist_per = rep(false, length1);
  int best_i = -1;
  
  // Load lbfgs and quic functions
  Environment lbfgs("package:lbfgs");
  Function lbfgsc = lbfgs["lbfgs"];
  Environment quic("package:QUIC");
  Function quicc = quic["QUIC"];
  
  // LoopFlag to control the order of optimaztion for B or Theta
  // flase optimize the Theta first
  // true optimize the B first
  bool loopFlag = false; 
  
  for(int j = 0; j < length2; ++j) {
    for(int i = 0; i < length1; ++i) {
      List solution(0);
      j > 0 ? loopFlag = true : loopFlag = false;
      //cout << j << " " << loopFlag << endl;
      if(!exist_best || exist_per[i]) {
        double lambda1 = lambda1_list[i];
        double lambda2 = lambda2_list[j];
        
        double solution_EBIC = 1e+30;
        
        solution = LognormalDirichletMultinomial(x, mScale, n, p, q, B_init, B0_init, Theta_init, Z_init, lambda1, lambda2, 
                      max_iteration, threshold, approx_num_Z, max_linesearch_Z, approx_num_B, 
                      max_linesearch_B, max_iteration_B, threshold_B, delta1_threshold_B, delta2_threshold_B, 
                      sy_threshold_B, max_iteration_B_coor, threshold_B_coor, verbose, lbfgsc, quicc, loopFlag);
        
        if(solution.size() != 0) {
          solution_EBIC = as<double>(solution[6]);
        }
        
        if(solution_EBIC > 0 && solution_EBIC < EBIC_min) {
          EBIC_min = as<double>(solution[6]);
          LDM_result = solution;
          best_i = i;
          if(j == 0 && i == 0) {
            B_init = as<MatrixXd>(LDM_result[0]);
            B0_init = as<VectorXd>(LDM_result[1]);
            Theta_init = as<MatrixXd>(LDM_result[2]);
            Z_init = as<MatrixXd>(LDM_result[3]);
          }
        }
      }
      
      LDM_result_all[count] = solution;
      count += 1;
    }
    
    if(best_i != -1) {
      exist_per[best_i] = true;
      exist_best = true;
    }
  }
  
  return List::create(Named("optimal") = LDM_result, Named("all") = LDM_result_all, Named("lambda1") = lambda1_list, Named("lambda2", lambda2_list));
}
