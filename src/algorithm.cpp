#include <RcppArmadillo.h>

RcppExport SEXP doEM_REML( SEXP curTheta_, SEXP curK_, SEXP y_, SEXP X_, SEXP yy_, SEXP Xy_, SEXP XX_, SEXP phi_, SEXP phiInv_, SEXP trphiInv_, SEXP S_, SEXP q_, SEXP threshold_){

  double curTheta = Rcpp::as<double>(curTheta_);
  double curK = Rcpp::as<double>(curK_);
  arma::colvec y = Rcpp::as<arma::colvec>(y_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  double yy = Rcpp::as<double>(yy_);
  arma::colvec Xy = Rcpp::as<arma::colvec>(Xy_);
  arma::mat XX = Rcpp::as<arma::mat>(XX_);
  arma::mat phi = Rcpp::as<arma::mat>(phi_);
  arma::mat phiInv = Rcpp::as<arma::mat>(phiInv_);
  double trphiInv = Rcpp::as<double>(trphiInv_);
  int S = Rcpp::as<int>(S_);
  int q = Rcpp::as<int>(q_);
  double threshold = Rcpp::as<double>(threshold_);

  arma::mat czz;
  arma::mat czzInv;
  double logL = std::numeric_limits<double>::infinity();
  double yVIy;
  double yPy;
  double newTheta;
  double newK;
  double newlogL;
  arma::colvec XVIy;
  arma::mat XVIX;
  arma::mat XVIXInv;
  arma::mat c22;
  arma::colvec blue;
  arma::colvec predicted;
  arma::colvec blup;
  arma::mat I = arma::eye<arma::mat>(S, S);
  arma::mat tX = arma::trans(X);
  double val1, val2, sign;
  double curK2;
  arma::vec res(2);

  int iter = 0;

  while(1){
    iter++;
    curK2 = pow(curK,2);
    czz = 1 / curK * I + phiInv;
    czzInv = arma::inv(czz);
    yVIy = yy / curK - arma::as_scalar(arma::trans(y)*czzInv*y) / curK2;
    XVIy = Xy / curK - tX * czzInv * y / curK2;
    XVIX = XX / curK - tX * czzInv * X / curK2;
    XVIXInv = arma::inv(XVIX);
    c22 = czzInv + czzInv * X * XVIXInv * tX * czzInv / curK2;

    blue = XVIXInv * XVIy;
    predicted = X * blue;
    blup = czzInv * (y / curK - predicted / curK);

    yPy = yVIy - arma::as_scalar( arma::trans(XVIy) * XVIXInv * XVIy );

    newTheta = yPy / (S-q);
    newK = curK  - ( trphiInv - arma::sum( arma::diagvec( phiInv * c22 * phiInv ) ) ) * curK2 / S + arma::sum( arma::square(y - predicted - blup) ) / curTheta / S;

    arma::log_det(val1, sign, XVIX);
    arma::log_det(val2, sign, phi + I * curK);

    newlogL = -0.5 * ( val1 +  val2 + (S - q) * log( curTheta ) + yPy / curTheta );

    //Rprintf("EMs: yVIy %f, yPy %f, newTheta %f, newK %f, newlogL %f \n", yVIy, yPy, newTheta, newK, newlogL);

    //    Rprintf("%f, %f  ", logL, newlogL);

    if(fabs(logL - newlogL) < threshold){
      res[0] = newTheta;
      res[1] = newK;
      return Rcpp::wrap(res);
    }

    logL = newlogL;
    curTheta = newTheta;
    curK = newK;

  }
}


RcppExport SEXP doEM_ML(SEXP snp_, SEXP cursig2g_, SEXP curK_, SEXP y_, SEXP X_, SEXP yy_, SEXP Xy_, SEXP XX_, SEXP phi_, SEXP phiInv_, SEXP trphiInv_, SEXP S_, SEXP q_, SEXP threshold_){

  arma::colvec snp = Rcpp::as<arma::colvec>(snp_);
  double cursig2g = Rcpp::as<double>(cursig2g_);
  double curK = Rcpp::as<double>(curK_);
  arma::colvec y = Rcpp::as<arma::colvec>(y_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  double yy = Rcpp::as<double>(yy_);
  arma::colvec Xy = Rcpp::as<arma::colvec>(Xy_);
  arma::mat XX = Rcpp::as<arma::mat>(XX_);
  arma::mat phi = Rcpp::as<arma::mat>(phi_);
  arma::mat phiInv = Rcpp::as<arma::mat>(phiInv_);
  double trphiInv = Rcpp::as<double>(trphiInv_);
  int S = Rcpp::as<int>(S_);
  int q = Rcpp::as<int>(q_);
  double threshold = Rcpp::as<double>(threshold_);

  double cursig2;
  double sig2g;
  double sig2;
  arma::mat czz;
  arma::mat czzInv;
  arma::mat V;
  arma::mat W;
  double logL = std::numeric_limits<double>::infinity();
  double bet;
  arma::colvec XVIy;
  arma::mat XVIX;
  arma::mat XVIXInv;
  arma::colvec blue;
  arma::colvec predicted;
  arma::colvec blup;
  arma::colvec err;
  arma::mat I = arma::eye<arma::mat>(S, S);
  arma::mat tX = arma::trans(X);
  double val, sign;
  double curK2;
  arma::colvec tmp;
  arma::vec res(2);

  cursig2 = curK * cursig2g;

  while(1){
    curK2 = pow(curK,2);
    V = cursig2 * I + cursig2g * phi;
    W = arma::inv(V);

    curK = cursig2 / cursig2g;
    czz = 1 / curK * I + phiInv;
    czzInv = arma::inv(czz);
    XVIy = Xy / curK - tX * czzInv * y / curK2;
    XVIX = XX / curK - tX * czzInv * X / curK2;
    XVIXInv = arma::inv(XVIX);

    blue = XVIXInv * XVIy;
    predicted = X *blue;
    blup = czzInv * ( y / curK - predicted / curK);
    err = y - predicted - blup;
    sig2 = cursig2 - pow(cursig2, 2) * arma::mean( arma::diagvec(W) ) + arma::mean( arma::square(err) );
    sig2g = cursig2g - pow(cursig2g, 2) * arma::mean( arma::diagvec( W * phi) ) + arma::as_scalar( arma::trans(blup) * phiInv * blup ) / S;
    bet = blue[blue.n_elem];
    tmp = y - snp*bet;
    arma::log_det(val, sign, sig2g * phi + sig2 * I);
    logL = - val - arma::as_scalar( arma::trans( tmp ) * W * tmp );

    if( sqrt( pow(cursig2g - sig2g, 2) + pow(cursig2 - sig2, 2) ) < threshold ){
      res[0] = sig2g;
      res[1] = sig2/sig2g;
      return Rcpp::wrap(res);
    }
    cursig2g = sig2g;
    cursig2 = sig2;
  }
}

RcppExport SEXP doNR_ML (SEXP curTheta_, SEXP curK_, SEXP y_, SEXP X_, SEXP phi_, SEXP S_, SEXP itrmax_, SEXP threshold_){

  double curTheta = Rcpp::as<double>(curTheta_);
  double curK =Rcpp::as<double>(curK_);
  arma::colvec y = Rcpp::as<arma::colvec>(y_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat phi = Rcpp::as<arma::mat>(phi_);
  int S = Rcpp::as<int>(S_);
  int itrmax = Rcpp::as<int>(itrmax_);
  double threshold = Rcpp::as<double>(threshold_);

  arma::mat V;
  arma::mat VI;
  arma::mat XVI;
  arma::mat XVIX;
  arma::mat XVIXInv;
  arma::colvec XVIy;
  arma::colvec beta;
  arma::colvec predicted;
  arma::colvec resid;
  arma::colvec blup;
  arma::mat P;
  arma::colvec Py;
  double sig2;
  double rVIr;
  double logL;
  double yPPy;
  double yPy;
  double trVI;
  double uTheta;
  double uK;
  double I11;
  double I12;
  double I22;
  double val;
  double sign;
   arma::mat II(2, 2);
   arma::mat IIInv;

  arma::colvec tmp1(2);
  arma::colvec tmp2(2);
  arma::colvec updates(2);

  arma::mat I = arma::eye<arma::mat>(S, S);
  int itr = 1;


  while(1){
    V = curK * I + phi;
    VI = arma::inv(V);
    XVI = arma::trans(X) * VI;
    XVIX = XVI * X;
    XVIXInv = arma::inv(XVIX);
    P = VI - arma::trans(XVI) * XVIXInv * XVI;
    Py = P * y;
    yPPy = arma::as_scalar( arma::trans(Py) * Py);
    yPy = arma::as_scalar( arma::trans(y) * Py);
    trVI = arma::as_scalar( arma::sum( arma::diagvec(VI) ) );

    uTheta = -0.5 * S / curTheta + 0.5 * yPy / pow(curTheta, 2);
    uK = -0.5 * trVI + 0.5 / curTheta * yPPy;
    I11 = S / (2 * pow(curTheta, 2) ) - 1 / pow(curTheta, 3) * yPy;
    I12 = -1 / (2 * pow(curTheta, 2)) * yPPy;
    I22 = 0.5 * arma::sum( arma::diagvec( VI * VI) ) - 1 / curTheta * arma::as_scalar( arma::trans( Py ) * P * Py );
    II(0, 0) = I11;
    II(0, 1) = I12;
    II(1, 0) = I12;
    II(1, 1) = I22;

    IIInv = arma::inv(II);

    tmp1[0] = curTheta;
    tmp1[1] = curK;
    tmp2 = IIInv * tmp1;
    updates = tmp1 - tmp2;
    curTheta = updates[0];
    curK = updates[1];

    if( sqrt( arma::sum( arma::pow( updates - tmp1, 2) ) ) < threshold ){
      V = curK * I + phi;
      VI = arma::inv(V);
      XVI = arma::trans(X) * VI;
      XVIX = XVI * X;
      XVIy = XVI * y;
      beta = arma::inv( XVIX ) * XVIy;
      predicted = X * beta;
      resid =  y - predicted;
      rVIr = arma::as_scalar(arma::trans( resid ) * VI * resid);
      sig2 = curK * curTheta;
      arma::log_det(val, sign, V);
      logL = -0.5*( val + S * log(rVIr) );
      blup = phi * VI * resid;
      return Rcpp::List::create(Rcpp::Named("para", beta), Rcpp::Named("sig2g", curTheta), Rcpp::Named("sig2", sig2), Rcpp::Named("logL", logL));
    }

    if(itr > itrmax){
      Rprintf("It does not converge!!!");
      return Rcpp::wrap(NA_REAL);
    }
    // Rprintf("itr: %d, sig2: %f, sig2g: %f\n", itr, curTheta * curK, curTheta);
    itr++;
  }
}


RcppExport SEXP doAI_REML(SEXP curTheta_, SEXP curK_, SEXP y_, SEXP X_, SEXP yy_, SEXP Xy_, SEXP XX_, SEXP phi_, SEXP phiInv_, SEXP trphiInv_, SEXP S_, SEXP q_, SEXP itrmax_, SEXP threshold_){
  double curTheta = Rcpp::as<double>(curTheta_);
  double curK = Rcpp::as<double>(curK_);
  arma::colvec y = Rcpp::as<arma::colvec>(y_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  double yy = Rcpp::as<double>(yy_);
  arma::colvec Xy = Rcpp::as<arma::colvec>(Xy_);
  arma::mat XX = Rcpp::as<arma::mat>(XX_);
  arma::mat phi = Rcpp::as<arma::mat>(phi_);
  arma::mat phiInv = Rcpp::as<arma::mat>(phiInv_);
  double trphiInv = Rcpp::as<double>(trphiInv_);
  int S = Rcpp::as<int>(S_);
  int q = Rcpp::as<int>(q_);
  int itrmax = Rcpp::as<int>(itrmax_);
  double threshold = Rcpp::as<double>(threshold_);

  arma::mat czz;
  arma::mat czzInv;
  double yVIy;
  arma::colvec XVIy;
  arma::mat XVIX;
  arma::mat XVIXInv;
  double yPy;
  arma::colvec Py;
  arma::mat Q;
  arma::mat QQ;
  arma::mat QX;
  arma::mat QVIQ;
  arma::mat QVIX;
  arma::mat QPQ;
  double Utheta;
  arma::mat c22;
  double UK;
  arma::colvec updates(2);
  arma::colvec tmp1(2);
  arma::colvec tmp2(2);

  arma::mat tX = arma::trans(X);
  arma::mat tQ;
  double logL;

  arma::colvec blue;
  arma::colvec predicted;
  arma::colvec blup;
  arma::mat I = arma::eye<arma::mat>(S, S);

  double val1, val2, sign;
  double curK2;
  arma::vec res(2);
  int itr = 1;

  while(1){
    curK2 = pow(curK,2);
    czz = 1 / curK * I + phiInv;
    czzInv = arma::inv(czz);
    XVIy = Xy / curK - tX * czzInv * y / curK2;
    XVIX = XX / curK - tX * czzInv * X / curK2;
    XVIXInv = arma::inv(XVIX);
    blue = XVIXInv * XVIy;
    predicted = X * blue;
    blup = czzInv * (y / curK - predicted / curK);
    Py = ( y - predicted - blup ) / curK;
    Q =  (y - predicted) / curTheta;
    Q.insert_rows(Q.n_elem, Py);
    Q.reshape(S, Q.n_elem / S);
    tQ = arma::trans(Q);
    QQ = tQ * Q;
    QX = tQ * X;
    QVIQ = QQ / curK - tQ * czzInv * Q / curK2;
    QVIX = QX / curK - tQ * czzInv * X / curK2;
    QPQ = QVIQ - QVIX * XVIXInv * arma::trans(QVIX);
    Utheta = - 0.5 * ( (S-q) / curTheta - QPQ(0, 0) );
    c22 = czzInv + czzInv * X * XVIXInv * tX * czzInv / curK2;
    UK = -0.5 * ( trphiInv - arma::sum( arma::diagvec( phiInv * c22 * phiInv ) ) - arma::as_scalar( arma::trans(Py) * Py ) / curTheta );
    tmp1[0] = curTheta;
    tmp1[1] = curK;
    tmp2[0] = Utheta;
    tmp2[1] = UK;
    updates = tmp1 + 2 * arma::vec( curTheta * arma::inv(QPQ) * tmp2 );
    if(updates[0] <=0 || updates[1] <= 0){
      Rprintf("Negative definite and EM is applied");
      updates[0] = QPQ(0,0) * pow(curTheta,2) / (S-q);
      updates[1] = curK - (trphiInv - arma::sum( arma::diagvec( phiInv * c22 * phiInv ) ) ) * curK2 / S + arma::sum( arma::pow( y - predicted - blup, 2) ) / curTheta / S;
    }

    if(sqrt( arma::sum( arma::pow(updates - tmp1, 2) ) ) < threshold ){
      curTheta = updates[0];
      curK = updates[1];
      czz = 1/curK * I + phiInv;
      czzInv = arma::inv(czz);
      yVIy = yy / curK - arma::as_scalar( arma::trans(y) * czzInv * y ) / curK2;
      XVIy = Xy / curK - tX * czzInv * y / curK2;
      XVIX = XX / curK - tX * czzInv * X / curK2;
      XVIXInv = arma::inv(XVIX);
      //yPy = yVIy - arma::scalar( arma::trans(XVIy) * XVIXInv * XVIy );

      arma::log_det(val1, sign, XVIX);
      arma::log_det(val2, sign, phi + I * curK);
      logL = -0.5 * ( val1 +  val2 + (S - q) * log( curTheta ) + QPQ(0, 0) * curTheta );
      blue = XVIXInv * XVIy;
      predicted = X * blue;
      blup = czzInv * (y / curK - predicted / curK);
      return Rcpp::List::create(Rcpp::Named("para",blue), Rcpp::Named("sig2g",curTheta), Rcpp::Named("sig2", curK*curTheta), Rcpp::Named("logL", logL) );
    }
    else{
      //Rprintf("itr: %d, sig2: %f, sig2g, %f \n", itr, updates[0] * updates[1], updates[0]);
      curTheta = updates[0];
      curK = updates[1];
    }
    if(itr > itrmax){
      Rprintf("It does not converge!!!");
      return Rcpp::wrap(NA_REAL);
    }
    itr++;
     }
}

