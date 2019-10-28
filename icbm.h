// ICBM - integrated combined baseline modification
// ICBM corrects baseline jumps of JMA acceleration waveform data by a segmented least squares fit on the velocity trace
// variable number of parameters is implemented
// requires -larmadillo when compiling

// COPYRIGHT 2019 Sebastian von Specht

// This software is distributed under the MIT license. See the license file in this repository.

#ifndef ICBM_H
#define ICBM_H

#include <cmath>
#include <fstream>
// #include <sstream>
#include <iostream>
#include <vector>
#include <armadillo>
#include <limits>

bool compare(double a, double b) {
  return (a < b);
}

std::vector<std::vector<double>> icbm(std::vector<std::vector<double>> x, int Nsmin = 1, int Nsmax = 6, double thrshld = 1e-5, int nmin = 500) {

  if (x.size() != 3) {
    std::cerr << "Wrong number of traces.\nTerminate.\n";
    // return VECTOR;
  }

  // std::ifstream ifs;
  std::ofstream ofs;


  int n = x[0].size();

  // ofs.open("diffinttest");

  double c[3];

  for (int i = 0; i < n; ++i) {

    double mc {0.};
    for (int j = 0; j < 3; ++j) {
      if (i == 0)
        c[j] = (x[j][i]);
      else
        c[j] += (x[j][i]);

      mc += c[j]/3.;

      // ofs << c[j] << " ";
    }

    double vc {0.};

    for (int j = 0; j < 3; ++j)
      vc += pow(c[j] - mc, 2.);

    vc /= 2.;
    // ofs << mc << " " << vc << "\n";
  }
  // ofs.close();

  double v1[3], ym[3] {0.,0.,0.}, tm[3] {0.,0.,0.};
  arma::vec y(3*n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == 0)
        // v1[j] = (x[j][i]); // rectangle rule
        v1[j] = 0.; // trapezoidal rule
      else
        // v1[j] += (x[j][i]); // rectangle rule
        v1[j] += (x[j][i]+x[j][i-1])/2.; // trapezoidal rule

      y(3*i+j) = v1[j];

      ym[j] += v1[j];
      tm[j] += (double)i;
    }
  }

  for (int j = 0; j < 3; ++j) {
    ym[j] /= (double)n;
    tm[j] /= (double)n;
  }

  double vxy[3] {0.,0.,0.}, vx[3] {0.,0.,0.};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < 3; ++j) {
      vxy[j] += (y(3*i+j) - ym[j])*((double)i - tm[j]);
      vx[j] += pow((double)i - tm[j],2.);
    }
  }

  double p0[3], p1[3];
  for (int j = 0; j < 3; ++j) {
    p1[j] = vxy[j] / vx[j];
    p0[j] = ym[j] - p1[j]*tm[j];

    std::cout << p0[j] << " " << p1[j] << "\n";
  }


  ofs.open("testaic");
  // std::ofstream ofs2("ibcparam");

  std::vector < double > vr;

  double minAIC {std::numeric_limits < double > :: infinity()};
  double minBIC {std::numeric_limits < double > :: infinity()};

  std::vector < std::vector < double > > x_corr(3, std::vector < double > (n));
  std::vector < std::vector < double > > v_corr(3, std::vector < double > (n));
  std::vector < std::vector < double > > v_mod(3, std::vector < double > (n));


  for (int Ns = Nsmin; Ns < Nsmax; ++Ns) {

    int N {Ns};

    arma::vec ym(3*n);
    arma::vec r(3*n), wl1(3*n);
    r.zeros();

    arma::vec p(6+4*(N-1), arma::fill::zeros);

    double sumres {0.}, sumresold {1.};
    int j {0};

    // std::ofstream ofs("testiter");


    while (((fabs(1.-sumres/sumresold) > thrshld) & (j < nmin))) {

      if (j == 0) {
        for (int i = 0; i < 3; ++i)
          p(i) = p0[i];
        for (int i = 3; i < 6; ++i)
          p(i) = p1[i-3];

        for (int i = 0; i < N-1; ++i) {
          double pt = (double)nmin/2.+ round((double)(i*(n-nmin)/(double)(N-2)));
          p(6+4*i+3) = pt;//round((double)((i)*n)/(double)(N-2));
        }
        // p.print();
      }

      // the Jacobian
      arma::mat J(3*n,6+4*(N-1), arma::fill::zeros);
      arma::umat Wind(2,3*n);

      for (int i = 0; i < n; ++i) {

        double t = (double) i;

        for (int k = 0; k < 3; ++k) {

          ym(3*i+k) = p(k) + p(3+k)*t;

          for (int l = 0; l < N-1; ++l) {
            ym(3*i+k) += t >= p(6+4*l+3) ? p(6+4*l+k)*(t - p(6+4*l+3)) : 0.;

            J(i*3+k,k) = 1.;
            J(i*3+k,3+k) = t;
            J(i*3+k,6+4*l+k) = t >= p(6+4*l+3) ? t - p(6+4*l+3) : 0.;
            J(i*3+k,6+4*l+3) = t >= p(6+4*l+3) ? -p(6+4*l+k) : 0.;
          }
          r(3*i+k) = y(3*i+k) - ym(3*i+k);
          wl1(3*i+k) = std::max(1./fabs(r(3*i+k)),thrshld);
          Wind(0,3*i+k) = Wind(1,3*i+k) = 3*i+k;
        }

        // ofs << i << " " << y(3*i) << " " << ym(3*i) << " " << j << "\n";
      }
      // ofs << "\n";

      // p += solve(J,r);

      // p.print();
      // std::cout << "===\n";


      arma::mat I(6+4*(N-1),6+4*(N-1),arma::fill::eye);
      arma::sp_mat W(Wind,wl1);

      arma::mat JtJ = J.t() * W * J;
      I.diag() = JtJ.diag()+1.; // the +one is for numeric reason, cause in the 1st step no segments are assumed, hence the time derivatives vanish and the Gramian has zero entries on the diagonal
      p += (JtJ + 1e-6*I).i() * J.t() *W* r;

      std::vector<unsigned int> idxvector;
      for (int k = 0; k < 6; ++k) {
        idxvector.push_back((unsigned int)k);
      }

      int Nnew {1};

      // round time steps and shift them by one to keep jump in correct location (because of t-T=0, but h(0)=1)
      for (int l = 0; l < N-1; ++l) {
        // p(6+4*l+3) = round(p(6+4*l+3)+1.);
        // std::cout << p(6+4*l+3) << "\n";
      }
      // std::cout << "==\n";


      for (int l = 0; l < N-1; ++l) {

        for (int m = l+1; m < N-1; ++m) {
          if (fabs(p(6+4*l+3) - p(6+4*m+3)) < (double)nmin) {
            for (int k = 0; k < 3; ++k) {

              // double dtl {0.}, dtm {0.};
              //
              // if (l < N-2)
              //   dtl = p(6+4*l+4)

              p(6+4*m+k) = (p(6+4*m+k) + p(6+4*l+k)) / 2.;
              p(6+4*l+k) = 0.; // This line and the next are necessary, as the lth parameter set is flushed with the next 'if' sequence
            }

            p(6+4*m+3) = round((p(6+4*m+3) + p(6+4*l+3)) / 2.);
            p(6+4*l+3) = (double)(2*n);
          }
        }



        if (p(6+4*l+3) < 0.) {
          for (int k = 0; k < 3; ++k) {
            p(k) += p(6+4*l+k);
            p(6+4*l+k) = 0.;
          }
        }
        else if (p(6+4*l+3) > (double)(n-1)) {
          for (int k = 0; k < 3; ++k)
            p(6+4*l+k) = 0.;
        }

        double sump {0.};
        for (int k = 0; k < 3; ++k)
          sump += fabs(p(6+4*l+k));

        if (sump < 1e-6)
          for (int k = 0; k < 3; ++k)
            p(6+4*l+k) = 0.;
        else {
          ++Nnew;
          for (int k = 0; k < 3; ++k)
            idxvector.push_back((unsigned int)(6+4*l+k));
          idxvector.push_back((unsigned int)(6+4*l+3));
        }
      }

      arma::uvec idxlist(idxvector.size());
      for (int l = 0; l < idxvector.size(); ++l)
        idxlist(l) = idxvector[l];

      // p.print();

      p = p.elem(idxlist);
      N = Nnew;
      // std::cout << "\n=====\n";
      //
      // p.print();


      sumresold = sumres;
      sumres = 0.;
      for (double i = 0; i < 3*n; ++i)
        sumres += pow(r(i),2.);

      ++j;
    }

    std::cout << "\n=====" << Ns << "=====\nPARAMETERS\n";
    for (int k = 0; k < 3; ++k)
      std::cout << p(k) << " ";
    std::cout << "\n";
    for (int k = 0; k < 3; ++k)
      std::cout << p(3+k) << " ";
    std::cout << "\n";
    for (int l = 0; l < N-1; ++l) {
      for (int k = 0; k < 3; ++k)
        std::cout << p(6+4*l+k) << " ";
      std::cout << p(6+4*l+3) << "\n";
      // ofs2 << p(6+4*l+3) << " ";
    }
    std::cout << "\n";
    // ofs2 << "\n";

    vr.push_back(0.);

    for (int i = 0; i < 3*n; ++i)
      vr.back() += pow(y(i)-ym(i), 2.);

    vr.back() /= (double) (3*n);

    // the residual sum (variance) is a free parameter as well

    double AIC = (double) (2*(6+4*(N-1)+1)) + (double) (3*n) * log(vr.back());
    double BIC = log((double)(3*n)) * (double)(6+4*(N-1)+1) + (double) (3*n) * log(vr.back());

    if (BIC < minBIC) {
      std::cout << "best model with N = " << N << "\n";
      minAIC = AIC;
      minBIC = BIC;
      for (int i = 0; i < n; ++i) {
        double t = (double) i;
        for (int k = 0; k < 3; ++k) {
          x_corr[k][i] = x[k][i] - p(3+k);
          v_mod[k][i] =  p(k) + p(3+k)*t;
          for (int l = 0; l < N-1; ++l) {
            x_corr[k][i] -= t >= p(6+4*l+3) ? p(6+4*l+k) : 0.;
            v_mod[k][i] += t >= p(6+4*l+3) ? p(6+4*l+k)*(t - p(6+4*l+3)) : 0.;
          }
          v_corr[k][i] = y(3*i+k) - v_mod[k][i];
        }
      }
    }

    std::cout << N << " " << vr.back() << " " << AIC << " " << BIC << "\n";
    ofs << N << " " << vr.back() << " " << AIC << " " << BIC << "\n";
  }

  ofs.close();

  ofs.open("testvcorr");
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < 3; ++k)
      ofs << v_corr[k][i] << " ";
    for (int k = 0; k < 3; ++k)
      ofs << v_mod[k][i] << " ";
    ofs << "\n";
  }

  ofs.close();

  return x_corr;
}

#endif
