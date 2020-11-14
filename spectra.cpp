
/*====================================================================

  Fourier transform of a real or complex function stored
    as two column data in a file using FFTW library.

  RUN COMMAND:
      ./spectra.exe <file_name> [<direction>] [<resol_factor>]
    where:
      <file_name> - two-column data file with real argument
        and real or complex value.
      <direction> - optional transform direction:
        'f' - forward, 'b' - backward. Forward is used by default.
      <resol_factor> - an optional integer factor to change result
        spectra resolution (default value resol_factor = 1).
        The time duration of the field loaded is incresed up
        to Tf = T*<resol_factor>. If <resol_factor> > 1
        the field is set to zero for time after T. Since
        frequency step of the result spectra is 2Pi/Tf, the
        spectra resolution multiplies by factor <resol_factor>.
        If <resol_factor> < 1 it is set to 1.

  ACKNOWLEDGEMENTS:

    Alexey D. Kondorskiy,
    P.N.Lebedev Physical Institute of the Russian Academy of Science.
    E-mail: kondorskiy@lebedev.ru, kondorskiy@gmail.com.


  Last modified: November 14, 2020.

====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <math.h>
#include <vector>
#include <cstdlib>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>
#include <fftw3.h>


/*--------------------------------------------------------------------
  File name prefixes:
--------------------------------------------------------------------*/

// For forward transform.
const std::string FWD_PREF = "dft_forward_";

// For backward transform.
const std::string BCK_PREF = "dft_backward_";


/*--------------------------------------------------------------------
  Complex numbers.
--------------------------------------------------------------------*/
const std::complex<double> ZERO(0.0, 0.0);
const std::complex<double> IRE(1.0, 0.0);
const std::complex<double> IIM(0.0, 1.0);


/*--------------------------------------------------------------------
  Check that the string is number.
--------------------------------------------------------------------*/
bool is_number(const std::string& s)
{
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it)) ++it;
  return !s.empty() && it == s.end();
}


/*--------------------------------------------------------------------
  Check that the string contains complex number.
--------------------------------------------------------------------*/
bool containsComplex(const std::string& s)
{
  return (s.find('(') != std::string::npos)
      && (s.find(',') != std::string::npos)
      && (s.find(')') != std::string::npos);
}


/*--------------------------------------------------------------------
  Remove complex number symbols.
--------------------------------------------------------------------*/
void cleanComplex(std::string& s)
{
  std::replace(s.begin(), s.end(), ',', ' ');
  std::replace(s.begin(), s.end(), '(', ' ');
  std::replace(s.begin(), s.end(), ')', ' ');
}


/*--------------------------------------------------------------------
  Read two column data from file.
--------------------------------------------------------------------*/
bool readTwoColumnData(const std::string &name,
  std::vector<double> &x, std::vector<std::complex<double> > &y)
{
  struct stat st;
  if (stat(name.c_str(), &st) != 0)
    return false;
  else {
    std::ifstream fin(name.c_str(), std::ios::in);
    while (!fin.eof()) {
      std::string str;
      std::getline(fin, str);
      double tmp;
      std::complex<double> ctmp;
      if (containsComplex(str)) {
        cleanComplex(str);
        std::stringstream ss(str);
        ss >> tmp;
        double y_re; ss >> y_re;
        double y_im; ss >> y_im;
        ctmp = IRE*y_re + IIM*y_im;
        if(!fin.eof()) x.push_back(tmp);
        if(!fin.eof()) y.push_back(ctmp);
      } else {
        std::stringstream ss(str);
        ss >> tmp;
        double y_re; ss >> y_re;
        ctmp = IRE*y_re;
        if(!fin.eof()) x.push_back(tmp);
        if(!fin.eof()) y.push_back(ctmp);
      }
      // std::cout << tmp << " " << ctmp << std::endl;
    }
    fin.close();
    return true;
  }
}


/*--------------------------------------------------------------------
  Make Fourier transform for given data file.
--------------------------------------------------------------------*/
void makeFourier(
  const std::string &file_name, // Input data file name.
  const double &resol_factor,   // Factor to increase resolution.
  const int &dir)           // Transformation direction
{
  std::vector<double> x;
  std::vector<std::complex<double> > y;
  if (!readTwoColumnData(file_name, x, y)) return;

  // Size of extended array.
  int n_time = x.size();
  int n_tot = n_time;
  if (resol_factor > 1.0) {
    n_tot = int(n_time*resol_factor);
    if (n_tot < n_time)
      n_tot = n_time;
  }

  // Make n_tot to be even.
  if (n_tot % 2 == 1)
    n_tot = n_tot - 1;

  // Create pointers and arrays.
  double *time, *freq;
  time = new double[n_tot];
  freq = new double[n_tot];
  fftw_complex *in, *out;
  in  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n_tot);
  out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n_tot);

  // Load data.
  double tstep = 0.0;
  for (int i = 0; i < n_tot; ++i) {

    if (i < n_time) {

      // Load data from file.
      time[i] = x[i];
      in[i][0] = real(y[i]);
      in[i][1] = imag(y[i]);

      // Prepare time grid.
      if (i > 0) {

        // Detect time step.
        if (tstep == 0.0) {
          tstep = time[1] - time[0];
          if(tstep <= 0.0) {
            std::cout << "Bad argument step in file " << file_name
              << "; t[1]-t[0] = " << time[1] - time[0] << std::endl;
            fftw_free(in);
            fftw_free(out);
            delete [] time, freq;
            return;
          }
        }

        // Check regular grid.
        double treg = time[0] + tstep*i;
        if ( std::abs(time[i] - treg) > (5.0e-2*tstep) ) {
          std::cout << "Irregular arguments in file. " << file_name
            <<  ". Index = " << i << ", ht = " << tstep
            << " t[i] - t[i-1] = " << std::abs(time[i] - treg)
            << std::endl;
          fftw_free(in);
          fftw_free(out);
          delete [] time, freq;
          return;
        }
      }
    } else {

      // Check time step.
      if(tstep <= 0.0) {
        std::cout << "Bad argument step in file " << file_name
          << "; t[1]-t[0] = " << time[1] - time[0] << std::endl;
        fftw_free(in);
        fftw_free(out);
        delete [] time, freq;
        return;
      }

      // Extrapolate data.
      time[i] = time[0] + tstep*i;
      in[i][0] = 0.0;
      in[i][1] = 0.0;
    }
  }

  // Prepare FFTW plans.
  fftw_plan p;
  if (dir == 1)
    p = fftw_plan_dft_1d(n_tot, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  else
    p = fftw_plan_dft_1d(n_tot, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  // Run FFTW.
  fftw_execute(p);

  // Parse the results.
  double x_step = (time[n_tot - 1] - time[0])/(n_tot - 1);
  double f_step = 2*M_PI/(n_tot*x_step);
  double norm = x_step/sqrt(2*M_PI);
  double sgn;
  if (dir == 1)
    sgn =  1.0;
  else
    sgn = -1.0;

  std::complex<double> *tmp;
  tmp = new std::complex<double>[n_tot];

  for(int i = 0; i < n_tot; ++i)
    tmp[i] = IRE*out[i][0] + IIM*out[i][1];

  for(int i = 0; i < n_tot; ++i)
    if(i < n_tot/2) {
      freq[i + n_tot/2] = i*f_step;
      std::complex<double> ctmp = tmp[i]*norm
        *exp(IIM*sgn*time[0]*freq[i + n_tot/2]);
      out[i + n_tot/2][0] = std::real(ctmp);
      out[i + n_tot/2][1] = std::imag(ctmp);
    } else {
      freq[i - n_tot/2] = (i - n_tot)*f_step;
      std::complex<double> ctmp = tmp[i]*norm
        *exp(IIM*sgn*time[0]*freq[i - n_tot/2]);
      out[i - n_tot/2][0] = std::real(ctmp);
      out[i - n_tot/2][1] = std::imag(ctmp);
    }
  delete [] tmp;

  // Print out the result.
  std::string fout_name;
  if (dir == 1)
    fout_name = BCK_PREF + file_name;
  else
    fout_name = FWD_PREF + file_name;
  std::ofstream fout(fout_name.c_str(), std::ios::out);
  // fout << std::scientific;
  for(int i = 0; i < n_tot; ++i) {
    std::complex<double> ctmp = IRE*out[i][0] + IIM*out[i][1];
    fout << freq[i] << " " << ctmp << std::endl;
    }
  fout.close();

  // Print separately real and imaginary parts.
  std::string fout_name_re = "re-" + fout_name;
  std::string fout_name_im = "im-" + fout_name;
  std::ofstream fout_re(fout_name_re.c_str(), std::ios::out);
  std::ofstream fout_im(fout_name_im.c_str(), std::ios::out);
  for(int i = 0; i < n_tot; ++i) {
    fout_re << freq[i] << " " << out[i][0] << std::endl;
    fout_im << freq[i] << " " << out[i][1] << std::endl;
    }
  fout_re.close();
  fout_im.close();

  // Clean up.
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
  delete [] time, freq;
}


/*********************************************************************
  MAIN
*********************************************************************/
int main(int argc, char **argv)
{
  // Read file name.
  std::string file_name;
  if (argc > 1) {
    std::string arg1(argv[1]);
    file_name = arg1;
  } else {
    return 0;
  }

  // Direction or resolution.
  int dir = 0;
  double resol_factor = 1.0;
  if (argc > 2) {
    std::string arg2(argv[2]);
    if (is_number(arg2)) {
      resol_factor = strtod(argv[2], NULL);
    } else {
      if ((*argv[2] == 'b') || (*argv[2] == 'B'))
        dir = 1;
      if(argc > 3)
        resol_factor = strtod(argv[3], NULL);
      if(resol_factor < 1.0)
        resol_factor = 1.0;
    }
  }

  // RUN.
  makeFourier(file_name, resol_factor, dir);

  return 0;
}


//====================================================================
