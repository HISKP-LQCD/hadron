#include <Rcpp.h>
#include <numeric>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>

using namespace Rcpp;

/**
 * @brief read 2*Nt doubles from a std::ifstream 
 *
 * @param ifs input stream
 * @param Nt unsigned integer: number of 'time slices' -> number of values to read from
 *                             file
 * @param real pre-allocated output vector for real part
 * @param imag pre-allocated output vector for imaginary part
 */
inline void read_correl(std::ifstream &ifs, const unsigned int Nt, std::vector<double> &real, std::vector<double> &imag)
{
  if( !ifs.good() ){
    stop("read_correl: input file stream not in a good state!");
  }
  for(unsigned int i=0; i<Nt; ++i)
  {
    ifs >> real[i];
    ifs >> imag[i];
  }
}

inline std::string make_key_2pt(unsigned int m1, unsigned int m2, unsigned int r1, unsigned int r2, std::string spin_comb)
{
  char ckey[100];
  snprintf(ckey, 100, "m1-%u_m2-%u_r1-%u_r2-%u_%s",
           m1, m2, r1, r2, spin_comb.c_str());
  return std::string(ckey);
}

inline void map_file(std::ifstream &ifs, std::map<std::string, std::iostream::pos_type> & filemap){
  if( !ifs.good() ){
    stop("map_file: input file stream not in a good state!");
  }

  std::string linebuf;
  std::string::size_type find_pos;

  // currently the key for a two-point function consists of four unsigned integer (two mass and two r indices)
  // and a generalisation of the reader will have to modify this as well as the key construction
  unsigned int key_components[4];

  while( ifs.good() ){
    std::getline(ifs, linebuf);
    // we search for commented lines
    find_pos = linebuf.find("#");
    if( find_pos != std::string::npos ){
      // in these commented lines, we extract either
      // the current set of mass / r parameter combinations
      // or the current spin combination
      // the line looks like so:
      // " # Contraction of S0_th0_m0_r0_ll ^ \dag and S0_th0_m0_r0_ll"
      std::string::size_type contr_pos = linebuf.find("Contraction");
      if( contr_pos != std::string::npos ){
        std::vector<char> lbcopy( linebuf.size() + 1 );
        lbcopy[ linebuf.size() ] = '\0';
        memcpy( lbcopy.data(), linebuf.c_str(), linebuf.size() );
        unsigned int key_components_counter = 0;
        char * token = strtok(lbcopy.data(), "_");
        while( key_components_counter != 4 | token != NULL ){
          if( token[0] == 'm' || token[0] == 'r'){
            key_components[key_components_counter] = atoi(token+1);
            key_components_counter++;
          }
          token = strtok(NULL, "_");
        }
      } else {
        // the line looks like so:
        // " # P5S0"
        std::string::size_type last_space_pos = linebuf.find_last_of(" ");
        std::string spin_comb = linebuf.substr(last_space_pos+1);
        // now we can build the key
        std::string key = make_key_2pt(key_components[0], key_components[2], key_components[1], key_components[3], spin_comb);
        // and store the position after the current newline as the starting point
        // of the present correlator
        filemap[ key ] = ifs.tellg();
      }
    }
  }
}

// [[Rcpp::export]]
NumericMatrix read_nissa_textcf_kernel(
    CharacterVector files_to_read,
    CharacterVector smear_combs_to_read,
    const unsigned int Nt,
    DataFrame combs_to_read)
{
  typedef NumericVector::iterator num_vec_iter;
  typedef CharacterVector::iterator char_vec_iter;

  const unsigned int n_correls = combs_to_read.nrows();
  const unsigned int n_smear_combs = smear_combs_to_read.size();
  const unsigned int n_files = files_to_read.size();

  // prepare some memory for output
  NumericMatrix cf_data(n_files, 2*Nt*n_correls*n_smear_combs );

  CharacterVector spin_combs = combs_to_read["spin_comb"];
  IntegerVector r1s = combs_to_read["r1"];
  IntegerVector r2s = combs_to_read["r2"];
  IntegerVector m1s = combs_to_read["m1"];
  IntegerVector m2s = combs_to_read["m2"];

  std::vector<double> realbuf(Nt);
  std::vector<double> imagbuf(Nt);

  for(unsigned int ifile = 0; ifile < n_files; ++ifile ){
    for(unsigned int ismear_comb = 0; ismear_comb < n_smear_combs; ++ismear_comb ){

      std::string filename = (std::string)(files_to_read[ifile]) + std::string("_") + 
                             (std::string)(smear_combs_to_read[ismear_comb]);

      std::ifstream ifs(filename.c_str(), std::ios::in);
      if( ! ifs.is_open() ){
        char message[200];
        snprintf(message, 200, "File %s could not be opened!", filename.c_str());
        stop(message);
      }
      std::map<std::string, std::iostream::pos_type> filemap;
      // create a map of the file which links a given observable
      // to a certain position in the file by reading the meta-data
      // this is a significant overhead incurred on a per-file basis,
      // but it is robust against changes in the file's structure
      map_file(ifs, filemap);
      // bring file back into good state
      ifs.clear();

      for( unsigned int icorrel = 0; icorrel < n_correls; ++icorrel ){
        // observable runs slowest, then smearing combination, then time, then real / imag
        // that should be the fastest way to split it into real and imaginary parts
        unsigned int out_idx = (icorrel * n_smear_combs + ismear_comb) * Nt * 2;

        // if we ever support other correlation functions from Nissa,
        // a simple point to generalise would be here in the way the key is
        // constructed
        std::string key = make_key_2pt((unsigned int)m1s[icorrel], 
                                       (unsigned int)m2s[icorrel], 
                                       (unsigned int)r1s[icorrel], 
                                       (unsigned int)r2s[icorrel],
                                       (std::string)spin_combs[icorrel]);
        if( filemap.count(key) != 1 ){
          char message[200];
          snprintf(message, 200, "correlator %s does not exist in file %s!", key.c_str(), filename.c_str());
          stop(message);
        }
        ifs.seekg( filemap[key] );
        read_correl(ifs, Nt, realbuf, imagbuf);

        for(int t=0; t<Nt; ++t){
          cf_data(ifile,out_idx  ) = realbuf[t];
          cf_data(ifile,out_idx+1) = imagbuf[t];
          out_idx += 2;
        }
      } // icorrel
    } // ismear_comb
  } // ifile
 return cf_data;
}
