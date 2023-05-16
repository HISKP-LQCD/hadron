#include <Rcpp.h>

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace Rcpp;

// corrtypes: what type of information is read from the file
// 1=2pt: m1, m2, r1, r2, spin_info
// 2=newcorr: op1, op2, spin_info

/**
 * @brief read 2*nts doubles from a std::ifstream 
 *
 * @param ifs input stream
 * @param nts unsigned integer: number of 'time slices' -> number of values to read from
 *                             file
 * @param real pre-allocated output vector for real part
 * @param imag pre-allocated output vector for imaginary part
 */
inline void read_correl(std::ifstream &ifs, const unsigned int nts, std::vector<double> &real, std::vector<double> &imag)
{
  char message[200];
  for(unsigned int i=0; i<nts; ++i)
  {
    if( !ifs.good() ){
      snprintf(message, 200, "read_correl: input file steram is not in a good state at time slice %u", i);
      stop(message);
    }
    ifs >> real[i];
    ifs >> imag[i];
  }
}

/**
 * @brief make a key string 
 *
 * @param m1_idx mass index of the backward (daggered) propagator in the contraction
 * @param m2_idx mass index of the forward propagator in the contraction
 * @param r1_idx r index (0=+, 1=-) for the backward (daggered) propagator in the contraction
 * @param r2_idx r index of the forward propagator in the contraction
 * @param spin_comb name of the spin combination, for example "P5P5" or "V1A0"
 *
 * @return key of the form "m1_idx-%u_m2_idx-%u_r1_idx-%u_r2_idx-%u_%s"
 */
inline std::string make_key_2pt(unsigned int m1_idx, unsigned int m2_idx, unsigned int r1_idx, unsigned int r2_idx, std::string spin_comb)
{
  char ckey[100];
  snprintf(ckey, 100, "m1_idx-%u_m2_idx-%u_r1_idx-%u_r2_idx-%u_%s",
           m1_idx, m2_idx, r1_idx, r2_idx, spin_comb.c_str());
  return std::string(ckey);
}

/**
 * @brief make a key string of two operators
 *
 * @param op1_idx name of the first correlator, for example C_H
 * @param op2_idx name of the second correlator, for example H_H_S_H
 * @param spin_comb name of the spin combination, for example "P5P5" or "V1A0"
 *
 * @return key of the form "op1_%s_op2_%s_%s"
 */
inline std::string make_key_newcorr(std::string op1_idx, std::string op2_idx, std::string spin_comb)
{
  char ckey[100];
  snprintf(ckey, 100, "op1_%s_op2_%s_%s",
           op1_idx.c_str(), op2_idx.c_str(), spin_comb.c_str());
  return std::string(ckey);
}

/**
 * @brief read correlator meta-data and map to file positions
 * @description Nissa text correlator files contain meta-data in the form
 *              of comment lines. One line gives information about the
 *              contracted propagators in the form
 *
 *                S0_th0_m0_r0_ll ^ \dag S0_th0_m0_r0_ll
 *
 *              where we are (currently) interested in the indices
 *              on the m and r tokens. Another comment line specifies
 *              the spin content of the correlator in the form
 *              "P5P5" or "V1A0", for example. This function
 *              reads through the whole file and extracts the two
 *              m and two r indices as well as the spin combination.
 *              It constructs a key based on this information and stores
 *              the file position of the character after the newline behind
 *              the spin combination tag. Thus it stores the position of the
 *              beginning of the correlator corresponding to the constructed key.
 *
 * @param ifs file to be parsed
 * @param filemap output, map of file positions for the various correlators in the file
 * @param corrtype: what kind of correlator to look for, either 1, with specified m, r, or 2, with any possible string containing an underscore specified
 */
inline void map_file(std::ifstream &ifs, std::map<std::string, std::iostream::pos_type> & filemap, const size_t corrtype=1){
  if( !ifs.good() ){
    stop("map_file: input file stream not in a good state!");
  }

  std::string linebuf;
  std::string::size_type comment_pos;

  if(corrtype==1){
      //~ std::cout << "corrtype 2pt" << std::endl;
    // currently the key for a two-point function consists of four unsigned integers (two mass and two r indices)
    // and the name of the spin combination stored in the correlator
    // a generalisation of the reader will have to modify this as well as the key construction
    // for now we store the numeric key components in this array
    unsigned int key_components[4];
    while( ifs.good() ){
      std::getline(ifs, linebuf);
      // we search for commented lines
      comment_pos = linebuf.find("#");
      if( comment_pos != std::string::npos ){
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
          while( key_components_counter != 4 || token != NULL ){
            if( token[0] == 'm' || token[0] == 'r'){
              key_components[key_components_counter] = atoi(token+1);
              key_components_counter++;
            }
            token = strtok(NULL, "_");
          }
          if(key_components_counter != 4 && token == NULL){
            char message[200];
            snprintf(message, 200, "map_file: unable to construct key components in parsing of '%s'", linebuf.c_str());
            stop(message);
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
      } // found commented line
    } //while good
  } // if corrtype==1
  if(corrtype==2){
      //~ std::cout << "corrtype newcorr" << std::endl;
    // currently the key for any correlator is the name of the first and the second correlator, separated by a "^ \dag and "
    // and the name of the spin combination stored in the correlator
    std::string key_components[2];
    while( ifs.good() ){
      std::getline(ifs, linebuf);
      // we search for commented lines
      comment_pos = linebuf.find("#");
      if( comment_pos != std::string::npos ){
          //~ std::cout << linebuf.c_str() << std::endl;
        // in these commented lines, we extract either
        // the current operators
        // or the current spin combination
        // the line looks like so:
        // " # Contraction of C_H ^ \dag and Dth0_A0_C_P_H_H_S_H "
        std::string::size_type contr_pos = linebuf.find("Contraction");
        if( contr_pos != std::string::npos ){
          unsigned int key_components_counter = 0;
          std::string token;
          std::string delimiter (" \0");
          size_t pos = linebuf.find_first_of(delimiter, 0);// = linebuf.find(' ');
          // The different parts of the string are separated by spaces
          // Take the part of the string until the first space.
          // See if it contains an underscore. 
          // If it does, this is the name of an operator
          // delete first part of line and restart search
        // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
          size_t last_part=0, j = 0;
          while( key_components_counter != 2 &&  last_part < 2 && j < 50) {
            token = linebuf.substr(0, pos);
            //~ std::cout << key_components_counter << " " << pos << " " << token << std::endl;
            if( token.find("_") != std::string::npos){
              key_components[key_components_counter] = (std::string)token;
              key_components_counter++;
            }
            linebuf.erase(0, pos+1);
            j++;
            pos = linebuf.find_first_of(delimiter);
            if (pos == std::string::npos) {
              last_part++;
              pos = linebuf.length();
            }
          } 
          if(key_components_counter != 2){
            char message[200];
            snprintf(message, 200, "map_file: unable to construct key components in parsing of '%s'", linebuf.c_str());
            stop(message);
          }
        } else {
          // the line looks like so:
          // " # P5S0"
          std::string::size_type last_space_pos = linebuf.find_last_of(" ");
          std::string spin_comb = linebuf.substr(last_space_pos+1);
          // now we can build the key
          std::string key = make_key_newcorr(key_components[0], key_components[1], spin_comb);
          
        //~ std::cout << key << std::endl;
          // and store the position after the current newline as the starting point
          // of the present correlator
          filemap[ key ] = ifs.tellg();
        }
      } // found commented line
    } //while good
  } // if corrtype==2pt
}

/**
 * @brief driver for the reading of Nissa text correlation functions
 * @description see R documentation of "readnissatextcf" for details on
 *              the parameters passed here
 */
// [[Rcpp::export]]
NumericMatrix read_nissa_textcf_kernel(
    CharacterVector file_basenames_to_read,
    CharacterVector smear_combs_to_read,
    const unsigned int nts,
    DataFrame combs_to_read, 
    size_t corrtype = 1)
{
//~ std::cout << "Hello\n";
  if (corrtype != 1 && corrtype != 2){
        char message[200];
        snprintf(message, 200, "No valid corrtype was given! You gave %lu, valid corrtypes are 1 for 2pt and 2 for newcorr", corrtype);
        stop(message);
  }

  const unsigned int n_correls = combs_to_read.nrows();
  const unsigned int n_smear_combs = smear_combs_to_read.size();
  const unsigned int n_files = file_basenames_to_read.size();

  // prepare some memory for output
  NumericMatrix cf_data(n_files, 2*nts*n_correls*n_smear_combs );

  CharacterVector spin_combs = combs_to_read["spin_comb"];
  CharacterVector op1_idcs, op2_idcs;
  IntegerVector r1_idcs, r2_idcs, m1_idcs, m2_idcs;
  // different correlator naming schemes require different input: type 1 needs info on all r, m, type 2 needs names, both need spin comb
  if (corrtype == 1) {
    r1_idcs = combs_to_read["r1_idx"];
    r2_idcs = combs_to_read["r2_idx"];
    m1_idcs = combs_to_read["m1_idx"];
    m2_idcs = combs_to_read["m2_idx"];
  }
  if (corrtype == 2){
    op1_idcs = combs_to_read["op1_idx"];
    op2_idcs = combs_to_read["op2_idx"];
  }

  std::vector<double> realbuf(nts);
  std::vector<double> imagbuf(nts);

  for(unsigned int ifile = 0; ifile < n_files; ++ifile ){
    for(unsigned int ismear_comb = 0; ismear_comb < n_smear_combs; ++ismear_comb ){

      std::string filename;
      
      //for type 2: use smeaar_comb to indicate file ending
      if(corrtype == 1) filename = (std::string)(file_basenames_to_read[ifile]) + std::string("_") + 
                             (std::string)(smear_combs_to_read[ismear_comb]);
      if(corrtype == 2) filename = (std::string)(file_basenames_to_read[ifile]) + std::string("/") + 
                             (std::string)(smear_combs_to_read[ismear_comb]);
      //~ std::cout << filename << std::endl;

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
      map_file(ifs, filemap, corrtype=corrtype);
      // bring file back into good state
      ifs.clear();

      for( unsigned int icorrel = 0; icorrel < n_correls; ++icorrel ){
        // observable runs slowest, then smearing combination, then time, then real / imag
        // that should be the fastest way to split it into real and imaginary parts
        unsigned int out_idx = (icorrel * n_smear_combs + ismear_comb) * nts * 2;

        // if we ever support other correlation functions from Nissa,
        // a simple point to generalise would be here in the way the key is
        // constructed
        std::string key;
        if(corrtype==1){ 
            key = make_key_2pt((unsigned int)m1_idcs[icorrel], 
                               (unsigned int)m2_idcs[icorrel], 
                               (unsigned int)r1_idcs[icorrel], 
                               (unsigned int)r2_idcs[icorrel],
                               (std::string)spin_combs[icorrel]);
        }
        if(corrtype==2){
            key = make_key_newcorr((std::string)op1_idcs[icorrel],
                                   (std::string)op2_idcs[icorrel],
                                   (std::string)spin_combs[icorrel]);
        }
        if( filemap.count(key) != 1 ){
          char message[200];
          snprintf(message, 200, "correlator %s does not exist in file %s!", key.c_str(), filename.c_str());
          stop(message);
        }
        ifs.seekg( filemap[key] );
        // read the single correlator at the key position
        read_correl(ifs, nts, realbuf, imagbuf);
        // and copy it to the output
        for(unsigned int t=0; t<nts; ++t){
          cf_data(ifile,out_idx  ) = realbuf[t];
          cf_data(ifile,out_idx+1) = imagbuf[t];
          out_idx += 2;
        }
      } // icorrel
    } // ismear_comb
  } // ifile
 return cf_data;
}
