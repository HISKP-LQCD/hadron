#include<Rcpp.h>
#include<numeric>
#include <fstream>
#include <iostream>
#include <string.h>

using namespace Rcpp;
using namespace std;

//Todo: Comments
//      Warning for not consistent order
inline std::vector<double> readFile(std::istream &is,double Nt)
{
  std::vector<double> result;
  for(int i=0;i<Nt;++i)
  {
    double number_real;
    double number_imag;
    is>>number_real;
    is>>number_imag;
    result.push_back(number_real);
    result.push_back(number_imag);
  }
  return result;
}

// [[Rcpp::export]]
NumericVector reading_nissa_corr(CharacterVector spinstructures,int confnumber,int Nt,string dirname,int nmasses,NumericVector nmasses1,NumericVector nmasses2,NumericVector r1,NumericVector r2){
   std::filebuf fb;
   int tindex;
   int confindex;
   int r1idx,r2idx;
   int N_muindex1,N_muindex2;
   char filename[100];
   CharacterVector gamma_names(61);
   NumericVector gammastructures(61);
   int sizem1=nmasses1.size();
   int sizem2=nmasses2.size();
   int sizer1=r1.size();
   int sizer2=r2.size();
   if ((sizem1 != sizem2) || (sizem1 != sizer1) || (sizem1 != sizer2)){
     Rprintf("Error the size of m1,m2,r2,r1 arrays has to be the same");
     exit(1);
   }
   for (int i=0; i<61; i++) { gammastructures[i] = i; };
   gamma_names[0]="P5S0";
   gamma_names[1]="P5V1";
   gamma_names[2]="P5V2";
   gamma_names[3]="P5V3";
   gamma_names[4]="P5V0";
   gamma_names[5]="P5P5";
   gamma_names[6]="P5A1";
   gamma_names[7]="P5A2";
   gamma_names[8]="P5A3";
   gamma_names[9]="P5A0";
   gamma_names[10]="P5T1";
   gamma_names[11]="P5T2";
   gamma_names[12]="P5T3";
   gamma_names[13]="P5B1";
   gamma_names[14]="P5B2";
   gamma_names[15]="P5B3";
   gamma_names[16]="S0P5";
   gamma_names[17]="V1P5";
   gamma_names[18]="V2P5";
   gamma_names[19]="V3P5";
   gamma_names[20]="V0P5";
   gamma_names[21]="P5P5";
   gamma_names[22]="A1P5";
   gamma_names[23]="A2P5";
   gamma_names[24]="A3P5";
   gamma_names[25]="A0P5";
   gamma_names[26]="T1P5";
   gamma_names[27]="T2P5";
   gamma_names[28]="T3P5";
   gamma_names[29]="B1P5";
   gamma_names[30]="B2P5";
   gamma_names[31]="B3P5";
   gamma_names[32]="S0S0";
   gamma_names[33]="V0V0";
   gamma_names[34]="A0A0";
   gamma_names[35]="V1V1";
   gamma_names[36]="V2V2";
   gamma_names[37]="V3V3";
   gamma_names[38]="A1A1";
   gamma_names[39]="A2A2";
   gamma_names[40]="A3A3";
   gamma_names[41]="T1T1";
   gamma_names[42]="T2T2";
   gamma_names[43]="T3T3";
   gamma_names[44]="V1T1";
   gamma_names[45]="V2T2";
   gamma_names[46]="V3T3";
   gamma_names[47]="T1V1";
   gamma_names[48]="T2V2";
   gamma_names[49]="T3V3";
   gamma_names[50]="B1B1";
   gamma_names[51]="B2B2";
   gamma_names[52]="B3B3";
   gamma_names[53]="A1B1";
   gamma_names[54]="A2B2";
   gamma_names[55]="A3B3";
   gamma_names[56]="B1A1";
   gamma_names[57]="B2A2";
   gamma_names[58]="B3A3";
   gamma_names[59]="V0S0";
   gamma_names[60]="S0V0";
   gammastructures.names()=gamma_names;
   vector<double> mat(2*Nt);
   NumericMatrix matR(confnumber,Nt*spinstructures.size()*nmasses1.size()*8);
   CharacterVector localsmeared(4);
   localsmeared[0]="ll";
   localsmeared[1]="ls";
   localsmeared[2]="sl";
   localsmeared[3]="ss";
   int summ=0;
   int indextable[nmasses][nmasses][2][2];
   for ( int mu2 =0; mu2 < nmasses ; mu2++) // Nmu is the number of mu values
     for ( int r1i = 0; r1i <=1; r1i++) // Wilson p a r a m e t e r r = -1 ,1
       for ( int mu1 =0; mu1 <= mu2 ; mu1++) // !!!! careful to the ending of the loop
         for ( int r2i = 0; r2i <=1; r2i++){
           indextable[mu2][mu1][r1i][r2i]=summ;
           summ++;
         }

   for (confindex=1;confindex<=confnumber;++confindex){
     int store_index=0;
     for(CharacterVector::iterator smearing = localsmeared.begin(); smearing != localsmeared.end(); ++smearing){
       string temporary=(string)(*smearing);
       sprintf(filename,"%s/%04d/mes_contr_2pts_%s",dirname.c_str(),confindex,temporary.c_str());
       for(CharacterVector::iterator spin = spinstructures.begin(); spin != spinstructures.end(); ++spin) {
         fb.open(filename,std::ios::in);
         std::istream is(&fb);
         int line_number=0;
         int spinindex=gammastructures[(string)*spin];
         string line;
         for (int mindex=0; mindex<nmasses1.size(); ++mindex){
           if (nmasses1[mindex] > nmasses2[mindex]){
             Rprintf("Order of masses is not correct\n");
             exit(1);
           }
           int index=indextable[(int)(nmasses2[mindex])][(int)(nmasses1[mindex])][(int)(r1[mindex])][(int)(r2[mindex])];
           int starting_point=(4+(Nt+2)*61)*index;
           for (;line_number<starting_point;++line_number)
             std::getline(is, line); 
           int endpoint=4;
           for(int j=0;j<endpoint;++j){
             getline(is,line);
             line_number++;
             if (confindex==1){
              cout<<line<<endl;
             }
           }
           if (confindex==1){
             Rprintf("Correlator to be readed\n");
           }
           for(int j=0;j<spinindex*(Nt+2);++j){
             line_number++;
             getline(is,line);
           }
           for(int j=0;j<1;++j){
             getline(is,line);
             line_number++;
             if (confindex==1){
               cout<<line<<endl;
             }
           }
           mat = readFile(is,Nt);
           line_number+=Nt;
           getline(is,line);
           line_number++;
           getline(is,line);
           for(int j=0;j<2*Nt;++j,++store_index){
             matR(confindex-1,store_index)=mat[j];
           }
         }
         fb.close();
       }
     }
   }
   return matR;
}
