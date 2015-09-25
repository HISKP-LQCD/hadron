#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static int pmodes_initialised = 0;
static int * degnrtDOF = NULL;
static int * arrayPmode = NULL;

//NPmode should be bigger than the biggest pmodeSqur we may meet in the iteration,
//DimMax is the biggest degenation degree within the range of NPmode
//Usually we should set big enough value for these two parameters or the program will crash.
//For precision of 1e-8, (NPmode=40, DimMAX=72) is good enough.
//Actually,by my test NPmode=40 can be selected up to precision 1e-11 for usual case.
//
//Other selections may be 
//(NPmode=70,DimMAX=96)
//NPmode=70 can be selected up to precision 1e-16 for usual case.
//
//(NPmode=100,DimMAX=120),
//(NPmode=145,DimMAX=168)


void get_npmode(int * pair, const int i) {

  if(i < 2) {
    pair[0] = 70;
    pair[1] = 96;
  }
  else if(i == 2) {
    pair[0] = 100;
    pair[1] = 120;
  }
  else if(i == 3){
    pair[0] = 145;
    pair[1] = 168;
  }
  else {
    pair[0] = 200;
    pair[1] = 240;    
  }
  return;
}

void pmode_free_arrays() {
  if(pmodes_initialised) {
    free(degnrtDOF);
    free(arrayPmode);
  }
  pmodes_initialised = 0;
  degnrtDOF = NULL;
  arrayPmode = NULL;
  return;
}


int gen_points_array(int ** _degnrtDOF, int ** _arrayPmode, const int NPmode, const int DimMAX) {	
  static int _NPmode = 0;
  static int _DimMAX = 0;

  if(!pmodes_initialised || NPmode > _NPmode || DimMAX >  _DimMAX) {
    pmode_free_arrays();

    if(NULL == (degnrtDOF = (int *)malloc(NPmode*sizeof(int)))) {
      printf("Malloc wrong for degnrtDOF!\n");
      return(-1);
    }
    
    if(NULL == (arrayPmode= (int *)malloc(NPmode*DimMAX*3*sizeof(int)))) {
      printf("Malloc wrong for arrayPmode!\n");
      return(-2);
    }
    pmodes_initialised = 1;
    _NPmode = NPmode;
    _DimMAX = DimMAX;

    for(int pmodeSqur = 0; pmodeSqur < _NPmode; pmodeSqur++) {
      int r = (int)(floor(sqrt(pmodeSqur)));
      int num = 0;
      
      for(int x = -r; x <= r; x++)
	for(int y = -r; y <= r; y++)
	  for(int z = -r; z <= r; z++){
	    
	    if(x*x+y*y+z*z == pmodeSqur){
	      arrayPmode[pmodeSqur*_DimMAX*3 + num*3 + 0] = x;
	      arrayPmode[pmodeSqur*_DimMAX*3 + num*3 + 1] = y;
	      arrayPmode[pmodeSqur*_DimMAX*3 + num*3 + 2] = z;
	      
	      num++;
	    }
	  }
      
      //degnrtDOF also record those pmodeSqur without any
      //corresponding points, in which case num = 0;
      degnrtDOF[pmodeSqur] = num;
    }
  }
  *_degnrtDOF = degnrtDOF;
  *_arrayPmode = arrayPmode;

  return 0;
}

/*
void Array_Alloc(int * arrayPoints, int npmode, int dimmax, int innerdim)
{
	printf("The volume of the array is %d.\n", npmode*dimmax*innerdim);
	if((arrayPoints = (int *)malloc(npmode*dimmax*innerdim*sizeof(int))) == NULL)
		{
			printf("Error: void Array_alloc(...\n");
			exit(-1);
		}
}
*/




