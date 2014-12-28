#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static int pmodes_initialised = 0;
static int * degnrtDOF = NULL;
static int * arrayPmode = NULL;

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
      int	num = 0;
      
      for(int x=-r; x<=r; x++)
	for(int y=-r; y<=r; y++)
	  for(int z=-r; z<=r; z++){
	    
	    if(x*x+y*y+z*z == pmodeSqur){
	      arrayPmode[pmodeSqur*_DimMAX*3 + num*3 + 0] = x;
	      arrayPmode[pmodeSqur*_DimMAX*3 + num*3 + 1] = y;
	      arrayPmode[pmodeSqur*_DimMAX*3 + num*3 + 2] = z;
	      
	      num = num + 1;
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




