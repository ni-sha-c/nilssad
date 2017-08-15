#include <stdio.h>
#include <stdlib.h>

void driver_fortran_primal_(double *, double *, int *, double *);

void main()
{

    FILE *fptr;
	double s;
	double y[3], X0[3];	
	char ch;
	int steps = 10;

    fptr = fopen("u0", "r");

    if (fptr == NULL)

    {

        printf("Cannot open file \n");

        exit(0);

    }

    ch = fgetc(fptr);
	s = 28.e0;
	y[0] = y[1] = y[2] = 0.1;
	X0[0] = X0[1] = 1.e0;
	X0[2] = 28.e0;
	driver_fortran_primal_(X0,&s,&steps,y);
	printf("y is: %f %f %f ", y[0], y[1], y[2]);
    while (ch != EOF)

    {

        //printf ("%c", ch);

        ch = fgetc(fptr);

    }

    fclose(fptr);

}
