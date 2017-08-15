#include <stdio.h>
#include <stdlib.h>

void driver_fortran_ht_(int *, double *, double *, double *, double *, int *);

void main()
{

    FILE *fptr;
	double ds, s;
	double y[3], X0[3], V0[3];	
	char ch;
	int i = 1, j = 10;

    /*  open the file for reading */

    fptr = fopen("u0", "r");

    if (fptr == NULL)

    {

        printf("Cannot open file \n");

        exit(0);

    }

    ch = fgetc(fptr);
	s = 5.e0;
	y[0] = y[1] = y[2] = 0.1;
	V0[0] = V0[1] = V0[2] = 0.1;
	y[0] = y[1] = y[2] = 0.1;
	driver_fortran_ht_(&i,&s,X0,V0,y,&j);
	printf("y is: %f %f %f ", y[0], y[1], y[2]);
    while (ch != EOF)

    {

        printf ("%c", ch);

        ch = fgetc(fptr);

    }

    fclose(fptr);

}
