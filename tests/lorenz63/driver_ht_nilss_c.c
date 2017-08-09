#include <stdio.h>
#include <stdlib.h>

 

void main()

{

    FILE *fptr;
	double ds, s;
	double y[3], X0[3], V0[3];	
	char ch;

    /*  open the file for reading */

    fptr = fopen("u0", "r");

    if (fptr == NULL)

    {

        printf("Cannot open file \n");

        exit(0);

    }

    ch = fgetc(fptr);

	our_rev_mode%tape= true;
	head_homogeneous_(ds,s,y,X0,V0,10);
    while (ch != EOF)

    {

        printf ("%c", ch);

        ch = fgetc(fptr);

    }

    fclose(fptr);

}
