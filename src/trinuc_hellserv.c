#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/*
This is the program to calculate the BENDING and CURVATURE
according the algorithm of Goodsell and Dickerson (NAR, 1994)
on the basis of known shift vector, and roll, tilt and twist matrixes

The program should read the nucleotide sequence (in FASTA format)
The program should also read three matrices
 */

double pi = 3.14159265358979323;


static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

double cover = 15.0, cunder = 3.0;
double binsize = 1.5;
double hist[32];

/* DEFAULTS */

/* rbend is a bracket across which the bend is calculated */
/* rcurve is a bracket for calculating curvature; The curvature
at base pair N is the angle between averaged normal vectors at base pairs
N-rcurve and N+rcurve.
 */



/*                        A      C      G      T             */


double mat_roll[4][4], mat_twist[4][4], mat_tilt[4][4];

/* DNaseI trinuc scale */

double DNaseI_bend_roll[4][4][4] = {
    {
        {0.1, 1.6, 4.2, 0.0},
        {5.8, 5.2, 5.2, 2.0},
        {6.5, 6.3, 4.7, 2.0},
        {9.7, 3.6, 8.7, 0.0}
    },
    {
        {6.2, 6.8, 9.6, 8.7},
        {0.7, 5.7, 3.0, 4.7},
        {5.8, 4.3, 3.0, 5.2},
        {7.8, 6.6, 9.6, 4.2}
    },
    {
        {5.1, 5.6, 6.6, 3.6},
        {7.5, 8.2, 4.3, 6.3},
        {6.2, 8.2, 5.7, 5.2},
        {6.4, 5.6, 6.8, 1.6}
    },
    {
        {7.3, 6.4, 7.8, 9.7},
        {10.0, 6.2, 5.8, 6.5},
        {10.0, 7.5, 0.7, 5.8},
        {7.3, 5.1, 6.2, 0.1}
    }
};



double Consensus_bend_roll[4][4][4] = {
    {
        {0.0633, 2.64115, 4.69915, 0.35},
        {5.4903, 5.31605, 5.3055, 3.9232},
        {4.8884, 6.8829, 5.0523, 3.9232},
        {6.2734, 4.44325, 7.7171, 0.35}
    },
    {
        {4.7618, 6.62555, 6.8996, 7.7171},
        {3.05865, 5.827, 3.869, 5.0523},
        {7.07195, 5.89135, 3.869, 5.3055},
        {5.00295, 5.9806, 6.8996, 4.69915}
    },
    {
        {4.0633, 5.51645, 5.9806, 4.44325},
        {6.75525, 9.0823, 5.89135, 6.8829},
        {4.9907, 9.0823, 5.827, 5.31605},
        {5.0673, 5.51645, 6.62555, 2.64115}
    },
    {
        {4.6709, 5.0673, 5.00295, 6.2734},
        {7.7, 4.9907, 7.07195, 4.8884},
        {7.7, 6.75525, 3.05865, 5.4903},
        {4.6709, 4.0633, 4.7618, 0.0633}
    }
};


double Nucleosome_bend_roll[4][4][4] = {
    {
        {0.0, 3.7, 5.2, 0.7},
        {5.2, 5.4, 5.4, 5.8},
        {3.3, 7.5, 5.4, 5.8},
        {2.8, 5.3, 6.7, 0.7}
    },
    {
        {3.3, 6.5, 4.2, 6.7},
        {5.4, 6.0, 4.7, 5.4},
        {8.3, 7.5, 4.7, 5.4},
        {2.2, 5.4, 4.2, 5.2}
    },
    {
        {3.0, 5.4, 5.4, 5.3},
        {6.0, 10., 7.5, 7.5},
        {3.8, 10., 6.0, 5.4},
        {3.7, 5.4, 6.5, 3.7}
    },
    {
        {2.0, 3.7, 2.2, 2.8},
        {5.4, 3.8, 8.3, 3.3},
        {5.4, 6.0, 5.4, 5.2},
        {2.0, 3.0, 3.3, 0.0}
    }
};

void swap(double to[4][4][4], double from[4][4][4]) {
    int i, j, k;

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {
                to[i][j][k] = from[i][j][k];
            }
        }
    }

}

void err(char *txt) {
    printf("Bad luck:\n%s\n", txt);
    exit(1);
}

void readmat(double mat[4][4], char name[50]) {

    FILE *fil;
    double a;
    int i, j;

    if ((fil = fopen(name, "r")) == NULL) {
        /*complain */
        printf("\nCheck the matrix filename.... %s\n", name);
        exit(1);
    }
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            fscanf(fil, "%lf", &a);
            mat[i][j] = a;
            printf("%6.2f", a);
            fclose(fil);
        }
        printf("\n");
    }

}

void deg2rad(double mat[4][4]) {

    int i, j;

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            mat[i][j] *= pi;
            mat[i][j] /= 180;
        }
    }
}

void load_seq(char *seq, int *seq_num) {

    int i = 0;

    while (seq[i++]) {
        if ((seq[i] == 'A') || (seq[i] == 'a')) seq_num[i] = 0;
        if ((seq[i] == 'C') || (seq[i] == 'c')) seq_num[i] = 1;
        if ((seq[i] == 'G') || (seq[i] == 'g')) seq_num[i] = 2;
        if ((seq[i] == 'T') || (seq[i] == 't')) seq_num[i] = 3;
    }
}

void curve_it(double *curve, double *bend, double *gc, char *seq, int rcurve, int rbend, char *matrix) {

    int p, v;
    int i, j;

    double a, b, c;
    double twistsum, curvescale, bendscale;
    double angle;
    double dx, dy, rxsum, rysum;
    double tmp;

    double *xave, *yave;
    double *x, *y;

    int *seq_num;

    long seqlen;

    char roll[50], tilt[50], twist[50];

    FILE *f;

    bendscale = 2.1;
    curvescale = 0.7;

    /* 
    these are proportionality coeffs.; bends are scaled by three relatives
    to the curves
     */


    /* default DNA parameters */
    /* geometry */
    a = 5.0;
    b = 6.0;
    c = 3.3;
    angle = acos(a / b);


    seqlen = strlen(seq);



    /* strcpy(roll, "/usr7/www/cgi-bin/dna/"); */
    strcpy(roll, matrix);
    strcat(roll, "_roll.mat");

    /* strcpy(tilt, "/usr7/www/cgi-bin/dna/"); */
    strcpy(tilt, matrix);
    strcat(tilt, "_tilt.mat");

    /* strcpy(twist, "/usr7/www/cgi-bin/dna/"); */
    strcpy(twist, matrix);
    strcat(twist, "_twist.mat");

    if ((f = fopen(roll, "r")) == NULL) {
        /*complain */
        printf("\nCheck the roll matrix filename.... %s\n", roll);
        exit(1);
    }
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            fscanf(f, "%lf", &a);
            mat_roll[i][j] = a * pi / 180.0;
        }
    }
    fclose(f);

    if ((f = fopen(tilt, "r")) == NULL) {
        /*complain */
        printf("\nCheck the tilt matrix filename....\n");
        exit(1);
    }

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            fscanf(f, "%lf", &a);
            mat_tilt[i][j] = a * pi / 180.0;
        }
    }
    fclose(f);

    if ((f = fopen(twist, "r")) == NULL) {
        /*complain */
        printf("\nCheck the twist matrix filename....\n");
        exit(1);
    }

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            fscanf(f, "%lf", &a);
            mat_twist[i][j] = a * pi / 180.0;
        }
    }

    /*
    deg2rad(mat_roll);
    deg2rad(mat_twist);
    deg2rad(mat_tilt);
     */




    /* program */



    xave = (double *) malloc(sizeof (double) * seqlen);
    if (!xave) err("XAVE allocation failure");
    yave = (double *) malloc(sizeof (double) * seqlen);
    if (!yave) err("YAVE allocation failure");
    x = (double *) malloc(sizeof (double) * seqlen);
    if (!x) err("X allocation failure");
    y = (double *) malloc(sizeof (double) * seqlen);
    if (!y) err("Y allocation failure");
    seq_num = (int *) malloc(sizeof (int) * seqlen);

    load_seq(seq, seq_num);

    /* Calculate trajectory of helix axis */

    twistsum = 0.0;
    x[0] = 0.0;
    y[0] = 0.0;
    for (i = 0; i < seqlen - 1; i++) {
        p = seq_num[i];
        v = seq_num[i + 1];
        twistsum += mat_twist[p][v];
        dx = mat_roll[p][v] * sin(twistsum) + mat_tilt[p][v] * sin(twistsum - pi / 2.0);
        dy = mat_roll[p][v] * cos(twistsum) + mat_tilt[p][v] * cos(twistsum - pi / 2.0);
        x[i + 1] = x[i] + dx;
        y[i + 1] = y[i] + dy;
    }

    /* Calculate average of trajectory over 10 base pairs */
    for (i = 5; i < seqlen - 5; i++) {
        rxsum = 0.0;
        rysum = 0.0;
        for (j = -4; j < 5; j++) {
            rxsum += x[i + j];
            rysum += y[i + j];
        }
        rxsum += x[i + 5] / 2.0 + x[i - 5] / 2.0;
        rysum += y[i + 5] / 2.0 + y[i - 5] / 2.0;
        xave[i] = rxsum / 10.0;
        yave[i] = rysum / 10.0;
    }

    /* now we should calculate GC-content */
    /* it will be calculated on the window [N-rcurve, N+rcurve] */


    /* printf ("Calculating GC-content...\n"); */


    for (i = rcurve; i < seqlen - rcurve + 1; i++) {
        gc[i] = 0;

        /*
        printf ("Window %d\n",i);
         */
        for (j = -rcurve; j < rcurve + 1; j++) {
            /*
            printf ("%d\t",seq_num[i+j]);
             */
            if ((seq_num[i + j] == 1) || (seq_num[i + j] == 2)) {
                /*
                printf ("!\n");
                 */
                gc[i] += 1;
            }

        }


    }

    for (i = rcurve; i < seqlen - rcurve + 1; i++) gc[i] /= (rcurve * 2.0 / 100);

    /* Calculate bend profile */

    for (i = rbend; i < seqlen - rbend; i++) {
        p = i + rbend;
        v = i - rbend;
        bend[i] = sqrt((x[p] - x[v])*(x[p] - x[v]) + (y[p] - y[v]) * (y[p] - y[v]));
    }

    /* calculate curve profile */

    for (i = rcurve + 5; i < seqlen - rcurve - 5; i++) {
        p = i + rcurve;
        v = i - rcurve;
        curve[i] = sqrt((xave[p] - xave[v])*(xave[p] - xave[v])+(yave[p] - yave[v])*(yave[p] - yave[v]));


    }




    for (i = rcurve + 5; i < seqlen - rcurve - 5; i++) {
        bend[i] *= ((10.5 / (rcurve * 2 + 1))*(180 / pi));
        curve[i] *= ((10.5 / (rcurve * 2 + 1))*(180 / pi));


    }

    free(seq_num);
    free(y);
    free(x);
    free(yave);
    free(xave);
}

void trinuc_curve_it(double *curve, double *bend, double *gc, char *seq, int rcurve, int rbend, char *matrix) {

    int p, v, q;
    int i, j;

    double a, b, c;
    double twistsum, curvescale, bendscale;
    double angle;
    double dx, dy, rxsum, rysum;
    double tmp;

    double *xave, *yave;
    double *x, *y;
    double bend_roll[4][4][4];


    int *seq_num;

    long seqlen;


    seqlen = strlen(seq);

    xave = (double *) malloc(sizeof (double) * seqlen);
    if (!xave) err("XAVE allocation failure");
    yave = (double *) malloc(sizeof (double) * seqlen);
    if (!yave) err("YAVE allocation failure");
    x = (double *) malloc(sizeof (double) * seqlen);
    if (!x) err("X allocation failure");
    y = (double *) malloc(sizeof (double) * seqlen);
    if (!y) err("Y allocation failure");
    seq_num = (int *) malloc(sizeof (int) * seqlen);

    load_seq(seq, seq_num);

    /* DEFAULTS */

    /* 
    these are proportionality coeffs.; bends are scaled by three relatives
    to the curves
     */

    bendscale = 2.1;
    curvescale = 0.7;


    /* default DNA parameters */
    /* geometry */
    a = 5.0;
    b = 6.0;
    c = 3.3;
    angle = acos(a / b);



    /* which scale are we using? */

    switch (matrix[0]) {
        case 'D': swap(bend_roll, DNaseI_bend_roll);
            break;
        case 'C': swap(bend_roll, Consensus_bend_roll);
            break;
        case 'N': swap(bend_roll, Nucleosome_bend_roll);
            break;
        default: swap(bend_roll, DNaseI_bend_roll);
    }



    /* Calculate trajectory of helix axis */

    twistsum = 0.0;
    x[0] = 0.0;
    y[0] = 0.0;
    for (i = 0; i < seqlen - 2; i++) {
        p = seq_num[i];
        v = seq_num[i + 1];
        q = seq_num[i + 2];

        /* 36 deg thingie
        twistsum += 0.62831853 ;
         */

        twistsum += 0.5983986;

        dx = bend_roll[p][v][q]*.017453292 * sin(twistsum);
        dy = bend_roll[p][v][q]*.017453292 * cos(twistsum);
        x[i + 1] = x[i] + dx;
        y[i + 1] = y[i] + dy;
    }

    /* Calculate average of trajectory over 10 base pairs */
    for (i = 5; i < seqlen - 5; i++) {
        rxsum = 0.0;
        rysum = 0.0;
        for (j = -4; j < 5; j++) {
            rxsum += x[i + j];
            rysum += y[i + j];
        }
        rxsum += x[i + 5] / 2.0 + x[i - 5] / 2.0;
        rysum += y[i + 5] / 2.0 + y[i - 5] / 2.0;
        xave[i] = rxsum / 10.0;
        yave[i] = rysum / 10.0;
    }

    /* now we should calculate GC-content */
    /* it will be calculated on the window [N-rcurve, N+rcurve] */


    for (i = rcurve; i < seqlen - (rcurve + 1); i++) {
        gc[i] = 0;
        for (j = -rcurve; j < rcurve + 1; j++) {
            if ((seq_num[i + j] == 1) || (seq_num[i + j] == 2)) {
                gc[i] += 1;
            }
        }
        gc[i] /= (rcurve * 2);
    }


    /* Calculate bend profile */

    for (i = rbend; i < seqlen - (rbend + 2); i++) {
        tmp = 0;
        bend[i] = 0;

        if (rbend) {
            for (j = -rbend; j < rbend; j++) {

                p = seq_num[i + j];
                v = seq_num[i + j + 1];
                q = seq_num[i + j + 2];

                tmp += bend_roll[p][v][q];
            }
            bend[i] = tmp / (rbend * 2);
        } else {
            p = seq_num[i];
            v = seq_num[i + 1];
            q = seq_num[i + 2];
            bend[i] = bend_roll[p][v][q];
        }
    }

    /* calculate curve profile */

    for (i = rcurve + 5; i < seqlen - rcurve - 5; i++) {
        p = i + rcurve;
        v = i - rcurve;
        curve[i] = sqrt((xave[p] - xave[v])*(xave[p] - xave[v])+(yave[p] - yave[v])*(yave[p] - yave[v]));
    }


    /*
    we want bend, and curvature profile now in deg/helical turn (10.5)
     */


    for (i = rcurve + 5; i < seqlen - rcurve - 5; i++) {
        curve[i] *= ((10.5 / (rcurve * 2 + 1))*(180 / pi));
        /*
        printf("%6.3f\t%6.3f\n",(double)gc[i]/rcurve/2.0, curve[i]); 
         */
    }


    /* Free what thou allocated have! */

    free(seq_num);
    free(y);
    free(x);
    free(yave);
    free(xave);
}





