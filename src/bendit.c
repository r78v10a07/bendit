#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/wait.h>

char *program_name;

void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: %s [options]\n", program_name);
    fprintf(stream, "\nTo run the program gnuplot have to be installed\n");
    fprintf(stream, "\nThe program will generate two output files: a text file used to run gnuplot and the resulting image\n");
    fprintf(stream, "\n\n%s options:\n\n", program_name);
    fprintf(stream, "-h,   --help                          Display this usage information.\n");
    fprintf(stream, "-s,   --seq                           File in FASTA with the sequence\n");
    fprintf(stream, "-o,   --output                        Output file (Text file used to create the gnuplot script)\n");
    fprintf(stream, "-S,   --Scale                         Scale to be used (values: Consensus OR DNaseI OR Nucleosome) Default: Consensus\n");
    fprintf(stream, "-c,   --cwindows                      Curvature window size (Default: 31)\n");
    fprintf(stream, "-b,   --bwindows                      Bendability/G+C content/Complexity window size (Default: 31)\n");
    fprintf(stream, "-t,   --type                          Plot type (values: 1D OR 2D) Default: 1D\n");
    fprintf(stream, "-1,   --smooth1                       Smooth1 (values: none, simple, cubic, bezier) Default: none\n");
    fprintf(stream, "-2,   --smooth2                       Smooth2 (values: none, simple, cubic, bezier) Default: none\n");
    fprintf(stream, "-g,   --complexity                    Complexity  (values: G (G+C content), B (Bendability), C (Complexity)) Default: G\n");
    fprintf(stream, "      --xmin                          X min\n");
    fprintf(stream, "      --xmax                          X max\n");
    fprintf(stream, "      --ymin                          Y min\n");
    fprintf(stream, "      --ymax                          Y max\n");
    fprintf(stream, "      --y2min                         Y2 min\n");
    fprintf(stream, "      --y2max                         Y2 max\n");
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Roberto Vera Alvarez (e-mail: Roberto.Vera@icgeb.org)\n\n");
    fprintf(stream, "********************************************************************************\n");
    exit(0);
}

#define MAX_ENTRIES 100000
#define MAXSEQ 5000000
#define PANIC(a) do { \
		perror(a); \
		if (temp_name) unlink(temp_name);\
		exit(1);\
	} while(0)

typedef struct {
    char *name;
    char *val;
} entry;

char *makeseq(char *seqr, int g1, int g2, int g3);
char *range(char *low, char *high);
char *makeword(char *line, char stop);
char *fmakeword(FILE *f, char stop, int *len);
char x2c(char *what);
void unescape_url(char *url);
void removespace(char *sq);
void plustospace(char *str);
void sendmail(char *email, char *attach, char *body);

/* int curve_it(double *curve, double *bend, double *gc, char *seq, int rcurve, int rbend, char *matrix); */
int trinuc_curve_it(double *curve, double *bend, double *gc, char *seq, int rcurve, int rbend, char *matrix);
int k1_complexity_n(char *seq, int len, double *k1values, int winlen, int word);

int izit(char *seq);

int main(int argc, char *argv[]) {
    int next_option;
    const char* const short_options = "hs:o:S:c:b:t:1:2:g:3:4:5:6:7:8:";

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "seq", 1, NULL, 's'},
        { "output", 1, NULL, 'o'},
        { "scale", 1, NULL, 'S'},
        { "cwindows", 1, NULL, 'c'},
        { "bwindows", 1, NULL, 'c'},
        { "type", 1, NULL, 't'},
        { "smooth1", 1, NULL, '1'},
        { "smooth2", 1, NULL, '2'},
        { "complexity", 1, NULL, 'g'},
        { "xmin", 1, NULL, '3'},
        { "xmax", 1, NULL, '4'},
        { "ymin", 1, NULL, '5'},
        { "ymax", 1, NULL, '6'},
        { "y2min", 1, NULL, '7'},
        { "y2max", 1, NULL, '8'},
        { NULL, 0, NULL, 0} /* Required at end of array.  */
    };

    register int i;
    int whichplot = 1;
    char c;
    char *smth1 = (char *) malloc(sizeof (char) * 50);
    char *smth2 = (char *) malloc(sizeof (char) * 50);
    double *curve, *bend, *gc;
    int seqlen = 0;
    int cwindow = 31;
    int bwindow = 31;
    FILE *data;
    char *seq = NULL;
    char *seq_name = NULL;
    char *scale = NULL;
    char *temp_name;
    char *complexity = (char *) malloc(sizeof (char) * 50);
    char *line = NULL;
    char *iline;
    char *xmin = (char *) malloc(sizeof (char) * 2);
    char *xmax = (char *) malloc(sizeof (char) * 2);
    char *ymin = (char *) malloc(sizeof (char) * 2);
    char *ymax = (char *) malloc(sizeof (char) * 2);
    char *y2min = (char *) malloc(sizeof (char) * 2);
    char *y2max = (char *) malloc(sizeof (char) * 2);
    size_t len = 0;
    ssize_t read;

    program_name = argv[0];

    strcpy(smth1, "  ");
    strcpy(smth2, "  ");
    strcpy(complexity, "G+C content");

    strcpy(xmin, "*");
    strcpy(xmax, "*");
    strcpy(ymin, "*");
    strcpy(ymax, "*");
    strcpy(y2min, "*");
    strcpy(y2max, "*");


    do {
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);

        switch (next_option) {
            case 'h':
                print_usage(stdout, 0);

            case 's':
                seq_name = (char *) malloc(sizeof ( char) * (strlen(optarg) + 1));
                strcpy(seq_name, optarg);
                break;

            case 'o':
                temp_name = (char *) malloc(sizeof ( char) * (strlen(optarg) + 1));
                strcpy(temp_name, optarg);
                break;

            case 'S':
                scale = (char *) malloc(sizeof ( char) * (strlen(optarg) + 1));
                strcpy(scale, optarg);
                break;

            case 'c':
                cwindow = atoi(optarg);
                break;

            case 'b':
                bwindow = atoi(optarg);
                break;

            case 't':
                if (strcmp(optarg, "1D") == 0) {
                    whichplot = 1;
                } else if (strcmp(optarg, "2D") == 0) {
                    whichplot = 2;
                }
                break;

            case '1':
                if (strcmp(optarg, "simple") == 0) {
                    strcpy(smth1, " smooth unique");
                } else if (strcmp(optarg, "cubic") == 0) {
                    strcpy(smth1, " smooth csplines");
                } else if (strcmp(optarg, "bezier") == 0) {
                    strcpy(smth1, " smooth bezier");
                } else {
                    strcpy(smth1, "  ");
                }
                break;

            case 'g':
                if (strcmp(optarg, "B") == 0) {
                    strcpy(complexity, "Bendability");
                } else if (strcmp(optarg, "C") == 0) {
                    strcpy(complexity, "Complexity");
                }
                break;

            case '3':
                xmin = (char *) malloc(sizeof ( char) * (strlen(optarg) + 1));
                strcpy(xmin, optarg);
                break;

            case '4':
                xmax = (char *) malloc(sizeof ( char) * (strlen(optarg) + 1));
                strcpy(xmax, optarg);
                break;

            case '5':
                ymin = (char *) malloc(sizeof ( char) * (strlen(optarg) + 1));
                strcpy(ymin, optarg);
                break;

            case '6':
                ymax = (char *) malloc(sizeof ( char) * (strlen(optarg) + 1));
                strcpy(ymax, optarg);
                break;

            case '7':
                y2min = (char *) malloc(sizeof ( char) * (strlen(optarg) + 1));
                strcpy(y2min, optarg);
                break;

            case '8':
                y2max = (char *) malloc(sizeof ( char) * (strlen(optarg) + 1));
                strcpy(y2max, optarg);
                break;

        }
    } while (next_option != -1);

    if (!seq_name) {
        printf("Give file with the sequence in FASTA format\n");
        print_usage(stdout, 0);
    } else {
        data = fopen(seq_name, "r");
        i = 0;
        while ((read = getline(&line, &len, data)) != -1) {
            iline = line;
            fflush(NULL);
            if (i == 0) {
                if (*line != '>') {
                    printf("Error!!! Give infile with the sequence in FASTA format\n");
                    exit(1);
                } else {
                    if (strlen(line) + 1 > strlen(seq_name)) {
                        seq_name = (char *) realloc(seq_name, sizeof (char) * (strlen(line)));
                    }
                    iline++;
                    while (*iline != '\0' && *iline != '\n' && *iline != '\r' && *iline != EOF && *iline != ' ') {
                        *(seq_name + i++) = *iline++;
                    }
                    *(seq_name + i) = '\0';
                    seq_name = (char *) realloc(seq_name, sizeof (char) * (strlen(seq_name) + 1));
                }
            } else {
                while (*iline != '\0' && *iline != '\n' && *iline != '\r' && *iline != EOF) {
                    if (*iline != ' ' && *iline != '\t') {
                        seq = (char *) realloc(seq, sizeof (char) * (seqlen + 1));
                        *(seq + seqlen++) = *iline;
                    }
                    iline++;
                }
                seq = (char *) realloc(seq, sizeof (char) * (seqlen + 1));
                *(seq + seqlen) = '\0';
            }
        }
        removespace(seq);
        if (line) free(line);
        fclose(data);
    }

    if (!temp_name) {
        printf("Give output file name\n");
        print_usage(stdout, 0);
    } else {
        data = fopen(temp_name, "w");
        if (!data) {
            printf("Can't open output file\n");
            print_usage(stdout, 1);
        }
    }

    if (!scale) {
        scale = strdup("Consensus");
    } else {
        if (strcmp(scale, "Consensus") != 0 && strcmp(scale, "DNaseI") != 0 && strcmp(scale, "Nucleosome") != 0) {
            printf("Scale values are:\n");
            printf("\tConsensus\n");
            printf("\tDNaseI\n");
            printf("\tNucleosome\n");
            print_usage(stdout, 1);
        }
    }

    if (!seqlen) {
        printf("Empty sequence!\n");
        printf("Really! You don't expect to get something without typing the seequence, do you?\n");
        printf("Try being a bit more intuitive!\n");
        exit(1);
    }
    if ((c = izit(seq))) {
        printf("Incorrect sequence!\n");
        printf("Error in character %c\n", c);
        printf("Please check your sequence for characters other than a t g c\
	          (case insensitive)\n");
        printf("Here's what You submitted:\n");
        printf("%s\n", makeseq(seq, 10, 70, i));
        exit(1);
    }

    curve = (double *) malloc(sizeof (double) * seqlen);
    bend = (double *) malloc(sizeof (double) * seqlen);
    gc = (double *) malloc(sizeof (double) * seqlen);

    /* Calculate what you have to */

    if (!strcmp(scale, "DNaseI")) trinuc_curve_it(curve, bend, gc, seq, cwindow / 2, bwindow / 2, scale);
    if (!strcmp(scale, "Consensus")) trinuc_curve_it(curve, bend, gc, seq, cwindow / 2, bwindow / 2, scale);
    if (!strcmp(scale, "Nucleosome")) trinuc_curve_it(curve, bend, gc, seq, cwindow / 2, bwindow / 2, scale);

    if (!strcmp(complexity, "Complexity")) k1_complexity_n(seq, seqlen, bend, cwindow, 2);


    printf("bend.it @ ICGEB Trieste, (c) 2015, R. Vera & S. Pongor\n\n");
    printf("This is a result for your query sequence \"%s\" (%d bp)\n", seq_name, seqlen);
    printf("You requested a tabulated output of %s and Predicted curvature\n", complexity);
    printf("Parameters:\n");
    printf("Profile: %s\n", scale);
    printf("Predicted curvature window size: %d\n", cwindow);
    printf("Bendability window size: %d\n", bwindow);
    printf("===================================================\n");
    printf("Column 1 = Position\n");
    printf("Column 2 = Sequence\n");
    printf("Column 3 = Predicted curvature\n");
    printf("Column 4 = %s\n", complexity);
    printf("---------------------------------------------------\n");
    printf("Ranges:\n");
    line = range(xmin, xmax);
    printf("\tNucleotide: %s\n", line);
    free(line);
    line = range(ymin, ymax);
    printf("\t%s:\t%s\n", complexity, line);
    free(line);
    line = range(y2min, y2max);
    printf("\tPredicted curvature: %s\n", line);
    free(line);

    /* Here we output to the data file (temp file) */
    //fprintf(data, "load \"params.gnu\"\n");
    fprintf(data, "set terminal png\n");
    fprintf(data, "set title \"%s\"\n", seq_name);
    fprintf(data, "set yrange [ %s : %s ] noreverse nowriteback\n", y2min, y2max);

    if (whichplot == 2) {
        fprintf(data, "set xlabel \"%s ", complexity);

        switch (complexity[0]) {
            case 'B': fprintf(data, "[a.u.]\"\n");
                break; /* Local Bend  */
            case 'G': fprintf(data, "[per cent]\"\n");
                break; /* G+C Content */
            case 'C': fprintf(data, "[a.u.]\"\n");
                break; /* Complexity  */
        }

        fprintf(data, "set xrange [ %s : %s ] noreverse nowriteback\n", ymin, ymax);
    } else {
        fprintf(data, "set xlabel \"Sequence\"\n");
        fprintf(data, "set y2tics border nomirror\n");
        fprintf(data, "set y2label \"%s", complexity);

        switch (complexity[0]) {
            case 'B': fprintf(data, "[a.u.]\"\n");
                break; /* Local Bend  */
            case 'G': fprintf(data, "[per cent]\"\n");
                break; /* G+C Content */
            case 'C': fprintf(data, "[a.u.]\"\n");
                break; /* Complexity  */
        }

        fprintf(data, "set xrange [ %s : %s ] noreverse nowriteback\n", xmin, xmax);
        fprintf(data, "set y2range [ %s : %s ] noreverse nowriteback\n", ymin, ymax);
    }
    fprintf(data, "set ylabel \"Predicted curvature [degrees/10.5 bp helical turn]\"\n");

    switch (whichplot) {
        case 1:
            fprintf(data, "plot \"-\" using 1:2 %s t \"Predicted curvature\" with lines, \"-\" using 1:3 axes x1y2 %s t \"%s\" with lines\n", smth2, smth1, complexity);
            for (i = cwindow; i < seqlen - cwindow + 1; i++) {
                fprintf(data, "%d\t%8.4lf\t%8.4lf\n", i, curve[i], (complexity[0] == 'G' ? gc[i] : bend[i]));
            }
            fprintf(data, "e\n");
            break;

        case 2:
            fprintf(data, "plot \"-\" using 3:2 title \"Predicted curvature vs. %s\"\n", complexity);
            break;
    }

    for (i = cwindow; i < seqlen - cwindow + 1; i++) {
        fprintf(data, "%d\t%8.4lf\t%8.4lf\n", i, curve[i], (complexity[0] == 'G' ? gc[i] : bend[i]));
    }
    fprintf(data, "e\n");
    
    fclose(data);
    
    printf("Running gnuplot\n");
    line  = (char *) malloc(sizeof(char) * (strlen(temp_name) * 2 + 30));
    sprintf(line, "gnuplot < %s > %s.png",temp_name, temp_name);
    system(line);
    free(line);
    
    if (temp_name) free(temp_name);
    if (seq_name) free(seq_name);
    if (seq) free(seq);
    if (scale) free(scale);
    if (curve) free(curve);
    if (bend) free(bend);
    if (gc) free(gc);
    if (smth1) free(smth1);
    if (smth2) free(smth2);
    if (complexity) free(complexity);
    if (xmin) free(xmin);
    if (xmax) free(xmax);
    if (ymin) free(ymin);
    if (ymax) free(ymax);
    if (y2min) free(y2min);
    if (y2max) free(y2max);
    exit(0);
}
