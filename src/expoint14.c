#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lnfac.h"
#include "lnnuc.h"

/* Minden ki van iratva! Altalanos a modszer! */
/* Egy javitas lehetne meg, nagyobb szomeretnel ideiglenes tomb! >2 eseten */

#define NNUM    4

void load_seq(char *seq, int *seq_num);

double k1_formula(unsigned int *state, unsigned int *signal, unsigned int winlen, unsigned int siglen, unsigned int word) {
    double total, totseq, ans, k1;
    unsigned int i;

    totseq = ((double) winlen) * lnnuc[word];
    ans = lnfac[winlen];
    for (i = 0; i < siglen; i++) ans -= lnfac[*(state + (*(signal + i)))];
    k1 = (double) ans / totseq;
    return ( k1);
}

void count_state_n(unsigned int *seq, unsigned int *state, unsigned int *signal, unsigned int winlen, unsigned int word, unsigned int *siglen, unsigned int poz) {
    unsigned int j, k, z, v, w, len;

    /* Kiszedes a gyujtobol. */
    for (j = 0, len = 0, z = 1, v = poz - 1; j < word; j++, z *= NNUM, v++) {
        if (seq[v] >= 2) {
            if (seq[v] > 2) len += 3 * z;
            else len += 2 * z;
        } else {
            if (seq[v] > 0) len += 1 * z;
        }
    }
    (*(state + len)) -= 1;
    if ((*(state + len)) == 0) {
        for (j = 0; (*(signal + j)) != len; j++);
        for (k = j; k < (*siglen); k++) {
            (*(signal + k)) = (*(signal + k + 1));
        }
        (*siglen)--;
    }
    /* Hozzaadas a gyujtohoz. */
    for (j = 0, len = 0, z = 1, v = poz + winlen - word; j < word; j++, z *= NNUM, v++) {
        if (seq[v] >= 2) {
            if (seq[v] > 2) len += 3 * z;
            else len += 2 * z;
        } else {
            if (seq[v] > 0) len += 1 * z;
        }
    }
    (*(state + len)) += 1;
    if ((*(state + len)) == 1) {
        (*(signal + (*siglen))) = len;
        (*siglen)++;
    }
}

void first_count_state_n(unsigned int *seq, unsigned int *state, unsigned int *signal, unsigned int winlen, unsigned int word, unsigned int *siglen) {
    unsigned int i, j, z, v, w, len;


    for (i = 0, w = 0; i < winlen - word + 1; i++) {
        for (j = 0, len = 0, z = 1, v = i; j < word; j++, z *= NNUM, v++) {
            if (seq[v] >= 2) {
                if (seq[v] > 2) len += 3 * z;
                else len += 2 * z;
            } else {
                if (seq[v] > 0) len += 1 * z;
            }
        }
        (*(state + len)) += 1;
        if ((*(state + len)) == 1) {
            (*(signal + w)) = len;
            w++;
        }
    }
    *siglen = w;
}

int k1_complexity_n(char *seqa, unsigned int seqlen, double *k1values, unsigned int winlen, unsigned int word) {

    unsigned int len, i, j, st, stlen, siglen;
    unsigned int *stateptr, *signalptr, *seq;
    double *work;

    seq = (unsigned int *) malloc(sizeof (unsigned int)*strlen(seqa));

    load_seq(seqa, (int *) seq);

    if ((winlen > seqlen || winlen < word)) return ( 1);
    for (i = 0, stlen = 1; i < word; i++) stlen *= NNUM;
    stateptr = (unsigned int *) malloc((stlen + 1) * sizeof (unsigned int));
    signalptr = (unsigned int *) malloc((stlen + 1) * sizeof (unsigned int));
    work = k1values;
    len = seqlen - winlen;
    for (i = 0; i < stlen; i++) {
        *(stateptr + i) = 0;
        *(signalptr + i) = 0;
    }
    first_count_state_n(seq, stateptr, signalptr, winlen, word, &siglen);
    *work = k1_formula(stateptr, signalptr, winlen, siglen, word);
    work++;
    for (i = 1; i <= len; i++) {
        count_state_n(seq, stateptr, signalptr, winlen, word, &siglen, i);
        *work = k1_formula(stateptr, signalptr, winlen, siglen, word);
        work++;
    }

    free(signalptr);
    free(stateptr);
    free(seq);

    return ( 0);
}

/*
double k1[200];

void main(void) {

 unsigned int winl,i,word;
 unsigned int aa[200]= { 0,1,2,3,0,1,2,3,2,3,1,2,3,2,3,2,3,1,0,0 };
 int *hm,t,len; 
 
 len=15;
 winl=8;
 word=4;
 hm=&t;
 *hm=k1_complexity_n(aa,len,k1,winl,word);
 printf ("UTAN IME:\n");
 for(i=0;aa[i]!='\0';i++) {
  printf ("%d",aa[i]);
 }
 printf ("\n");
 printf ("UTAN A K1:\n");
 for(i=0;i<200;i++) {
  printf ("%f\n",k1[i]);
 }
 printf ("\n");
}

 */
