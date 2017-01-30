#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CR 13

char izit(char *kiki) {
    register int x, d;

    d = 1;
    x = 0;
    while (kiki[x] && d) {
        d = (kiki[x] == 'a' || kiki[x] == 'A' ||
                kiki[x] == 't' || kiki[x] == 'T' ||
                kiki[x] == 'g' || kiki[x] == 'G' ||
                kiki[x] == 'c' || kiki[x] == 'C');
        x++;
    }
    return (d ? 0 : kiki[x]);
}

char *makeseq(char *seqr, int g1, int g2, int g3) {
    char *temp = (char *) malloc(sizeof (char) * (strlen(seqr) + 50));
    register int x, y;

    strcpy(temp, "<pre>");
    for (x = 5, y = 0; seqr[y]; x++, y++) {
        temp[x] = seqr[y];
        if (x == g3) {
            temp[x] = '<';
            temp[++x] = 'B';
            temp[++x] = '>';
            temp[++x] = seqr[y];
            temp[++x] = '<';
            temp[++x] = '/';
            temp[++x] = 'B';
            temp[++x] = '>';
        }
        if (!((y + 1) % g1)) temp[++x] = ' ';
        if (!((y + 1) % g2)) temp[++x] = '\n';
    }
    strcat(temp, "</pre>\0");
    return temp;
}

char *range(char *low, char *high) {
    char *temp = (char *) malloc(sizeof (char) * 20);


    if (low[0] == '*') {
        if (high[0] == '*') {
            strcpy(temp, "Automatic");
        } else {
            strcpy(temp, "Max ");
            strcat(temp, high);
        }
    } else {
        if (high[0] == '*') {
            strcpy(temp, "Min ");
            strcat(temp, low);
        } else {
            strcpy(temp, low);
            strcat(temp, " to ");
            strcat(temp, high);
        }
    }

    return temp;
}

void getword(char *word, char *line, char stop) {
    int x = 0, y;

    for (x = 0; ((line[x]) && (line[x] != stop)); x++)
        word[x] = line[x];

    word[x] = '\0';
    if (line[x]) ++x;
    y = 0;

    while (line[y++] = line[x++]);
}

char *makeword(char *line, char stop) {
    int x = 0, y;
    char *word = (char *) malloc(sizeof (char) * (strlen(line) + 1));

    for (x = 0; ((line[x]) && (line[x] != stop)); x++)
        word[x] = line[x];

    word[x] = '\0';
    if (line[x]) ++x;
    y = 0;

    while (line[y++] = line[x++]);
    return word;
}

char *fmakeword(FILE *f, char stop, int *cl) {
    int wsize;
    char *word;
    int ll;

    wsize = 102400;
    ll = 0;
    word = (char *) malloc(sizeof (char) * (wsize + 1));

    while (1) {
        word[ll] = (char) fgetc(f);
        if (ll == wsize) {
            word[ll + 1] = '\0';
            wsize += 102400;
            word = (char *) realloc(word, sizeof (char)*(wsize + 1));
        }
        --(*cl);
        if ((word[ll] == stop) || (feof(f)) || (!(*cl))) {
            if (word[ll] != stop) ll++;
            word[ll] = '\0';
            return word;
        }
        ++ll;
    }
}

char x2c(char *what) {
    register char digit;

    digit = (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A') + 10 : (what[0] - '0'));
    digit *= 16;
    digit += (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A') + 10 : (what[1] - '0'));
    return (digit);
}

void unescape_url(char *url) {
    register int x, y;

    for (x = 0, y = 0; url[y]; ++x, ++y) {
        if ((url[x] = url[y]) == '%') {
            url[x] = x2c(&url[y + 1]);
            y += 2;
        }
    }
    url[x] = '\0';
}

void removespace(char *sq) {
    register int x, y;

    for (x = 0, y = 0; sq[y]; ++y) {
        if (((sq[y] >= 'A') && (sq[y] <= 'Z')) || ((sq[y] >= 'a') && (sq[y] <= 'z'))) {
            sq[x++] = tolower(sq[y]);
        }
    }
    sq[x] = '\0';
}

void plustospace(char *str) {
    register int x;

    for (x = 0; str[x]; x++) if (str[x] == '+') str[x] = ' ';
}

int rind(char *s, char c) {
    register int x;
    for (x = strlen(s) - 1; x != -1; x--)
        if (s[x] == c) return x;
    return -1;
}

/* commented, because interferes with NAB getline 

int getline(char *s, int n, FILE *f) {
    register int i=0;

    while(1) {
        s[i] = (char)fgetc(f);

        if(s[i] == CR)
            s[i] = fgetc(f);

        if((s[i] == 0x4) || (s[i] == LF) || (i == (n-1))) {
            s[i] = '\0';
            return (feof(f) ? 1 : 0);
        }
        ++i;
    }
}
 */


void send_fd(FILE *f, FILE *fd) {
    int num_chars = 0;
    char c;

    while (1) {
        c = fgetc(f);
        if (feof(f))
            return;
        fputc(c, fd);
    }
}

int ind(char *s, char c) {
    register int x;

    for (x = 0; s[x]; x++)
        if (s[x] == c) return x;

    return -1;
}

void escape_shell_cmd(char *cmd) {
    register int x, y, l;

    l = strlen(cmd);
    for (x = 0; cmd[x]; x++) {
        if (ind("&;`'\"|*?~<>^()[]{}$\\", cmd[x]) != -1) {
            for (y = l + 1; y > x; y--)
                cmd[y] = cmd[y - 1];
            l++; /* length has been increased */
            cmd[x] = '\\';
            x++; /* skip the character */
        }
    }
}

