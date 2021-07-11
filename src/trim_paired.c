#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include "sickle.h"
#include "kseq.h"
#include "print_record.h"

#include "kthread.h"

__KS_GETC(gzread, BUFFER_SIZE)
__KS_GETUNTIL(gzread, BUFFER_SIZE)
__KSEQ_READ

int paired_qual_threshold = 20;
int paired_length_threshold = 20;

static struct option paired_long_options[] = {
    {"qual-type", required_argument, 0, 't'},
    {"pe-file1", required_argument, 0, 'f'},
    {"pe-file2", required_argument, 0, 'r'},
    {"pe-combo", required_argument, 0, 'c'},
    {"output-pe1", required_argument, 0, 'o'},
    {"output-pe2", required_argument, 0, 'p'},
    {"output-single", required_argument, 0, 's'},
    {"output-combo", required_argument, 0, 'm'},
    {"qual-threshold", required_argument, 0, 'q'},
    {"length-threshold", required_argument, 0, 'l'},
    {"no-fiveprime", no_argument, 0, 'x'},
    {"truncate-n", no_argument, 0, 'n'},
    {"gzip-output", no_argument, 0, 'g'},
    {"output-combo-all", required_argument, 0, 'M'},
    {"quiet", no_argument, 0, 'z'},
    {"n_thread", required_argument, 0, 'u'},
    {"p_thread", required_argument, 0, 'U'},
    {"block_size", required_argument, 0, 'b'},    
    {GETOPT_HELP_OPTION_DECL},
    {GETOPT_VERSION_OPTION_DECL},
    {NULL, 0, NULL, 0}
};

void paired_usage (int status, char *msg) {

    fprintf(stderr, "\nIf you have separate files for forward and reverse reads:\n");
    fprintf(stderr, "Usage: %s pe [options] -f <paired-end forward fastq file> -r <paired-end reverse fastq file> -t <quality type> -o <trimmed PE forward file> -p <trimmed PE reverse file> -s <trimmed singles file>\n\n", PROGRAM_NAME);
    fprintf(stderr, "If you have one file with interleaved forward and reverse reads:\n");
    fprintf(stderr, "Usage: %s pe [options] -c <interleaved input file> -t <quality type> -m <interleaved trimmed paired-end output> -s <trimmed singles file>\n\n\
If you have one file with interleaved reads as input and you want ONLY one interleaved file as output:\n\
Usage: %s pe [options] -c <interleaved input file> -t <quality type> -M <interleaved trimmed output>\n\n", PROGRAM_NAME, PROGRAM_NAME);
    fprintf(stderr, "Options:\n\
Paired-end separated reads\n\
--------------------------\n\
-f, --pe-file1, Input paired-end forward fastq file (Input files must have same number of records)\n\
-r, --pe-file2, Input paired-end reverse fastq file\n\
-o, --output-pe1, Output trimmed forward fastq file\n\
-p, --output-pe2, Output trimmed reverse fastq file. Must use -s option.\n\
-u, --n_thread,   Number of worker threads to used [default:1]\n\
-U, --p_thread,   Number of pipeline threads to used [default:3]\n\
-b, --block_size, Number of sequences to read per time per thread [default:20000000]\n\n\
Paired-end interleaved reads\n\
----------------------------\n");
    fprintf(stderr,"-c, --pe-combo, Combined (interleaved) input paired-end fastq\n\
-m, --output-combo, Output combined (interleaved) paired-end fastq file. Must use -s option.\n\
-M, --output-combo-all, Output combined (interleaved) paired-end fastq file with any discarded read written to output file as a single N. Cannot be used with the -s option.\n\n\
Global options\n\
--------------\n\
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)\n");
    fprintf(stderr, "-s, --output-single, Output trimmed singles fastq file\n\
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
-x, --no-fiveprime, Don't do five prime trimming.\n\
-n, --truncate-n, Truncate sequences at position of first N.\n");


    fprintf(stderr, "-g, --gzip-output, Output gzipped files.\n--quiet, do not output trimming info\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

    if (msg) fprintf(stderr, "%s\n\n", msg);
    exit(status);
}

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

typedef struct { // global data structure for kt_pipeline()
    int block_len; // read at most block_len for each pipe
    int n_thread;  // use n_threads for trimming reads
    kseq_t *ks[2]; // open 2 kseq_t files for reading PE
    int qualtype;
    int paired_length_threshold;
    int paired_qual_threshold;
    int no_fiveprime;
    int trunc_n;
    int debug;
    gzFile pec;

    int combo_all;
    int gzip_output;
    FILE *outfile1; /* forward output file handle */
    FILE *outfile2; /* reverse output file handle */
    FILE *combo;    /* combined output file handle */
    FILE *single;   /* single output file handle */
    gzFile outfile1_gzip;
    gzFile outfile2_gzip;
    gzFile combo_gzip;
    gzFile single_gzip;

    //Stats
    int total;
    int kept_p;
    int kept_s1;
    int kept_s2;
    int discard_p;
    int discard_s1;
    int discard_s2;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
    pldat_t *p;
    int n, m, sum_len;
    kseq_t *kseq; // Sequence data
    cutsites **pcut; // trim results
    //buf_c4_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
    stepdat_t *s = (stepdat_t*)data;

    // printf("%s\n", s->kseq[i].seq.s);
    //printf("%i working on %li\n", tid, i);
    s->pcut[i] = sliding_window(&s->kseq[i], s->p->qualtype, s->p->paired_length_threshold, s->p->paired_qual_threshold, s->p->no_fiveprime, s->p->trunc_n, s->p->debug);
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    pldat_t *p = (pldat_t*)data;
    if (step == 0) { // step 1: read a block of sequences
        int ret1, ret2;
        stepdat_t *s;
        CALLOC(s, 1);
        s->p = p;
        while ((ret1 = kseq_read(p->ks[0])) >= 0) {
            
            ret2 = kseq_read(p->ks[1]);
            if (ret2 < 0) {
                fprintf(stderr, "Warning: PE file 2 is shorter than PE file 1. Disregarding rest of PE file 1.\n");
                break;
            }

            // Check if we have enough space left for at least 2 entries
            if (s->n >= s->m - 2) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->kseq, s->m);
            }
            //Store PE1
            //Copy name
            MALLOC(s->kseq[s->n].name.s, p->ks[0]->name.l);
            memcpy(s->kseq[s->n].name.s, p->ks[0]->name.s, p->ks[0]->name.l);
            s->kseq[s->n].name.l = p->ks[0]->name.l;
            s->kseq[s->n].name.s[p->ks[0]->name.l] = '\0';

            //Copy comment
            if (p->ks[0]->comment.l) {
                MALLOC(s->kseq[s->n].comment.s, p->ks[0]->comment.l);
                memcpy(s->kseq[s->n].comment.s, p->ks[0]->comment.s, p->ks[0]->comment.l);
                s->kseq[s->n].comment.l = p->ks[0]->comment.l;
                s->kseq[s->n].comment.s[p->ks[0]->comment.l] = '\0';
            }

            //Copy sequence
            MALLOC(s->kseq[s->n].seq.s, p->ks[0]->seq.l);
            memcpy(s->kseq[s->n].seq.s, p->ks[0]->seq.s, p->ks[0]->seq.l);
            s->kseq[s->n].seq.l = p->ks[0]->seq.l;
            s->kseq[s->n].seq.s[p->ks[0]->seq.l] = '\0';

            s->sum_len += p->ks[0]->seq.l;

            //Copy quality
            MALLOC(s->kseq[s->n].qual.s, p->ks[0]->qual.l);
            memcpy(s->kseq[s->n].qual.s, p->ks[0]->qual.s, p->ks[0]->qual.l);
            s->kseq[s->n].qual.l = p->ks[0]->qual.l;
            s->kseq[s->n].qual.s[p->ks[0]->qual.l] = '\0';

            //Increment to store PE read
            s->n++;

            //Store PE2
            //Copy name
            MALLOC(s->kseq[s->n].name.s, p->ks[1]->name.l);
            memcpy(s->kseq[s->n].name.s, p->ks[1]->name.s, p->ks[1]->name.l);
            s->kseq[s->n].name.l = p->ks[1]->name.l;
            s->kseq[s->n].name.s[p->ks[1]->name.l] = '\0';

            //Copy comment
            if (p->ks[0]->comment.l) {
                MALLOC(s->kseq[s->n].comment.s, p->ks[1]->comment.l);
                memcpy(s->kseq[s->n].comment.s, p->ks[1]->comment.s, p->ks[1]->comment.l);
                s->kseq[s->n].comment.l = p->ks[1]->comment.l;
                s->kseq[s->n].comment.s[p->ks[1]->comment.l] = '\0';
            }
            
            //Copy sequence
            MALLOC(s->kseq[s->n].seq.s, p->ks[1]->seq.l);
            memcpy(s->kseq[s->n].seq.s, p->ks[1]->seq.s, p->ks[1]->seq.l);
            s->kseq[s->n].seq.l = p->ks[1]->seq.l;
            s->kseq[s->n].seq.s[p->ks[1]->seq.l] = '\0';

            s->sum_len += p->ks[1]->seq.l;

            //Copy quality
            MALLOC(s->kseq[s->n].qual.s, p->ks[1]->qual.l);
            memcpy(s->kseq[s->n].qual.s, p->ks[1]->qual.s, p->ks[1]->qual.l);
            s->kseq[s->n].qual.l = p->ks[1]->qual.l;
            s->kseq[s->n].qual.s[p->ks[1]->qual.l] = '\0';

            //Increment to store next PE read
            s->n++;

            // Filled block_size with sequence data
            if (s->sum_len >= p->block_len) {
                break;
            }
        }
        // No sequences left
        if (s->sum_len == 0) {
            free(s);
        // Filled a block_size with data to continue with step 2
        } else {
            //fprintf(stderr, "%s\n", s->kseq[s->n-1].name.s);
            s->pcut = calloc(s->n, sizeof(struct cutsites *));
            s->p->total += s->n; // add to total reads
            fprintf(stderr, "\rRead: PE %i sequences (total %i)", s->n/2, s->p->total/2);
            return s;
        }
    } else if (step == 1) { // step 2: trim sequences
        stepdat_t *s = (stepdat_t*)in;

        kt_for(p->n_thread, worker_for, s, s->n);

        return s;
    } else if (step == 2) { // step 3: Output files
        stepdat_t *s = (stepdat_t*)in;

        for (int i = 0; i < s->n - 1; i += 2) {

            if (s->p->debug) {
                printf("p1cut: %d,%d\n", s->pcut[i]->five_prime_cut, s->pcut[i]->three_prime_cut);
                printf("p2cut: %d,%d\n", s->pcut[i+1]->five_prime_cut, s->pcut[i+1]->three_prime_cut);
            }

            /* The sequence and quality print statements below print out the sequence string starting from the 5' cut */
            /* and then only print out to the 3' cut, however, we need to adjust the 3' cut */
            /* by subtracting the 5' cut because the 3' cut was calculated on the original sequence */

            /* if both sequences passed quality and length filters, then output both records */
            if (s->pcut[i]->three_prime_cut >= 0 && s->pcut[i+1]->three_prime_cut >= 0) {
                if (!s->p->gzip_output) {
                    if (s->p->pec) {
                        print_record(s->p->combo, &s->kseq[i], s->pcut[i]);
                        print_record(s->p->combo, &s->kseq[i+1], s->pcut[i+1]);
                    } else {
                        print_record(s->p->outfile1, &s->kseq[i], s->pcut[i]);
                        print_record(s->p->outfile2, &s->kseq[i+1], s->pcut[i+1]);
                    }
                } else {
                    if (s->p->pec) {
                        print_record_gzip (s->p->combo_gzip, &s->kseq[i], s->pcut[i]);
                        print_record_gzip (s->p->combo_gzip, &s->kseq[i+1], s->pcut[i+1]);
                    } else {
                        print_record_gzip (s->p->outfile1_gzip, &s->kseq[i], s->pcut[i]);
                        print_record_gzip (s->p->outfile2_gzip, &s->kseq[i+1], s->pcut[i+1]);
                    }
                }

                s->p->kept_p += 2;
            /* if only one sequence passed filter, then put its record in singles and discard the other */
            /* or put an "N" record in if that option was chosen. */
            } else if (s->pcut[i]->three_prime_cut >= 0 && s->pcut[i+1]->three_prime_cut < 0) {
                    if (!s->p->gzip_output) {
                        if (s->p->combo_all) {
                            print_record (s->p->combo, &s->kseq[i], s->pcut[i]);
                            print_record_N (s->p->combo, &s->kseq[i+1], s->p->qualtype);
                        } else {
                            print_record (s->p->single, &s->kseq[i], s->pcut[i]);
                        }
                    } else {
                        if (s->p->combo_all) {
                            print_record_gzip (s->p->combo_gzip, &s->kseq[i], s->pcut[i]);
                            print_record_N_gzip (s->p->combo_gzip, &s->kseq[i+1], s->p->qualtype);
                        } else {
                            print_record_gzip (s->p->single_gzip, &s->kseq[i], s->pcut[i]);
                        }
                    }

                    s->p->kept_s1++;
                    s->p->discard_s2++;
            } else if (s->pcut[i]->three_prime_cut < 0 && s->pcut[i+1]->three_prime_cut >= 0) {
                if (!s->p->gzip_output) {
                    if (s->p->combo_all) {
                        print_record_N (s->p->combo, &s->kseq[i], s->p->qualtype);
                        print_record (s->p->combo, &s->kseq[i+1], s->pcut[i+1]);
                    } else {
                        print_record (s->p->single, &s->kseq[i+1], s->pcut[i+1]);
                    }
                } else {
                    if (s->p->combo_all) {
                        print_record_N_gzip (s->p->combo_gzip, &s->kseq[i], s->p->qualtype);
                        print_record_gzip (s->p->combo_gzip, &s->kseq[i+1], s->pcut[i+1]);
                    } else {
                        print_record_gzip (s->p->single_gzip, &s->kseq[i+1], s->pcut[i+1]);
                    }
                }

                s->p->kept_s2++;
                s->p->discard_s1++;
            } else {

                /* If both records are to be discarded, but the -M option */
                /* is being used, then output two "N" records */
                if (s->p->combo_all) {
                    if (!s->p->gzip_output) {
                        print_record_N (s->p->combo, &s->kseq[i], s->p->qualtype);
                        print_record_N (s->p->combo, &s->kseq[i+1], s->p->qualtype);
                    } else {
                        print_record_N_gzip (s->p->combo_gzip, &s->kseq[i], s->p->qualtype);
                        print_record_N_gzip (s->p->combo_gzip, &s->kseq[i+1], s->p->qualtype);
                    }
                }

                s->p->discard_p += 2;
            }
        }

        //Clean up
        for (int i = 0; i < s->n; ++i) {
            free(s->kseq[i].name.s);
            free(s->kseq[i].comment.s);
            free(s->kseq[i].seq.s);
            free(s->kseq[i].qual.s);

            free(s->pcut[i]);
        }
        free(s->kseq);

        free(s->pcut);

        free(s);
    }
    return 0;
}

int paired_main(int argc, char *argv[]) {

    gzFile pe1 = NULL;          /* forward input file handle */
    gzFile pe2 = NULL;          /* reverse input file handle */
    gzFile pec = NULL;          /* combined input file handle */
    kseq_t *fqrec1 = NULL;
    kseq_t *fqrec2 = NULL;
    int l1, l2;
    FILE *outfile1 = NULL;      /* forward output file handle */
    FILE *outfile2 = NULL;      /* reverse output file handle */
    FILE *combo = NULL;         /* combined output file handle */
    FILE *single = NULL;        /* single output file handle */
    gzFile outfile1_gzip = NULL;
    gzFile outfile2_gzip = NULL;
    gzFile combo_gzip = NULL;
    gzFile single_gzip = NULL;
    int debug = 0;
    int optc;
    extern char *optarg;
    int qualtype = -1;
    cutsites *p1cut;
    cutsites *p2cut;
    char *outfn1 = NULL;        /* forward file out name */
    char *outfn2 = NULL;        /* reverse file out name */
    char *outfnc = NULL;        /* combined file out name */
    char *sfn = NULL;           /* single/combined file out name */
    char *infn1 = NULL;         /* forward input filename */
    char *infn2 = NULL;         /* reverse input filename */
    char *infnc = NULL;         /* combined input filename */
    int kept_p = 0;
    int discard_p = 0;
    int kept_s1 = 0;
    int kept_s2 = 0;
    int discard_s1 = 0;
    int discard_s2 = 0;
    int quiet = 0;
    int no_fiveprime = 0;
    int trunc_n = 0;
    int gzip_output = 0;
    int combo_all=0;
    int combo_s=0;
    int total=0;
    int n_thread = 1; // worker threads
    int p_thread = 3; // pipeline threads
    int block_size = 20000000; // Total size of the sequence loaded per thread

    while (1) {
        int option_index = 0;
        optc = getopt_long(argc, argv, "df:r:c:t:o:p:m:M:u:U:b:s:q:l:xng", paired_long_options, &option_index);

        if (optc == -1)
            break;

        switch (optc) {
            if (paired_long_options[option_index].flag != 0)
                break;

        case 'f':
            infn1 = (char *) malloc(strlen(optarg) + 1);
            strcpy(infn1, optarg);
            break;

        case 'r':
            infn2 = (char *) malloc(strlen(optarg) + 1);
            strcpy(infn2, optarg);
            break;
        
        case 'u':
            n_thread = atoi(optarg);
            break;
        case 'U':
            p_thread = atoi(optarg);
            break;
        case 'b':
            block_size = atoi(optarg);
            break;

        case 'c':
            infnc = (char *) malloc(strlen(optarg) + 1);
            strcpy(infnc, optarg);
            break;

        case 't':
            if (!strcmp(optarg, "illumina")) qualtype = ILLUMINA;
            else if (!strcmp(optarg, "solexa")) qualtype = SOLEXA;
            else if (!strcmp(optarg, "sanger")) qualtype = SANGER;
            else {
                fprintf(stderr, "Error: Quality type '%s' is not a valid type.\n", optarg);
                return EXIT_FAILURE;
            }
            break;

        case 'o':
            outfn1 = (char *) malloc(strlen(optarg) + 1);
            strcpy(outfn1, optarg);
            break;

        case 'p':
            outfn2 = (char *) malloc(strlen(optarg) + 1);
            strcpy(outfn2, optarg);
            break;

        case 'm':
            outfnc = (char *) malloc(strlen(optarg) + 1);
            strcpy(outfnc, optarg);
            combo_s = 1;
            break;

        case 'M':
            outfnc = (char *) malloc(strlen(optarg) + 1);
            strcpy(outfnc, optarg);
            combo_all = 1;
            break;

        case 's':
            sfn = (char *) malloc(strlen(optarg) + 1);
            strcpy(sfn, optarg);
            break;

        case 'q':
            paired_qual_threshold = atoi(optarg);
            if (paired_qual_threshold < 0) {
                fprintf(stderr, "Quality threshold must be >= 0\n");
                return EXIT_FAILURE;
            }
            break;

        case 'l':
            paired_length_threshold = atoi(optarg);
            if (paired_length_threshold < 0) {
                fprintf(stderr, "Length threshold must be >= 0\n");
                return EXIT_FAILURE;
            }
            break;

        case 'x':
            no_fiveprime = 1;
            break;

        case 'n':
            trunc_n = 1;
            break;

        case 'g':
            gzip_output = 1;
            break;

        case 'z':
            quiet = 1;
            break;

        case 'd':
            debug = 1;
            break;

        case_GETOPT_HELP_CHAR(paired_usage);
        case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

        case '?':
            paired_usage(EXIT_FAILURE, NULL);
            break;

        default:
            paired_usage(EXIT_FAILURE, NULL);
            break;
        }
    }

    /* required: qualtype */
    if (qualtype == -1) {
        paired_usage(EXIT_FAILURE, "****Error: Quality type is required.");
    }

    /* make sure minimum input filenames are specified */
    if (!infn1 && !infnc) {
        paired_usage(EXIT_FAILURE, "****Error: Must have either -f OR -c argument.");
    }

    if (infnc) {      /* using combined input file */

        if (infn1 || infn2 || outfn1 || outfn2) {
            paired_usage(EXIT_FAILURE, "****Error: Cannot have -f, -r, -o, or -p options with -c.");
        }

        if ((combo_all && combo_s) || (!combo_all && !combo_s)) {
            paired_usage(EXIT_FAILURE, "****Error: Must have only one of either -m or -M options with -c.");
        }

        if ((combo_s && !sfn) || (combo_all && sfn)) {
            paired_usage(EXIT_FAILURE, "****Error: -m option must have -s option, and -M option cannot have -s option.");
        }

        /* check for duplicate file names */
        if (!strcmp(infnc, outfnc) || (combo_s && (!strcmp(infnc, sfn) || !strcmp(outfnc, sfn)))) {
            fprintf(stderr, "****Error: Duplicate filename between combo input, combo output, and/or single output file names.\n\n");
            return EXIT_FAILURE;
        }

        /* get combined output file */
        if (!gzip_output) {
            combo = fopen(outfnc, "w");
            if (!combo) {
                fprintf(stderr, "****Error: Could not open combo output file '%s'.\n\n", outfnc);
                return EXIT_FAILURE;
            }
        } else {
            combo_gzip = gzopen(outfnc, "w");
            if (!combo_gzip) {
                fprintf(stderr, "****Error: Could not open combo output file '%s'.\n\n", outfnc);
                return EXIT_FAILURE;
            }
        }

        pec = gzopen(infnc, "r");
        if (!pec) {
            fprintf(stderr, "****Error: Could not open combined input file '%s'.\n\n", infnc);
            return EXIT_FAILURE;
        }

    } else {     /* using forward and reverse input files */

        if (infn1 && (!infn2 || !outfn1 || !outfn2 || !sfn)) {
            paired_usage(EXIT_FAILURE, "****Error: Using the -f option means you must have the -r, -o, -p, and -s options.");
        }

        if (infn1 && (infnc || combo_all || combo_s)) {
            paired_usage(EXIT_FAILURE, "****Error: The -f option cannot be used in combination with -c, -m, or -M.");
        }

        if (!strcmp(infn1, infn2) || !strcmp(infn1, outfn1) || !strcmp(infn1, outfn2) ||
            !strcmp(infn1, sfn) || !strcmp(infn2, outfn1) || !strcmp(infn2, outfn2) || 
            !strcmp(infn2, sfn) || !strcmp(outfn1, outfn2) || !strcmp(outfn1, sfn) || !strcmp(outfn2, sfn)) {

            fprintf(stderr, "****Error: Duplicate input and/or output file names.\n\n");
            return EXIT_FAILURE;
        }

        pe1 = gzopen(infn1, "r");
        if (!pe1) {
            fprintf(stderr, "****Error: Could not open input file '%s'.\n\n", infn1);
            return EXIT_FAILURE;
        }

        pe2 = gzopen(infn2, "r");
        if (!pe2) {
            fprintf(stderr, "****Error: Could not open input file '%s'.\n\n", infn2);
            return EXIT_FAILURE;
        }

        if (!gzip_output) {
            outfile1 = fopen(outfn1, "w");
            if (!outfile1) {
                fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn1);
                return EXIT_FAILURE;
            }

            outfile2 = fopen(outfn2, "w");
            if (!outfile2) {
                fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn2);
                return EXIT_FAILURE;
            }
        } else {
            outfile1_gzip = gzopen(outfn1, "w");
            if (!outfile1_gzip) {
                fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn1);
                return EXIT_FAILURE;
            }

            outfile2_gzip = gzopen(outfn2, "w");
            if (!outfile2_gzip) {
                fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn2);
                return EXIT_FAILURE;
            }

        }
    }

    /* get singles output file handle */
    if (sfn && !combo_all) {
        if (!gzip_output) {
            single = fopen(sfn, "w");
            if (!single) {
                fprintf(stderr, "****Error: Could not open single output file '%s'.\n\n", sfn);
                return EXIT_FAILURE;
            }
        } else {
            single_gzip = gzopen(sfn, "w");
            if (!single_gzip) {
                fprintf(stderr, "****Error: Could not open single output file '%s'.\n\n", sfn);
                return EXIT_FAILURE;
            }
        }
    }

    if (pec) {
        fqrec1 = kseq_init(pec);
        fqrec2 = (kseq_t *) malloc(sizeof(kseq_t));
        fqrec2->f = fqrec1->f;
    } else {
        fqrec1 = kseq_init(pe1);
        fqrec2 = kseq_init(pe2);
    }

    pldat_t pl;
    if (pec) {
        pl.ks[0] = kseq_init(pec);
        pl.ks[1] = pl.ks[0];
    } else {
        pl.ks[0] = kseq_init(pe1);
        pl.ks[1] = kseq_init(pe2);
    }
    pl.n_thread = n_thread;
    pl.block_len = block_size;
    pl.qualtype = qualtype;
    pl.paired_length_threshold = paired_length_threshold;
    pl.paired_qual_threshold = paired_qual_threshold;
    pl.no_fiveprime = no_fiveprime;
    pl.trunc_n = trunc_n;
    pl.debug = debug;
    pl.pec = pec;

    pl.combo_all = combo_all;
    pl.gzip_output = gzip_output;
    pl.outfile1 = outfile1;
    pl.outfile2 = outfile2;
    pl.combo = combo;
    pl.single = single;
    pl.outfile1_gzip = outfile1_gzip;
    pl.outfile2_gzip = outfile2_gzip;
    pl.combo_gzip = combo_gzip;
    pl.single_gzip = single_gzip;

    //Stats
    pl.total = 0;
    pl.kept_p = 0;
    pl.kept_s1 = 0;
    pl.kept_s2 = 0;
    pl.discard_p = 0;
    pl.discard_s1 = 0;
    pl.discard_s2 = 0;

    kt_pipeline(p_thread, worker_pipeline, &pl, 3);

    if (kseq_read(pl.ks[0]) < 0) {
        l2 = kseq_read(pl.ks[1]);
        if (l2 >= 0) {
            fprintf(stderr, "Warning: PE file 1 is shorter than PE file 2. Disregarding rest of PE file 2.\n");
        }
    }

    if (pec) {
        kseq_destroy(pl.ks[0]);
    } else {
        kseq_destroy(pl.ks[0]);
        kseq_destroy(pl.ks[1]);
    }

    if (!quiet) {
        if (infn1 && infn2) fprintf(stdout, "\nPE forward file: %s\nPE reverse file: %s\n", infn1, infn2);
        if (infnc) fprintf(stdout, "\nPE interleaved file: %s\n", infnc);
        fprintf(stdout, "\nTotal input FastQ records: %d (%d pairs)\n", pl.total, (pl.total / 2));
        fprintf(stdout, "\nFastQ paired records kept: %d (%d pairs)\n", pl.kept_p, (pl.kept_p / 2));
        if (pec) fprintf(stdout, "FastQ single records kept: %d\n", (pl.kept_s1 + pl.kept_s2));
        else fprintf(stdout, "FastQ single records kept: %d (from PE1: %d, from PE2: %d)\n", (pl.kept_s1 + pl.kept_s2), pl.kept_s1, pl.kept_s2);

        fprintf(stdout, "FastQ paired records discarded: %d (%d pairs)\n", pl.discard_p, (pl.discard_p / 2));

        if (pec) fprintf(stdout, "FastQ single records discarded: %d\n\n", (pl.discard_s1 + pl.discard_s2));
        else fprintf(stdout, "FastQ single records discarded: %d (from PE1: %d, from PE2: %d)\n\n", (pl.discard_s1 + pl.discard_s2), pl.discard_s1, pl.discard_s2);
    }

    kseq_destroy(fqrec1);
    if (pec) free(fqrec2);
    else kseq_destroy(fqrec2);

    if (sfn && !combo_all) {
        if (!gzip_output) fclose(single);
        else gzclose(single_gzip);
    }

    if (pec) {
        gzclose(pec);
        if (!gzip_output) fclose(combo);
        else gzclose(combo_gzip);
    } else {
        gzclose(pe1);
        gzclose(pe2);
        if (!gzip_output) {
            fclose(outfile1);
            fclose(outfile2);
        } else {
            gzclose(outfile1_gzip);
            gzclose(outfile2_gzip);
        }
    }

    return EXIT_SUCCESS;
}                               /* end of paired_main() */
