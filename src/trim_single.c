#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include "sickle.h"
#include "kseq.h"
#include "print_record.h"

#include "kthread.h"

__KS_GETC(gzread, BUFFER_SIZE)
__KS_GETUNTIL(gzread, BUFFER_SIZE)
__KSEQ_READ

int single_qual_threshold = 20;
int single_length_threshold = 20;

static struct option single_long_options[] = {
    {"fastq-file", required_argument, 0, 'f'},
    {"output-file", required_argument, 0, 'o'},
    {"qual-type", required_argument, 0, 't'},
    {"qual-threshold", required_argument, 0, 'q'},
    {"length-threshold", required_argument, 0, 'l'},
    {"no-fiveprime", no_argument, 0, 'x'},
    {"discard-n", no_argument, 0, 'n'},
    {"gzip-output", no_argument, 0, 'g'},
    {"quiet", no_argument, 0, 'z'},
    {"n_thread", required_argument, 0, 'u'},
    {"p_thread", required_argument, 0, 'U'},
    {"block_size", required_argument, 0, 'b'},    
    {GETOPT_HELP_OPTION_DECL},
    {GETOPT_VERSION_OPTION_DECL},
    {NULL, 0, NULL, 0}
};

void single_usage(int status, char *msg) {

    fprintf(stderr, "\nUsage: %s se [options] -f <fastq sequence file> -t <quality type> -o <trimmed fastq file>\n\
\n\
Options:\n\
-f, --fastq-file, Input fastq file (required)\n\
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)\n\
-u, --n_thread,   Number of worker threads to used [default:1]\n\
-U, --p_thread,   Number of pipeline threads to used [default:3]\n\
-b, --block_size, Number of sequences to read per time per thread [default:20000000]\n\n\
-o, --output-file, Output trimmed fastq file (required)\n", PROGRAM_NAME);

    fprintf(stderr, "-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
-x, --no-fiveprime, Don't do five prime trimming.\n\
-n, --trunc-n, Truncate sequences at position of first N.\n\
-g, --gzip-output, Output gzipped files.\n\
--quiet, Don't print out any trimming information\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

    if (msg) fprintf(stderr, "%s\n\n", msg);
    exit(status);
}

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

typedef struct {   // global data structure for kt_pipeline()
    int block_len; // read at most block_len for each pipe
    int n_thread;  // use n_threads for trimming reads
    kseq_t *ks;    // open kseq_t file for reading SE
    int qualtype;
    int single_length_threshold;
    int single_qual_threshold;
    int no_fiveprime;
    int trunc_n;
    int debug;
    gzFile pec;

    int combo_all;
    int gzip_output;
    FILE *outfile; /* forward output file handle */
    gzFile outfile_gzip;

    //Stats
    int total;
    int kept;
    int discard;
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
    s->pcut[i] = sliding_window(&s->kseq[i], s->p->qualtype, s->p->single_length_threshold, s->p->single_qual_threshold, s->p->no_fiveprime, s->p->trunc_n, s->p->debug);
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    pldat_t *p = (pldat_t*)data;
    if (step == 0) { // step 1: read a block of sequences
        int ret1, ret2;
        stepdat_t *s;
        CALLOC(s, 1);
        s->p = p;
        while ((ret1 = kseq_read(p->ks)) >= 0) {

            // Check if we have enough space left for at least 2 entries
            if (s->n >= s->m - 2) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->kseq, s->m);
            }
            //Store PE1
            //Copy name
            MALLOC(s->kseq[s->n].name.s, p->ks->name.l);
            memcpy(s->kseq[s->n].name.s, p->ks->name.s, p->ks->name.l);
            s->kseq[s->n].name.l = p->ks->name.l;
            s->kseq[s->n].name.s[p->ks->name.l] = '\0';

            //Copy comment
            if (p->ks->comment.l) {
                MALLOC(s->kseq[s->n].comment.s, p->ks->comment.l);
                memcpy(s->kseq[s->n].comment.s, p->ks->comment.s, p->ks->comment.l);
                s->kseq[s->n].comment.l = p->ks->comment.l;
                s->kseq[s->n].comment.s[p->ks->comment.l] = '\0';
            }
            
            //Copy sequence
            MALLOC(s->kseq[s->n].seq.s, p->ks->seq.l);
            memcpy(s->kseq[s->n].seq.s, p->ks->seq.s, p->ks->seq.l);
            s->kseq[s->n].seq.l = p->ks->seq.l;
            s->kseq[s->n].seq.s[p->ks->seq.l] = '\0';

            s->sum_len += p->ks->seq.l;

            //Copy quality
            MALLOC(s->kseq[s->n].qual.s, p->ks->qual.l);
            memcpy(s->kseq[s->n].qual.s, p->ks->qual.s, p->ks->qual.l);
            s->kseq[s->n].qual.l = p->ks->qual.l;
            s->kseq[s->n].qual.s[p->ks->qual.l] = '\0';

            //Increment to store next SE read
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
            fprintf(stderr, "\rRead: SE %i sequences (total %i)", s->n, s->p->total);
            return s;
        }
    } else if (step == 1) { // step 2: trim sequences
        stepdat_t *s = (stepdat_t*)in;

        kt_for(p->n_thread, worker_for, s, s->n);

        return s;
    } else if (step == 2) { // step 3: Output files
        stepdat_t *s = (stepdat_t*)in;

        for (int i = 0; i < s->n; i++) {
            
            if (s->p->debug) {
                printf("P1cut: %d,%d\n", s->pcut[i]->five_prime_cut, s->pcut[i]->three_prime_cut);
            }

            /* if sequence quality and length pass filter then output record, else discard */
            if (s->pcut[i]->three_prime_cut >= 0) {
                if (!s->p->gzip_output) {
                    /* This print statement prints out the sequence string starting from the 5' cut */
                    /* and then only prints out to the 3' cut, however, we need to adjust the 3' cut */
                    /* by subtracting the 5' cut because the 3' cut was calculated on the original sequence */

                    print_record(s->p->outfile, &s->kseq[i], s->pcut[i]);
                } else {
                    print_record_gzip(s->p->outfile_gzip, &s->kseq[i], s->pcut[i]);
                }

                s->p->kept++;
            } else {
                s->p->discard++;
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

int single_main(int argc, char *argv[]) {

    gzFile se = NULL;
    kseq_t *fqrec;
    int l;
    FILE *outfile = NULL;
    gzFile outfile_gzip = NULL;
    int debug = 0;
    int optc;
    extern char *optarg;
    int qualtype = -1;
    cutsites *p1cut;
    char *outfn = NULL;
    char *infn = NULL;
    int kept = 0;
    int discard = 0;
    int quiet = 0;
    int no_fiveprime = 0;
    int trunc_n = 0;
    int gzip_output = 0;
    int total=0;
    int n_thread = 1; // worker threads
    int p_thread = 3; // pipeline threads
    int block_size = 20000000; // Total size of the sequence loaded per thread

    while (1) {
        int option_index = 0;
        optc = getopt_long(argc, argv, "df:t:u:U:b:o:q:l:zxng", single_long_options, &option_index);

        if (optc == -1)
            break;

        switch (optc) {
            if (single_long_options[option_index].flag != 0)
                break;

        case 'f':
            infn = (char *) malloc(strlen(optarg) + 1);
            strcpy(infn, optarg);
            break;

        case 't':
            if (!strcmp(optarg, "illumina"))
                qualtype = ILLUMINA;
            else if (!strcmp(optarg, "solexa"))
                qualtype = SOLEXA;
            else if (!strcmp(optarg, "sanger"))
                qualtype = SANGER;
            else {
                fprintf(stderr, "Error: Quality type '%s' is not a valid type.\n", optarg);
                return EXIT_FAILURE;
            }
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

        case 'o':
            outfn = (char *) malloc(strlen(optarg) + 1);
            strcpy(outfn, optarg);
            break;

        case 'q':
            single_qual_threshold = atoi(optarg);
            if (single_qual_threshold < 0) {
                fprintf(stderr, "Quality threshold must be >= 0\n");
                return EXIT_FAILURE;
            }
            break;

        case 'l':
            single_length_threshold = atoi(optarg);
            if (single_length_threshold < 0) {
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

        case_GETOPT_HELP_CHAR(single_usage)
        case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

        case '?':
            single_usage(EXIT_FAILURE, NULL);
            break;

        default:
            single_usage(EXIT_FAILURE, NULL);
            break;
        }
    }


    if (qualtype == -1 || !infn || !outfn) {
        single_usage(EXIT_FAILURE, "****Error: Must have quality type, input file, and output file.");
    }

    if (!strcmp(infn, outfn)) {
        fprintf(stderr, "****Error: Input file is same as output file.\n\n");
        return EXIT_FAILURE;
    }

    se = gzopen(infn, "r");
    if (!se) {
        fprintf(stderr, "****Error: Could not open input file '%s'.\n\n", infn);
        return EXIT_FAILURE;
    }

    if (!gzip_output) {
        outfile = fopen(outfn, "w");
        if (!outfile) {
            fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn);
            return EXIT_FAILURE;
        }
    } else {
        outfile_gzip = gzopen(outfn, "w");
        if (!outfile_gzip) {
            fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn);
            return EXIT_FAILURE;
        }
    }

    pldat_t pl;
    pl.ks = kseq_init(se);
    pl.n_thread = n_thread;
    pl.block_len = block_size;
    pl.qualtype = qualtype;
    pl.single_length_threshold = single_length_threshold;
    pl.single_qual_threshold = single_qual_threshold;
    pl.no_fiveprime = no_fiveprime;
    pl.trunc_n = trunc_n;
    pl.debug = debug;

    pl.gzip_output = gzip_output;
    pl.outfile = outfile;
    pl.outfile_gzip = outfile_gzip;

    //Stats
    pl.total = 0;
    pl.kept = 0;
    pl.discard = 0;

    kt_pipeline(p_thread, worker_pipeline, &pl, 3);

    kseq_destroy(pl.ks);

    if (!quiet) fprintf(stdout, "\nSE input file: %s\n\nTotal FastQ records: %d\nFastQ records kept: %d\nFastQ records discarded: %d\n\n", infn, pl.total, pl.kept, pl.discard);

    kseq_destroy(fqrec);
    gzclose(se);
    if (!gzip_output) fclose(outfile);
    else gzclose(outfile_gzip);

    return EXIT_SUCCESS;
}
