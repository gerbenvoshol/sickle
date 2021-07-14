#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include "sickle.h"
#include "kseq.h"

__KS_BASIC(gzFile, 16384)
__KS_GETC(gzread, 16384)
__KS_GETUNTIL(gzread, 16384)

extern kseq_t *kseq_init(gzFile fd);
void kseq_destroy(kseq_t *ks);
int kseq_read(kseq_t *seq);

// pre is the refix to search in str, only requiring n matches (Note n should be shorter than strlen(pre))
int prefix(const char *str, const char *pre, int n)
{
    return strncmp(pre, str, n) == 0;
}

// Check from the name if the reads are from a two color system
// Note: This function resets/rewinds the input stream, so check BEFORE when no sequences are read yet
int check_two_color_system(kseq_t *read)
{
    int retval = 0;
    int ret = 0;

    // Get read
    if ((ret = kseq_read(read)) < 0) {
        switch (ret) {
            case -1:
                fprintf(stderr, "Unexpected end-of-file\n");
                break;
            case -2:
                fprintf(stderr, "Truncated quality string\n");
                break;
            case -3:
                fprintf(stderr, "Error while reading stream\n");
                break;
            default:
                break;
        }
    }

    // NEXTSEQ500, NEXTSEQ 550/550DX, NOVASEQ
    if(prefix(read->name.s, "NS", 2) || prefix(read->name.s, "NB", 2) || prefix(read->name.s, "NDX", 3) || prefix(read->name.s, "A0", 2)) {
        retval = 1;
    }

    // rewind the stream
    ks_rewind(read->f);

    return retval;
}

// Trim 3prime (tail) of sequence if it contains a poly X of min_size length
int trim_polyX(kseq_t *read, const char X, const int min_size)
{
    const int one_mismatch_per = 8; // Allow one mismatch for each 8 base pairs
    const int max_mismatch = 5;     // Allow maximum of 5 mismatches

    int slen = read->seq.l - 1;
    int mismatch = 0;
    int i, allowed_mismatch, firstX = slen;
    for(i = slen; i >= 0; i--) {
        // Found another X update position
        if(read->seq.s[i] == X) {
            firstX = i;
        // Not an X, add mismatch and check if we continue 
        } else {
        	mismatch++;
	        allowed_mismatch = (slen - i + 1) / one_mismatch_per;
	        if ((mismatch > max_mismatch) || (mismatch > allowed_mismatch && (slen - i - 1) >= min_size)) {
	            break;
	        }        
        }
    }

    // If the poly X tail is long enough, trim the sequence
    if ((slen - firstX) >= min_size) {
        read->seq.l = firstX + 1;
        read->qual.l = firstX + 1;
        return 1;
    }

    return 0;
}

int get_quality_num (char qualchar, int qualtype, kseq_t *fqrec, int pos) {
  /* 
     Return the adjusted quality, depending on quality type.

     Note that this uses the array in sickle.h, which *approximates*
     the SOLEXA (pre-1.3 pipeline) qualities as linear. This is
     inaccurate with low-quality bases.
  */

  int qual_value = (int) qualchar;

  if (qual_value < quality_constants[qualtype][Q_MIN] || qual_value > quality_constants[qualtype][Q_MAX]) {
	fprintf (stderr, "ERROR: Quality value (%d) does not fall within correct range for %s encoding.\n", qual_value, typenames[qualtype]);
	fprintf (stderr, "Range for %s encoding: %d-%d\n", typenames[qualtype], quality_constants[qualtype][Q_MIN], quality_constants[qualtype][Q_MAX]);
	fprintf (stderr, "FastQ record: %s\n", fqrec->name.s);
	fprintf (stderr, "Quality string: %s\n", fqrec->qual.s);
	fprintf (stderr, "Quality char: '%c'\n", qualchar);
	fprintf (stderr, "Quality position: %d\n", pos+1);
	exit(1);
  }

  return (qual_value - quality_constants[qualtype][Q_OFFSET]);
}


cutsites* sliding_window (kseq_t *fqrec, int qualtype, int length_threshold, int qual_threshold, int no_fiveprime, int trunc_n, int debug) {

	int window_size = (int) (0.1 * fqrec->seq.l);
	int i,j;
	int window_start=0;
	int window_total=0;
	int three_prime_cut = fqrec->seq.l;
	int five_prime_cut = 0;
	int found_five_prime = 0;
	double window_avg;
	cutsites* retvals;
    char *npos;

    retvals = malloc(sizeof(cutsites));

	/* discard if the length of the sequence is less than the length threshold */
    if (fqrec->seq.l < length_threshold) {
		retvals->three_prime_cut = -1;
		retvals->five_prime_cut = -1;
		return (retvals);
	}

	/* if the seq length is less then 10bp, */
	/* then make the window size the length of the seq */
	if (window_size == 0) window_size = fqrec->seq.l;

	for (i=0; i<window_size; i++) {
		window_total += get_quality_num (fqrec->qual.s[i], qualtype, fqrec, i);
	}

	for (i=0; i <= fqrec->qual.l - window_size; i++) {

		window_avg = (double)window_total / (double)window_size;

        if (debug) printf ("no_fiveprime: %d, found 5prime: %d, window_avg: %f\n", no_fiveprime, found_five_prime, window_avg);

		/* Finding the 5' cutoff */
		/* Find when the average quality in the window goes above the threshold starting from the 5' end */
		if (no_fiveprime == 0 && found_five_prime == 0 && window_avg >= qual_threshold) {

        if (debug) printf ("inside 5-prime cut\n");

			/* at what point in the window does the quality go above the threshold? */
			for (j=window_start; j<window_start+window_size; j++) {
				if (get_quality_num (fqrec->qual.s[j], qualtype, fqrec, j) >= qual_threshold) {
					five_prime_cut = j;
					break;
				}
			}

            if (debug) printf ("five_prime_cut: %d\n", five_prime_cut);

			found_five_prime = 1;
		}

		/* Finding the 3' cutoff */
		/* if the average quality in the window is less than the threshold */
		/* or if the window is the last window in the read */
		if ((window_avg < qual_threshold || 
			window_start+window_size >= fqrec->qual.l) && (found_five_prime == 1 || no_fiveprime)) {

			/* at what point in the window does the quality dip below the threshold? */
			for (j=window_start; j<window_start+window_size; j++) {
				if (get_quality_num (fqrec->qual.s[j], qualtype, fqrec, j) < qual_threshold) {
					three_prime_cut = j;
					break;
				}
			}

			break;
		}

		/* instead of sliding the window, subtract the first qual and add the next qual */
		window_total -= get_quality_num (fqrec->qual.s[window_start], qualtype, fqrec, window_start);
		if (window_start+window_size < fqrec->qual.l) {
			window_total += get_quality_num (fqrec->qual.s[window_start+window_size], qualtype, fqrec, window_start+window_size);
		}
		window_start++;
	}


    /* If truncate N option is selected, and sequence has Ns, then */
    /* change 3' cut site to be the base before the first N */
    if (trunc_n && ((npos = strstr(fqrec->seq.s, "N")) || (npos = strstr(fqrec->seq.s, "n")))) {
        three_prime_cut = npos - fqrec->seq.s;
    }

    /* if cutting length is less than threshold then return -1 for both */
    /* to indicate that the read should be discarded */
    /* Also, if you never find a five prime cut site, then discard whole read */
    if ((found_five_prime == 0 && !no_fiveprime) || (three_prime_cut - five_prime_cut < length_threshold)) {
        three_prime_cut = -1;
        five_prime_cut = -1;

        if (debug) printf("%s\n", fqrec->name.s);
    }

    if (debug) printf ("\n\n");

	retvals->three_prime_cut = three_prime_cut;
	retvals->five_prime_cut = five_prime_cut;
	return (retvals);
}
