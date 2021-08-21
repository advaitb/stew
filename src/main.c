#include <stdio.h>
#include <omp.h>
#include <log.h>
#include <ketopt.h>
#include <kseq.h>
#include <ascii.h>


#include <hll.h>
#include <hll_private.h>
#include <zlib.h>

#define FILE_LOG_LEVEL 0
#define CONSOLE_LOG_LEVEL 2
#define LOG_FILE "stew.log"
#define _VERSION_ "0.1.0"

KSEQ_INIT(gzFile , gzread);

// longopts params
static ko_longopt_t main_longopts[] = {
        { "threads", ko_required_argument, 't' },
        { "platters", ko_required_argument, 'p' },
        { "cups", ko_required_argument, 'c' },
        { "kmers", ko_required_argument, 'k' },
        { "threshold", ko_required_argument, 'x' },
        { "momentum", ko_required_argument, 'm' },
        { "version", ko_no_argument, 'v' },
        { NULL, 0, 0 }
};

// log params and setup
void log_setup(const char* fname, int f_log_lvl, int c_log_lvl)
{
    log_set_level(c_log_lvl);
    FILE *lfp = fopen(fname,"w+");
    log_add_fp(lfp,f_log_lvl);
}

// write the output
void stew_write(kseq_t *seq, bool is_fastq, FILE *fp_o)
{
    if (is_fastq)
    {
        fprintf(fp_o, "@%s %s\n", seq->name.s, seq->comment.s);
        fprintf(fp_o, "%s\n", seq->seq.s);
        fprintf(fp_o, "+\n");
        fprintf(fp_o, "%s\n", seq->qual.s);
    }
    else
    {
        fprintf(fp_o, ">%s\n", seq->name.s);
        fprintf(fp_o, "%s\n", seq->seq.s);
    }
}


int main(int argc, char *argv[])
{

    // usage string
    char* usage = ascii_art
                  "\nStew v"
                  _VERSION_
                  "\nDeveloper: Advait Balaji (advait@rice.edu)\n"
                  "\n"
                  "***Diversify and subsample reads into a stew!***\n"
                  "\n\nUsage stew [Subcommand] [options] [input.*|input1.*|input2.*] "
                  "[out.*...]\n"
                  "Subcommands:\n"
                  "\tS - Single end read mode\n"
                  "\tP - Paired end read mode\n"
                  "\n"
                  "Main options:\n"
                  "\t-t (--threads) - Number of threads [Default: 1]\n"
                  "\t-p (--platters) - Number of platters (arrays) of HLL structures "
                  "[Default: 10, Max: 50]\n"
                  "\t-c (--cups) - Number of cups (bits) in each HLL platter (array) "
                  "[Default: 8, Min: 4, Max: 16]\n"
                  "\t-k (--kmers) - Kmer size [Default: 23, Max: 100]\n"
                  "\t-x (--threshold) - Threshold for similarity "
                  "[Default: 0.5, Min: 0 (least selective), Max: 1 (most selective)]\n"
                  "\t-m (--momentum) - Momentum applied to boost score (Useful in bigger "
                  "datasets) [Default: 0.000001, Max 0.001]\n"
                  "\t-v (--version) - Print version\n"
                  "\n"
                  "Subcommand postitional options:\n"
                  "\tS\n"
                  "\t\tS [input.*] [out.*]\n"
                  "\tP\n"
                  "\t\tP [input1.*] [input2.*] [out1.*] [out2.*]\n"
                  "\n";

    // set logging
    log_setup(LOG_FILE, FILE_LOG_LEVEL, CONSOLE_LOG_LEVEL);

    // argument parsing
    ketopt_t om = KETOPT_INIT, os = KETOPT_INIT;
    int i, j, c;
    char *sf[2], *pf[4], *params;
    int t = 1, p = 10, cps = 16, k = 23;
    float x = 0.5, m = 0.001;
    while ((c = ketopt(&om, argc, argv, 1, "t:p:k:c:x:m:v", main_longopts)) >= 0)
    {
        if (c == 't')
        {
            t = om.arg ?  atoi(om.arg) : 1;
        }
        else if (c == 'p')
        {
            p = om.arg ? atoi(om.arg) : 10;
        }
        else if (c == 'c')
        {
            cps = om.arg ? atoi(om.arg) : 8;
        }
        else if (c == 'k')
        {
            k = om.arg ? atoi(om.arg) : 23;
        }
        else if (c == 'x')
        {
            x = om.arg ? atof(om.arg) : 0.5;
        }
        else if (c == 'm')
        {
            m  = om.arg ? atof(om.arg) : 0.000001;
        }
        else if (c == 'v')
        {
            log_info("stew version: %s", _VERSION_);
            return 0;
        }
    }
    // check args
    if (om.ind == argc)
    {
        log_error("No subcommand provided!");
        log_debug(usage);
        return 1;
    }

    // check subcommand
    char *sub = argv[om.ind];
    if (strcmp(sub,"S") && strcmp(sub,"P"))
    {
        log_error("No subcommand provided!");
        log_debug(usage);
        return 1;
    }


    // reset args
    t = (t > omp_get_max_threads()  || t <= 0) ?
    log_warn("Threads out of bounds, setting threads to maximum available"),
    omp_get_max_threads() : t;
    p = p > 50 ? log_warn("Platters more than allowed, setting to 50"), 50 : p;
    cps = (cps < 4 || cps > 16) ? log_warn("Cups out of bounds, setting to 8"), 8 : cps;
    k = ( k > 100 || k <= 0) ? log_warn("Kmers out of bound, setting to 100"), 100 : k;
    x = (x > 1 || x < 0) ?
            log_warn("Threshold out of bounds, all sequences will be preserved!"), 0 : x;
    m = (m > 0.001) ? log_warn("Momentum out of bounds, setting to 0.001"), 0.001 : m;

    // set threads
    //omp_set_num_threads(t);

    log_info(ascii_art);
    log_info("Preparing stew!...");

    if (!strcmp(sub,"S"))
    {
        if (argc - (os.ind + om.ind) != 2)
        {
            log_error("Positional arguments should (only) include "
                      " input.* and output.*");
            log_debug(usage);
            return 1;
        }

        for (i = os.ind + om.ind, j = 0; i < argc; ++i, ++j)
        {
            sf[j] = (char *)malloc(strlen(argv[i]) + 1);
            strcpy(sf[j],argv[i]);
        }
        params = "Stew params:\n"
                 "\tMode: %s\n"
                 "\tThreads: %d\n"
                 "\tPlatters: %d\n"
                 "\tCups: %d\n"
                 "\tKmers: %d\n"
                 "\tInput: %s\n"
                 "\tOutput: %s";
        log_debug(params,sub, t, p, cps, k, sf[0], sf[1]);
    }
    else
    {
        if (argc - (os.ind + om.ind) != 4)
        {
            log_error("Positional arguments should (only) include "
                      "input1.* input2.* out1.* out2.*");
            log_debug(usage);
            return 1;
        }
        for (i = os.ind + om.ind, j = 0; i < argc; ++i, ++j)
        {
            pf[j] = (char *)malloc(strlen(argv[i]) + 1);
            strcpy(pf[j],argv[i]);
        }
        params = "Stew params:\n"
                 "\tMode: %s\n"
                 "\tThreads: %d\n"
                 "\tPlatters: %d\n"
                 "\tCups: %d\n"
                 "\tKmers: %d\n"
                 "\tInput1: %s\n"
                 "\tInput2: %s\n"
                 "\tOutput1: %s\n"
                 "\tOutput1: %s";
        log_debug(params,sub, t, p, cps, k, pf[0], pf[1], pf[2], pf[3]);
    }

    log_info("Ingredients check completed! Firing up the stove!...");

    // create hll arrays
    hll_t *hll[p];

    for (int i = 0; i < p; i++)
    {
        hll[i] = hll_create(cps);
    }

    log_info("Cups and Platters are ready!...");

    int count = 1, sel_count = 0;
    if (!strcmp(sub,"S"))
    {
        int l;
        gzFile *sfp = gzopen(sf[0], "r");
        FILE *sfp_o = fopen(sf[1],"w+");
        if (!sfp)
        {
            log_error("Couldn't open file");
            return 1;
        }

        int *prev_cnt = (int *)calloc(p, sizeof(int));
        int *curr_cnt = (int *)calloc(p, sizeof(int));
        int *avg = (int *)calloc(p, sizeof(int));


        kseq_t *seq = kseq_init(sfp);
        int max_nk = 0, corr = 0, diff_cnt = 0;
        float corr_cnt = 0.0;

        log_debug("Reading the recipe!...");

        while ((l = kseq_read(seq)) >= 0)
        {
            bool is_fastq = false;
            if (seq->qual.l && seq->comment.l) is_fastq = true;
            float score = 0.0;
            long sum_curr = 0;

            int _nk = (seq->seq.l - k + 1) / p; // kmer per bucket
            int _throw = (seq->seq.l - k + 1) % p; // extra kmers
            int _effk = seq->seq.l - k + 1 - _throw; // effective kmers


            if (_nk < max_nk) // is this the largest number of kmers?
            {
                corr = max_nk - _nk; // no? apply corrections
            }
            else
            {
                max_nk = _nk; // yes? assign max
                corr = 0; // no correction needed
            }

            int _p = -1;
            char *kmer = (char *)malloc(sizeof(char)*k);


            for (int _s = 0; _s < _effk; _s++) // kmerize and add to HLL
            {
                if (!(_s % _nk)) _p++;
                strncpy(kmer, seq->seq.s+_s, k);
                hll_add(hll[_p], kmer, k);
            }

            //#pragma omp parallel for num_threads(t)
            for (int i = 0; i < p; i++) // estimate the count and calculate the uniqueness score
            {
                hll_estimate_t estimate;
                hll_get_estimate(hll[i], &estimate);
                curr_cnt[i] = estimate.estimate;
            }

            // split loop - may lead to lesser cache misses

            //#pragma omp parallel for num_threads(t)
            for (int i = 0; i < p; i++)
            {
                diff_cnt = curr_cnt[i] - prev_cnt[i];
                corr_cnt  = diff_cnt + (1-x)*(corr+m*avg[i]+m*count); // corrections added to unique kmers
                avg[i] = (avg[i]*(count-1) + corr_cnt) / count;
                score += (corr_cnt / _nk)*curr_cnt[i];
                sum_curr += curr_cnt[i];
                prev_cnt[i] = curr_cnt[i];
            }

            score /= sum_curr; // normalize

            if (score > x) // yup! we need this sequence.
            {
                sel_count++;
                stew_write(seq,is_fastq,sfp_o); // write this
            }
            count++;
        }
        log_debug("Finished processing the recipe!...");

        // clean up
        kseq_destroy(seq);
        gzclose(sfp);
        fclose(sfp_o);

        free(prev_cnt);
        free(curr_cnt);
        free(avg);
    }

    // release HLL allocs
    for (int i = 0; i < p; i++)
    {
            hll_release(hll[i]);
    }

    // that's all folks!
    log_info("Selected %d out of %d sequences!..", sel_count, --count);
    log_info("Piping hot Stew served! Bon appetit!...");

    return 0;
}
