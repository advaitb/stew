#include <stdio.h>
#include <omp.h>
#include <log.h>
#include <ketopt.h>
#include <ascii.h>

#define FILE_LOG_LEVEL 0
#define CONSOLE_LOG_LEVEL 2
#define LOG_FILE "stew.log"
#define _VERSION_ "0.1.0"


static ko_longopt_t main_longopts[] = {
        { "threads", ko_optional_argument, 't' },
        { "platters", ko_optional_argument, 'p' },
        { "kmers", ko_optional_argument, 'k' },
        { "version", ko_optional_argument, 'v'},
        { NULL, 0, 0 }
};


void log_setup(const char* fname, int f_log_lvl, int c_log_lvl)
{
    log_set_level(c_log_lvl);
    FILE *lfp = fopen(fname,"w+");
    log_add_fp(lfp,f_log_lvl);
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
                  "\n\nUsage stew [Subcommand] [options] [read.*|read1.*|read2.*] [out.*...]\n"
                  "Subcommands:\n"
                  "\tS - Single end read mode\n"
                  "\tP - Paired end read mode\n"
                  "\n"
                  "Main options:\n"
                  "\t-t (--threads) - Number of threads [Default: 1]\n"
                  "\t-p (--platters) - Number of platters (arrays) of HLL structures "
                  "[Default: 10]\n"
                  "\t-k (--kmers) - Kmer size [Default: 23]\n"
                  "\t-v (--version) - Print version\n"
                  "\n"
                  "Subcommand postitional options:\n"
                  "\tS\n"
                  "\t\tS [read.*] [out.*]\n"
                  "\tP\n"
                  "\t\tP [read1.*] [read2.*] [out1.*] [out2.*]\n"
                  "\n";

    // set logging
    log_setup(LOG_FILE, FILE_LOG_LEVEL, CONSOLE_LOG_LEVEL);

    // get threads in system
    //int t = omp_get_max_threads();
    //log_trace("Num threads available %d",t);

    // argument parsing
    ketopt_t om = KETOPT_INIT, os = KETOPT_INIT;
    int i, c;
    int t,p,k;
    while ((c = ketopt(&om, argc, argv, 0, "tpkv:", main_longopts)) >= 0)
    {
        if (c == 't')
        {
            t = om.arg ? om.arg : 1;
        }
        else if (c == 'p')
        {
            p = om.arg ? om.arg : 10;
        }
        else if (c == 'k')
        {
            k =  om.arg  ? om.arg : 23;
        }
        else if (c == 'v')
        {
            log_info("stew version: %s", _VERSION_);
            return 0;
        }
        else log_warn("Unrecognized argument %c skipped", om.opt);
    }
    if (om.ind == argc)
    {
        log_error("No subcommand provided!");
        log_info(usage);
        return 1;
    }

    char *sub = argv[om.ind];
    if (strcmp(sub,"S") && strcmp(sub,"P"))
    {
        log_error("Wrong subcommand provided!");
        log_info(usage);
        return 1;
    }

    fprintf("%u", argc-(os.ind+om.ind));
    for (i = os.ind + om.ind; i < argc; ++i)
        fprintf("%s\n", argv[i]);


    // check if subcommand valid
    return 0;
}
