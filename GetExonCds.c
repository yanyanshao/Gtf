/******************************************************************
    > File Name: GetExonCds.c
    >  Author: yys
    >  mail: shayy0919@163.com
    >  Created Time: 2019年04月21日 星期日 14时00分52秒
******************************************************************/

#include <getopt.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <regex.h>
#include <stdbool.h>

#define PATH_MAX 256
#define GENE_MAX 128

enum SEQTYPE { CDS = 0, EXON = 1 };

typedef struct __exon_id exon_id;
struct __exon_id{
    uint64_t start;
    uint64_t end;
    exon_id *link;
};

typedef struct __cds_id cds_id;
struct __cds_id{
    uint64_t start;
    uint64_t end;
    cds_id *link;
};

typedef struct __gene_e_t {
    char gene[0x80];
    int exon_n;
    exon_id *exon;
}gene_e_t;

typedef struct __gene_c_t {
    char gene[0x80];
    int cds_n;
    cds_id *cds;
}gene_c_t;

typedef struct __chr_t {
    char chr[0x80];
    int gene_n;
    gene_e_t *gene_e;
}chr_t;

typedef struct __gtf_t {
    int chrnum;
    chr_t *chr;
}gtf_t;

typedef struct __arg_t {
    int help;
    int region; // CDS(cds) or EXON(exon)
    char gtfile[PATH_MAX];
    char gene[PATH_MAX];
    char outdir[PATH_MAX];
} arg_t;

#define err_open(_fp, _fn, _mode) do {\
    _fp = fopen(_fn, _mode); \
    if (!_fp) { \
        fprintf(stderr, "\nerr: failed to open %s!\n", _fn); exit(-1); \
    }\
}while(0)

#define err_realloc(_p, _n, _type) do { \
    _type *tem = (_type *)realloc((_p), (_n)*sizeof(_type)); \
    if (!tem) { \
        fprintf(stderr, "\nerr: failed to realloced memory!\n"); exit(-1); \
    } (_p) = tem; \
}while(0)

#define err_calloc(_p, _n, _type) do { \
    _type *tem = (_type *)calloc(_n,sizeof(_type)); \
    if (!tem) { \
        fprintf(stderr, "\nerr: failed to calloc memory!\n"); exit(-1); \
    } (_p) = tem; \
}while(0)

void Usage(void)
{
    char *usage =
        "\nUsage: GetExonCds [options]\n"
        "Version: 1.0\n"
        "\n"
        "Options:\n"
        "       -h|--help        print help infomation\n"
        "       -r|--region      [required] the region type [exon|cds]\n"
        "       -f|--gtfile      [required] input gtf file [.gtf]\n"
        "       -g|--gene        [required] required gene name [\"chr7-EGFR, chr17-TP53\"]\n"
        "       -o|--outdir      [required] output directory for trimed fastq file [dir]\n\n";

    fprintf(stderr, "%s", usage);
    exit(-1);
}

static const struct option long_options[] =
{
    { "help", no_argument, NULL, 'h' },
    { "region", required_argument, NULL, 'r' },
    { "gtfile", required_argument, NULL, 'f' },
    { "gene", required_argument, NULL, 'g' },
    { "outdir", required_argument, NULL, 'o' },
    { NULL, 0, NULL, 0 }
};

static void ArgInit(arg_t *Arg)
{
    Arg->region = -1;
}

arg_t *ParseOpt( int argc, char **argv )
{
    int opt =0, opterr =0;
    arg_t *Arg;

    err_calloc(Arg, 1, arg_t);
    ArgInit(Arg);
    while ( (opt = getopt_long(argc, argv, "r:f:g:o:h", long_options, NULL)) != -1 )
    {
        switch (opt) {
            case 'h': Arg->help = 1; break;
            case 'r': if (!strcmp(optarg, "cds")) Arg->region = CDS;
                      else if (!strcmp(optarg, "exon")) Arg->region = EXON;
                      else Arg->region = -1; break;
            case 'f': strcpy(Arg->gtfile, optarg); break;
            case 'g': strcpy(Arg->gene, optarg); break;
            case 'o': strcpy(Arg->outdir, optarg); break;
            case '?': fprintf(stderr, \
                            "[Err::%s::%d]  Option error occour!.\n", __func__, __LINE__);
                      Arg->help = 1;
        }
    }
    if (!Arg->gtfile[0] || !Arg->gene[0] || (Arg->region == -1) || !Arg->outdir[0]) {
        fprintf(stderr, \
            "[Err::%s::%d]  Please give the [requied] parmeters!\n", __func__, __LINE__);
        Arg->help = 1;
    }

    return Arg;
}

void substr(char *des, char *src, uint64_t start, uint64_t end)
{
    char *p = src+start+11;

    int nlen = start;
}

int InsertNodeExon(exon_id **linkp, uint64_t start, uint64_t end)
{
    exon_id *current;
    exon_id *new;

    while((current = *linkp) != NULL) {
        //next node address
        linkp = &current->link;
    }

    new = (exon_id *)calloc(1,sizeof(exon_id));
    if(new == NULL) return false;
    new->start = start;
    new->end = end;

    new->link = current;
    *linkp = new;
    return true;
}

void PrintList(exon_id *S, char *chr, char *gene, char *outdir)
{
    exon_id *p;
    p = S;

    FILE *fp;
    char outfile[PATH_MAX]={'\0'};
    
    strcat(outfile, outdir); 
    strcat(outfile, gene);
    strcat(outfile, ".bed");
    fp = fopen(outfile, "w");
 
    while(p != NULL) {
        fprintf(fp,"%s\t%s\t%d\t%d\n", chr, gene, p->start, p->end);
        p = p->link;
    }
}

int CalCommaCount(char *genelist, char (*buf)[PATH_MAX])
{
    int n = 0;
    char *token;
    //char buf[GENE_MAX][PATH_MAX] = {'\0'};

    token = strtok(genelist, ",");
    while(token != NULL) {
        strcpy(buf[n], token);n++;
        token = strtok(NULL, ",");
    }
    return n;
}

gtf_t *GetExonCds(char *filename)
{
    FILE *gtfp;

    int status = 0;
    regex_t reg;
    regmatch_t pmatch[1];
    int flag = REG_EXTENDED;
    const size_t nmatch = 1;
    const char *pattern ="gene_name \"\\w*\\W?\\w*?\\W?\\w*?\\W?\\w?\";";

    uint64_t start, end;
    char buf[0x1000];
    char chr[0x80], info[0x80], suffix[0x1000];

    gtf_t *gtf;
    chr_t *chromsome;
    gene_e_t *gene_e;
    char curgene[0x80], postgene[0x80]="\0";
    char curchr[0x80], postchr[0x80]="\0";

    err_open(gtfp, filename, "r");
    gtf = (gtf_t *)calloc(1, sizeof(gtf_t));
    
    while(fgets(buf, 0x1000, gtfp)) {
        if(buf[0] == '#') continue;
        sscanf(buf, "%s%*s%s%ld%ld", curchr, info, &start, &end);
        regcomp(&reg, pattern, flag);
        status = regexec(&reg, buf, nmatch, pmatch, 0);
        if(status == REG_NOMATCH){
            printf("no match\n");
        }
        else if(status == 0){
            if (!strcmp(info, "exon")) {
                strncpy(curgene,buf+pmatch[0].rm_so+11,pmatch[0].rm_eo-pmatch[0].rm_so-13);
                curgene[pmatch[0].rm_eo-pmatch[0].rm_so-13]='\0';
         
                if (strcmp(curchr, postchr)) {
                    if (gtf->chrnum % 0x20 == 0) {
                        err_realloc(gtf->chr, gtf->chrnum+0x20, chr_t);
                    }
                    chromsome = &gtf->chr[gtf->chrnum++];
                    memset(chromsome, 0, sizeof(chr_t));
                    strcpy(chromsome->chr, curchr);
                    strcpy(postchr, curchr);
                }

                if (strcmp(curgene, postgene)) {
                    if (chromsome->gene_n % 0x20 == 0) {
                        err_realloc(chromsome->gene_e, chromsome->gene_n+0x20, gene_e_t);
                    }
                    gene_e = &chromsome->gene_e[chromsome->gene_n++];
                    memset(gene_e, 0, sizeof(gene_e_t));
                    strcpy(gene_e->gene, curgene);
                    strcpy(postgene, curgene);
                }
                gene_e->exon_n++;
                InsertNodeExon(&gene_e->exon, start, end);
            }
            regfree(&reg);
        }
    }
    return gtf;
}

int main(int argc, char **argv)
{
    arg_t *args = ParseOpt(argc, argv);
    if (args->help) Usage();

    int gene_count;
    char ch[0x80],ge[0x80];
    char buf[GENE_MAX][PATH_MAX] = {'\0'};

    gtf_t *gtf;
    chr_t *chrom;
    gene_e_t *gene;

    gtf = GetExonCds(args->gtfile);
    gene_count = CalCommaCount(args->gene,buf);

    for(int i=0; i<gene_count; i++) {
        strcpy(ch,strtok(buf[i],"-"));
        strcpy(ge,strtok(NULL,"-"));
        for (int chr_n=0; chr_n<gtf->chrnum; chr_n++) {
            chrom = &gtf->chr[chr_n];
            if (!strcmp(chrom->chr, ch)) {
                for (int gene_n=0; gene_n<chrom->gene_n; gene_n++) {
                    gene = &chrom->gene_e[gene_n];
                    if (!strcmp(gene->gene,ge)) {
                        if (args->region)
                            PrintList(gene->exon, chrom->chr, gene->gene, args->outdir);
                    }
                }
            }
        }
    }
}
