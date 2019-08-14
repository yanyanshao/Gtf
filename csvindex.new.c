#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*! @typedef cidx_t
 @abstract structure for the chrom idex
 @field chr          chr name [char *]
 @field bnum         bin count [ uint32_t *]
 @field bstart       position for bin start [uint64_t *]
 @field offset       offset of the bin [uint64_t *]
*/
typedef struct {
    char chr[0x20];
    uint32_t bnum;
    uint32_t *bstart;
    uint64_t *offset;
} cidx_t;

/*! @typedef vidx_t
 @abstract structure for vcf index
 @field bsize        the size for each bin [uint32_t *]
 @field chrnum       number of chr [ uint32_t *]
 @field cidx         index for the chr [cidx_t *]
*/
typedef struct {
    uint32_t bsize;
    uint32_t chrnum;
    cidx_t *cidx;
} vidx_t;

/*! @typedef loc_t
 @abstract structure for the output of search
 @field n            count of snp line within the region [uint32_t *]
 @field offset       offset for the first snp line in the region [ uint64_t *]
*/
typedef struct {
    uint32_t n;
    uint64_t offset;
} loc_t;


#define err_realloc(_p, _n, _type) do { \
    _type *tem = (_type *)realloc((_p), (_n)*sizeof(_type)); \
    if (!tem) { \
        fprintf(stderr, "\nerr: failed to realloced memory!\n"); exit(-1); \
    } (_p) = tem; \
} while(0)

#define err_open(_fp, _fn, _mode) do {\
    _fp = fopen(_fn, _mode); \
    if (!_fp) { \
        fprintf(stderr, "\nerr: failed to open %s!\n", _fn); exit(-1); \
    } \
} while(0)

#define err_read(_ptr, _size, _n, _fp) do {\
    size_t ret = fread((_ptr), (_size), (_n), _fp);\
    if (ret != (_n)) {\
        fprintf(stderr, "\nerr: the index file is truncated!\n"); exit(-2);\
    }\
} while(0)

#define err_write(_ptr, _size, _n, _fp) do {\
    size_t ret = fwrite((_ptr), (_size), (_n), _fp);\
    if (ret != (_n)) {\
        fprintf(stderr, "\nerr: failed to write the index file!\n"); exit(-2);\
    }\
} while(0)

#define err_order(_pos, _check) do { \
    if ((_pos) < (_check)) { \
        fprintf(stderr, "\nerr: vcffile must be ordered by chr and pos!\n"); \
        exit(-3);\
    } \
    (_check) = (_pos);\
} while(0)

static vidx_t *IndexBuild(char *vcf, uint32_t bsize)
{
    FILE *fpvcf;
    vidx_t *vidx; cidx_t *cidx;
    char buf[0x1000], cur[0x20], pre[0x20] = "dummy";
    uint32_t pos, prepos, check=0;

    vidx = (vidx_t *)calloc(1, sizeof(vidx_t));
    vidx->bsize = bsize;
    err_open(fpvcf, vcf, "r");

    while (fgets(buf, 0x1000, fpvcf)) {
        if (buf[0] == '#') continue;
        sscanf(buf, "%s%d", cur, &pos);
        if (strcmp(cur, pre)) {
            if (!(vidx->chrnum % 0x20))
                err_realloc(vidx->cidx, vidx->chrnum+0x20, cidx_t);
            cidx = &vidx->cidx[vidx->chrnum++];
            memset(cidx, 0, sizeof(cidx_t));
            strcpy(cidx->chr, cur); strcpy(pre, cur); check = 0;
        } 
        err_order(pos, check); // check whther the pos are ordered
        if (!(cidx->bnum % 0x400)) {
            int n = cidx->bnum + 0x400;
            err_realloc(cidx->bstart, n, uint32_t);
            err_realloc(cidx->offset, n, uint64_t);
        }
        if (!cidx->bnum) {
            cidx->bstart[cidx->bnum] = prepos = pos - pos%vidx->bsize;
            cidx->offset[cidx->bnum++] = ftell(fpvcf) - strlen(buf); continue;
        }
        if (pos - prepos > vidx->bsize) {
            prepos += vidx->bsize; 
            cidx->bstart[cidx->bnum] = prepos;
            cidx->offset[cidx->bnum++] = ftell(fpvcf) - strlen(buf);
        }
    } fclose(fpvcf);
    
    return vidx;
}


static int IndexSave(vidx_t *vidx, char *vcfidx)
{
    FILE *fp;
    cidx_t *cidx;
    char buf[4] = "IDX\1";

    err_open(fp, vcfidx, "w");
    err_write(buf, 1, 4, fp); // magic info
    err_write(&vidx->bsize, 4, 1, fp);
    err_write(&vidx->chrnum, 4, 1, fp);

    for (int i=0; i < vidx->chrnum; ++i) {
        cidx = &vidx->cidx[i];
        err_write(cidx->chr, 1, 0x20, fp);
        err_write(&cidx->bnum, 4, 1, fp);
        err_write(cidx->bstart, 4, cidx->bnum, fp);
        err_write(cidx->offset, 8, cidx->bnum, fp);
    } fclose(fp);

    return 0;
}

static vidx_t *IndexLoad(char *vcfidx)
{
    FILE *fp;
    char buf[4];
    cidx_t *cidx;

    vidx_t *vidx = calloc(1, sizeof(vidx_t));
    err_open(fp, vcfidx, "r"); err_read(buf, 1, 4, fp);

    if (strncmp(buf, "IDX\1", 4)) {
        fprintf(stderr, "err: %s is not a valid index file\n", vcfidx); exit(-2);
    }
    err_read(&vidx->bsize, 4, 1, fp);
    err_read(&vidx->chrnum, 4, 1, fp);
    vidx->cidx = (cidx_t*)calloc(vidx->chrnum, sizeof(cidx_t));

    for (int i=0; i < vidx->chrnum; ++i) {
        cidx = &vidx->cidx[i];
        err_read(cidx->chr, 1, 0x20, fp);
        err_read(&cidx->bnum, 4, 1, fp);

        err_realloc(cidx->bstart, cidx->bnum, uint32_t);
        err_realloc(cidx->offset, cidx->bnum, uint64_t);

        err_read(cidx->bstart, 4, cidx->bnum, fp);
        err_read(cidx->offset, 8, cidx->bnum, fp);
    } fclose(fp);

    return vidx;
}

static void IndexDestroy(vidx_t *vidx)
{
    cidx_t *cidx;

    for (int i=0; i < vidx->chrnum; ++i) {
        cidx = &vidx->cidx[i];
        free(cidx->bstart); free(cidx->offset);
    }
    free(vidx->cidx); free(vidx);

    return ;
}

static cidx_t *GetChr(vidx_t *vidx, char *chr)
{
    cidx_t *cidx;

    for (int i=0; i < vidx->chrnum; ++i) {
        cidx = &vidx->cidx[i];
        if (!strcmp(cidx->chr, chr)) return cidx;
    }
    return NULL;
}

static loc_t IndexSearch(FILE *vcf, vidx_t *vidx, char *chr, int32_t start, int32_t end)
{
    cidx_t *cidx;
    uint32_t index, pos, flag=1;
    char buf[0x1000], cur[0x20];
    loc_t loc = {0,0};

    cidx = GetChr(vidx, chr);
    if (!cidx) {
        fprintf(stderr, "err: can't find chr %s in index file\n", chr); 
        return loc;
    }

    int rlen = start - cidx->bstart[0];
    index = rlen > 0 ? rlen/vidx->bsize : 0;
    loc.offset = cidx->offset[index];
    fseek(vcf, loc.offset, SEEK_SET);

    while (fgets(buf, 0x1000, vcf)) {
        sscanf(buf, "%s%d", cur, &pos);
        if (strcmp(cur, chr) || pos > end) break;
        if (pos < start) continue;
        if (flag) {loc.offset = ftell(vcf) - strlen(buf); flag = 0;} 
        ++loc.n;
    }

    return loc;
}

/*****************************
 *  1. Index build and save  *
 *  2. Index load            *
 *  3. Binary search         *
*****************************/

void Build(char *vcf, uint32_t bsize)
{
    vidx_t *vidx;
    clock_t t;
    char buf[0x1000];

    t = clock();
    strcpy(buf, vcf); strcat(buf, ".idx");
    fprintf(stderr, "[indexbuild] start to build index ... ");
    vidx = IndexBuild(vcf, bsize);
    fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    t = clock();
    fprintf(stderr, "[indexbuild] start to save index ... ");
    IndexSave(vidx, buf);
    fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

    IndexDestroy(vidx);

    return ;
}

void *Load(char *vcf)
{
    vidx_t *vidx;
    clock_t t;
    char buf[0x1000];
   
    t = clock();
    strcpy(buf, vcf); strcat(buf, ".idx");
    fprintf(stderr, "[indexload] start to load the index ... ");
    vidx = IndexLoad(buf);
    fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

    return (void*)vidx;
}

loc_t Search(char *vcf, void *vidx, char *chr, int32_t start, int32_t end)
{
    FILE *fp;
    loc_t loc = {0,0};

    if (start > end) {
        fprintf(stderr, "err: invalid parmeters (start > end)\n"); 
        return loc;
    }
    err_open(fp, vcf, "r");
    loc = IndexSearch(fp, (vidx_t*)vidx, chr, start, end); fclose(fp);

    return loc;
}

void Destory(void *vidx)
{
    return IndexDestroy((vidx_t*)vidx);
}


int main(void)
{
    char vcf[0x100] = "./yys.snp";
    //loc_t loc;
    Build(vcf, 1000);
    //void *idx = Load(vcf);
    //loc = Search(vcf, idx, "chr5", 11814, 12112);
    //printf("%ld\n%d\n", loc.offset, loc.n);
}
