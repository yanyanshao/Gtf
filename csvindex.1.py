#!/usr/bin/python2.7

import sys
import os.path
from ctypes import *
import pdb

class loc_t(Structure):
    ''' @Structure: the return value of method search
        @field:
            n[c_uint]: line count between the start and end
            offset[c_ulonglong]: the offset of the first snp line
    '''
    _fields_ = [ ("n", c_uint), ("offset", c_ulonglong) ]


class Libidx(object):
    def __init__(self, idxpath):
        ''' @Initialization: load the dynamic library 
                    and set the return type of their function
            @parameter
                idxpath: the path to dynamic library(*.so)
        '''
        self.idx = CDLL(idxpath)
        self.idx.Load.restype = c_void_p
        self.idx.Search.restype = loc_t

    def build(self, vcf, bsize=1000):
        ''' @func: build index for *.vcf file
            @return: None
            @parameter
                vcf[c_char_p]: path to *.vcf file (vcf must be sorted by chr and pos)
                bsize[c_uint]: bin size [default 1000]
        '''
        return self.idx.Build(vcf, bsize)

    def load(self, vcf):
        ''' @func: load the index(.vcf.idx) of *.cvf file
            @return: (c_void_p)vidx
            @parameter
                vcf[c_char_p]: path to *.vcf file (vcf must be sorted by chr and pos)
        '''
        return self.idx.Load(vcf)

    def search(self, vcf, vidx, chrom, start, end):
        ''' @func: search for the snp line between the start and end
            @return: (loc_t)loc
                loc.offset[c_ulonglong]: the offset of the first snp line
                loc.n[c_uint]: line count between the start and end
            @parameter
                vcf[c_char_p]: path to *.vcf file (vcf must be sorted by chr and pos)
                vidx[c_void_p]: index of *.vcf
        '''
        return self.idx.Search(vcf, vidx, chrom, start, end)

    def free(self, vidx):
        ''' @func: free the memory allocated for vidx
            @return: None
            @parameter
                vidx[c_void_p]: index of *.vcf file
        '''
        return self.idx.Destory(vidx)


def test():
    vcf = "/home/yys/Project/Lilab/Primer/PrimerMotify/Data/vcf/dbsnp.sort.vcf"
    Idx = Libidx("./libidx.so")

    if not os.path.exists(vcf + ".idx"):
        Idx.build(vcf, 1000)
    pdb.set_trace()
    ll = [("chr5", 11814, 12112), ("chr2", 11814, 12112), ("chr3", 11814, 12112)]
    idx = Idx.load("/home/yys/Project/Lilab/Primer/PrimerMotify/Data/vcf/dbsnp.sort.vcf")
    
    for e in ll:
        loc = Idx.search(vcf, idx, e[0], e[1], e[2])
        print "%d\t%d" %(loc.n, loc.offset)

    Idx.free(idx)
    

if __name__ == '__main__':
    test()
