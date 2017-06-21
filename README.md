# lcfg_reptools

repGrep is a simple tool for finding many-to-many T-B.cr repertoire intersections with approximate matches (hamming)

usage:

> shared.counts(cls_list,tgs,nsubs,mod)

cls_list - list of cloneset data.frames, derived from parse.folder (or parse.file) functions from tcR package

tgs - target cdr3 aminoacid sequences to look for in the cloneset list. If ommited the function takes as input all of aa cdr3 sequences from the cloneset list, so the operation of intersection inside the cloneset list is performed  

nsubs - maximum number of aa substitutions allowed by search. Default - 1

mod - may take two possible values - "samples" and "clones". In "samples" mode the function returns as an output the number of samples, in which corresponding cdr3 sequence is present with supplied (nsubs) number of substitutions. "clones" counts the number of individual clones, which are similar with the target sequences with supplied (nsubs) number of substitutions



> getCdrs(cls_list) - gets melted table of all cdr3 aminoacid sequences in the list of clonest data.frames.
