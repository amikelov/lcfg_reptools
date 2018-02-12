# lcfg_reptools


mixcr_TRB is a simple R wrapper function for MIXCR (Milaboratory) for standard processing of multiple TRB repseq files (5'RACE derived)
 usage:
 
 mixcr_TRB(samples, outputFolder, new_meta_features, return_cloneset_list=T )
 
 Arguments:
      samples              character string, path to a tab-delimited metafile. Headers are required, columns id, R1, R2 are mandatory.
                           Can also take a data.frame as input.
                           
                          | id            | R1                   | R2                   | Any_optional_featute|
                          | ------------- |:--------------------:| --------------------:|---------------------|
                          | Sample_X      | /data/SX_R1.fastq.gz | /data/SX_R2.fastq.gz |        Healthy      |
                          | Sample_Y      | /data/SX_R1.fastq.gz | /data/SX_R2.fastq.gz |          AS         |
                          | Sample_5      |        ...           |         ...          |          ...        |
      
      
      outputFolder           path to outputFolder
      
      new_meta_features      character vector, containing colnames for features to be saved in new meta file. By default - id, paths to cloneset files                              and other optional features are used, while R1 and R2 are dropped.
      
      return_cloneset_list   boolean, if T all clonesets are returned as a list of data.tables. Else empty vector is returned
      
      
  Output:
  In outputFolder following subdirectories are created: alignments/, clonesets/, export/, corresponding to mixcr conventional steps
  Additionally, new metafile with ids and paths to clonesets is created.
      
  