# dexseq

  pjdir=/nfs/med-bfx-activeprojects/Soleimanpour_RS1_weishwu_damki_NovaA-225
  alndir=/nfs/med-bfx-activeprojects/Soleimanpour_RS1_weishwu_damki_NovaA-225/analysis_20200318/alignment_results/04-rsem_star_align
  dexseqbin=/nfs/med-bfx-activeprojects/weishwu/common/Anaconda3/envs/dexseq/bin
  cd ${pjdir}/dexseq

  # prepare annotation
  # by default, genes with overlapping exons are merged. If disabled, these exons will be skipped.
  [ -f Mus_musculus.GRCm38.98.gff ] || python ${dexseqbin}/dexseq_prepare_annotation.py /nfs/med-bfx-common/ENSEMBL_references/Mus_musculus/GRCm38/Mus_musculus.GRCm38.98.gtf Mus_musculus.GRCm38.98.gff

  # count reads
  # PE: -p yes
  # sam needs to be sorted either by name or by coordinate. By default, it assumes sorted by name
  # by default it assumes strand-specific (first sequence pass is on the same strand as the gene). Otherwise use "-s no" or "-s reverse"
  # by default all reads with a lower quality than specified (with default -a 10) are skipped.
  smps="Sample_716-AS-18 Sample_716-AS-19 Sample_716-AS-20 Sample_716-AS-4 Sample_716-AS-5 Sample_716-AS-6"

  for f in ${smps}
  do
  [ -f ${f}.genome.sorted.bam ] || /nfs/med-bfx-activeprojects/weishwu/common/GATK/gatk-4.1.4.0/gatk --java-options "-Xmx60G" SortSam -I ${alndir}/${f}.genome.bam -O ${f}.genome.sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true --CREATE_MD5_FILE false --TMP_DIR /nfs/med-bfx-activeprojects/weishwu/temp 2>&1|tee >sortbam_${f}.log &
  done

  wait

  for f in ${smps}
  do
  [ -f ${f}.txt ] || python ${dexseqbin}/dexseq_count.py -p yes -r pos -f bam Mus_musculus.GRCm38.98.gff ${f}.genome.sorted.bam ${f}.txt 2>&1|tee >dexseqcount_${f}.log &
  done


  ## in R

  # create a sample table
  smps=read.table('../../samples.csv',header=T,sep=',')
  #smps=smps[smps$sample%in%c("Sample_716-AS-18","Sample_716-AS-19","Sample_716-AS-20","Sample_716-AS-4","Sample_716-AS-5","Sample_716-AS-6"),]smps=smps[smps$sample%in%c("Sample_716-AS-7","Sample_716-AS-8","Sample_716-AS-9","Sample_716-AS-4","Sample_716-AS-5","Sample_716-AS-6"),]
  sampleTable=smps[,-1]
  rownames(sampleTable)=smps[,1]

  # load data
  library( "DEXSeq" )

  countFiles = list.files(".", pattern=".txt", full.names=TRUE)
  flattenedFile = "../Mus_musculus.GRCm38.98.gff"

  # if there are two factors: design= ~ sample + exon + libType:exon + condition:exon
  dxd = DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design= ~ sample + exon + Condition:exon, flattenedfile=flattenedFile )
  dxd = estimateSizeFactors( dxd )
  dxd = estimateDispersions( dxd )
  dxd = testForDEU( dxd )
  dxd = estimateExonFoldChanges( dxd, fitExpToVar="Condition")save(dxd, file='dexseq_Soleimanpour')
  dxr1 = DEXSeqResults( dxd )save(dxr1, file='dexseq_Soleimanpour_res')
  plotDEXSeq( dxr1, "ENSMUSG00000027668", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, fitExpToVar='Condition') # Mfn1: ENSMUSG00000027668; Mfn2: ENSMUSG00000029020
wh = (dxr1$groupID=="ENSMUSG00000027668")


