library(ShortRead)

samtools <- '~/software/samtools/bin/samtools'
bwa2 <- '~/software/bwa-mem2/bwa-mem2'

system('cat WZ592V_1_gRVA_1.fastq WZ592V_2_gRVA_2.fastq > gRVA.fastq')
system('cat WZ592V_3_gRVB_1.fastq WZ592V_4_gRVB_2.fastq > gRVB.fastq')
system('cat WZ592V_5_gRVC_1.fastq WZ592V_6_gRVC_2.fastq > gRVC.fastq')

alignReads <- function(f){
  o <- readFastq(f)
  o <- trimTails(o, 2, '+', 5)
  o <- o[width(o) >= 30]
  writeFastq(o, 'f.fastq', mode = 'w')
  system(paste0(bwa2, ' mem -t 10 ref f.fastq > f.sam'))
  system(paste0(samtools, ' view -S -b f.sam > f.bam'))
  system(paste0(samtools, ' view -b -f 4 f.bam > f.unmapped.bam'))
  system(paste0(samtools, ' view -b -F 4 f.bam > f.mapped.bam'))
  system(paste0(samtools, ' sort -o f.mapped.sorted.bam f.mapped.bam'))
  system(paste0(samtools, ' index f.mapped.sorted.bam'))
  system(paste0(samtools, ' view f.mapped.sorted.bam | cut -f1 | sort | uniq > f.mapped.sorted.readIDs'))
  message((length(o) - length(readLines('f.mapped.sorted.readIDs'))), ' trimmed reads did not alingn to the expected proviral sequence.')
  message(sprintf('%.2f%%', (length(readLines('f.mapped.sorted.readIDs')) / length(o)) * 100), ' of ', length(readLines('f.mapped.sorted.readIDs')), ' trimmed reads aligned.')
}

alignReads('gRVA.fastq')
system('mv f.mapped.sorted.bam gRVA.mapped.sorted.bam')
system('mv f.mapped.sorted.bam.bai gRVA.mapped.sorted.bam.bai')
system('rm f.*')

alignReads('gRVB.fastq')
system('mv f.mapped.sorted.bam gRVB.mapped.sorted.bam')
system('mv f.mapped.sorted.bam.bai gRVB.mapped.sorted.bam.bai')
system('rm f.*')

alignReads('gRVC.fastq')
system('mv f.mapped.sorted.bam gRVC.mapped.sorted.bam')
system('mv f.mapped.sorted.bam.bai gRVC.mapped.sorted.bam.bai')
system('rm f.*')
