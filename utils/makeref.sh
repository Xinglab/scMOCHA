

# cellranger mkref --config=/home/liuc9/data/refdata/mgatk_index/reference.json


# UPDATE GTF file by removing contigs that not in fasta then makeref again,

#filter chromosome




cellranger mkref --genome=Human \
  --fasta=/mnt/isilon/xing_lab/liuc9/refdata/mgatk_index/genome.fa \
  --genes=/mnt/isilon/xing_lab/liuc9/refdata/mgatk_index/Homo_sapiens.GRCh38.107.new.gtf \
  --nthreads=50 \
  --memgb=200
