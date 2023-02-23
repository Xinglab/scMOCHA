

# cellranger mkref --config=/home/liuc9/data/refdata/mgatk_index/reference.json


cellranger mkref --genome=Human --fasta=/mnt/isilon/xing_lab/liuc9/refdata/mgatk_index/genome.fa --genes=/mnt/isilon/xing_lab/liuc9/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf

# TODO UPDATE GTF file by removing contigs that not in fasta then makeref again,