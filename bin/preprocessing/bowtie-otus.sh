# extract subset of query seqs to align with reference
tail -n 2000 seqs.fna > seqs-2000.fna

#1. align subset with the gg 97% OTUs
bsub -o align.lsf -q hour -R "rusage[mem=8]" "align_seqs.py -i seqs-2000.fna -t /home/unix/dknights/lib/bowtie2-2.0.2/bowtie2-build /seq/microbiome/RESOURCES/greengenes/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011_aligned.fasta -o align"

#2. trim ref seqs to match aligned query seqs, remove gaps
time python ~/bin/trim_ref_alignment.py $1 align/seqs-2000_aligned.fasta

#4. build bowtie db from trimmed ref seqs
mkdir bowtie
cd bowtie
bsub -o btbuild -q priority -R "rusage[mem=8]" "/home/unix/dknights/lib/bowtie2-2.0.2/bowtie2-build /seq/microbiome/RESOURCES/greengenes/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011_aligned.fasta refset"
cd ..

#5. run alignment
bsub -o btalign -q priority -R "rusage[mem=8]" "/home/unix/dknights/lib/bowtie2-2.0.2/bowtie2 --very-fast -x bowtie/refset -U seqs.fna --no-head --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-.03 -f > btout.txt"

#6. convert to QIIME otu table
time python ~/bin/sam2uclust.py btout.txt > btclusters.txt
make_otu_table.py -i btclusters.txt -t /seq/microbiome/RESOURCES/greengenes/gg_otus_4feb2011/taxonomies/greengenes_tax.txt -o otus.txt
convert_biom.py -i otus.biom -b -o otus.txt --header_key taxonomy --output_metadata_id "Consensus Lineage" --process_obs_metadata taxonomy
sed -i -e 's:; :;:g' otus.txt
