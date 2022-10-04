## HOMER analysis
module load homer/4.8

id='AC_uF_IRF2-44'
size=200
bedfile="${id}.bed"
genome='/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/genome.fasta'
outdir="homer_${id}"

findMotifsGenome.pl \
${bedfile} \
${genome} \
${outdir} \
-size ${size}
