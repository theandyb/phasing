# Prepping data for TOPMed phasing

# The commands below assume that each individual VCF file has been indexed

# Example: 50 samples
# Generate list of files for merging
start=396
end=400

for i in {$start..$end}; do
echo "pair_${i}_test.bcf" >> file_list.txt
done

# Merge VCFs using bcftools
bcftools merge -0 --file-list file_list.txt | bcftools annotate -Oz -x "FILTER,INFO" > topmed_input/pair${start}_${end}.vcf.gz
rm file_list.txt

for i in `seq 63 79`; do
(( a = $i*5 + 1 ))
(( b = $i*5 + 5 ))
for i in `seq $a $b`; do
echo "pair_${i}_test.bcf" >> file_list.txt
done
bcftools merge -0 -Oz --file-list file_list.txt > topmed_input/pair${a}_${b}.vcf.gz
rm file_list.txt
done

# sftp://snowwhite.sph.umich.edu/home/beckandy/research/phasing/output/filter_switch_errors/vcf/topmed_input/pair391_395.vcf.gz

# After we have phased the above vcf, we will need to merge the corresponding
# "true" VCFs as well:

## Generate site list
for i in `seq 301 400`; do
echo "pair_${i}_true.bcf" >> file_list.txt
done

## Merge VCFs
bcftools merge -0 -Ov --file-list file_list.txt > vcf_gz/topmed/pair301_400_true.vcf

# To run whatshap, we'll have to unzip our VCFs first

# And then we let'r rip
whatshap compare --names truth,topmed --tsv-pairwise whatshap/pair1_50 pair1_50_true.vcf topmed_1_50.phased.vcf

# Also do vcftools as well
vcftools --vcf topmed_201_205.phased.vcf \
  --diff pair201_205_true.vcf \
  --diff-switch-error \
  --out  vcftools/pair_201_205
