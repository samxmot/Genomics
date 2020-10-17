from snps import SNPs
import os
import errno
import gzip



'''gene_path = input("Input path to user genomic file: ") 
trait_path = input("Input path to PGS Catalog trait file: ") 
if not os.path.isfile(gene_path):
	raise FileNotFoundError(
    	errno.ENOENT, os.strerror(errno.ENOENT), gene_path)
if not os.path.isfile(trait_path):
	raise FileNotFoundError(
    	errno.ENOENT, os.strerror(errno.ENOENT), trait_path)'''


def build_snp_weight_dictionary(file_path):
	snp_weigh_dictionary = {}
	with gzip.open(trait_path) as f:
		for line in f:
			line = line.decode("utf-8")
			field = line.split("\t")
			if field[0] == 'rsID':
				fields = field
				break

		rsID_idx = (fields.index('rsID'))
		effect_allele_idx = (fields.index('effect_allele'))
		effect_weight_idx = (fields.index('effect_weight'))
	with gzip.open(trait_path) as f:
		for line in f:
			line = line.decode("utf-8")
			field = line.split("\t")
			if field[0][:2] == 'rs':
				snp_weigh_dictionary[field[0]] = (field[effect_allele_idx], field[effect_weight_idx])
	return snp_weigh_dictionary


trait_path = "trait/bmi/PGS000298.txt.GZ"
gene_path = "data/9876.23andme.8136"
s = SNPs(gene_path)
snp_weigh_dictionary = {}
df = s.snps
weight_dictionary = build_snp_weight_dictionary(trait_path)
print(weight_dictionary)