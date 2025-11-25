import argparse
import os
import gzip
import json
import logging
import sys
sys.path.append("/mnt/rawdata4/tests/linsh2/scripts/work/")

from collections import OrderedDict
from variant_wkf.settings import client1 as client

class parameterrs():
	def __init__(self):
		doc = '血液病背噪文件更新'
		self.parser = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)

		self.parser.add_argument('-i', '--input', dest='input', metavar='input file', type=str, help="the input file")

		self.parser.add_argument('-m', '--mode', dest='mode', metavar='leukemia/lymphoma', type=str, help="the input file")

		self.parser.add_argument('-o', '--output', dest='output', metavar='path/to/output', type=str, help="output file")




		###参数解析
		self.args = self.parser.parse_args()

	def log_file_func(self):
		logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s', filemode='a',
							 stream=sys.stdout)

	def repeat(self):
		logging.info('output      : ' + self.args.output)
		logging.info('mode      : ' + self.args.mode)
		if self.args.input:
			logging.info('input       : ' + self.args.input)

def mongo_find_NGS_Report():
	db = client.NGS_Report
	db.authenticate(name='admin', password='ngs2020!')
	collection = db.bioinfo_result
	info = collection.find({'bio_version': '"blood_pipeline_v1.2"'})
	return info

project_panel = {'血液肿瘤180基因配对检测': 'leukaemia180_spikein', '血液肿瘤180基因单肿瘤检测': 'leukaemia180_spikein', '髓系肿瘤128基因配对检测': 'leukaemia180_spikein', '髓系肿瘤128基因单肿瘤检测': 'leukaemia180_spikein', 'AML_96基因配对检测': 'leukaemia180_spikein', 'AML_96基因单肿瘤检测': 'leukaemia180_spikein', 'MDS_68基因配对检测': 'leukaemia180_spikein', 'MDS_68基因单肿瘤检测': 'leukaemia180_spikein', 'ALL_96基因配对检测': 'leukaemia180_spikein', 'ALL_96基因单肿瘤检测': 'leukaemia180_spikein', 'MPN_56基因配对检测': 'leukaemia180_spikein', 'MPN_56基因单肿瘤检测': 'leukaemia180_spikein', 'MDS_MPN_48基因配对检测': 'leukaemia180_spikein', 'MDS_MPN_48基因单肿瘤检测': 'leukaemia180_spikein', '淋巴瘤_218基因配对检测': 'Lymphoma218_spikein', '淋巴瘤_218基因单肿瘤检测': 'Lymphoma218_spikein', 'B细胞淋巴瘤164基因配对检测': 'Lymphoma218_spikein', 'B细胞淋巴瘤164基因单肿瘤检测': 'Lymphoma218_spikein', 'DLBCL_140基因配对检测': 'Lymphoma218_spikein', 'DLBCL_140基因单肿瘤检测': 'Lymphoma218_spikein', 'TNK细胞淋巴瘤95基因配对检测': 'Lymphoma218_spikein', 'TNK细胞淋巴瘤95基因单肿瘤检测': 'Lymphoma218_spikein', 'CLL_100基因配对检测':'Lymphoma218_spikein', 'CLL_100基因单肿瘤检测':'Lymphoma218_spikein'}

if __name__ == '__main__':
	params = parameterrs()
	params.log_file_func()
	params.repeat()

	###读入老版背噪数据
	inputfile = params.args.input
	mode = params.args.mode
	outfile = params.args.output

	fa = open(inputfile,'r')
	data = json.load(fa)
	fa.close()
	for item in data:
		print (f'{item}:{data[item]}')
		break

	samples_mongo = mongo_find_NGS_Report()
	variant_report = {}
	for item in samples_mongo:
		if 'test_info' not in item:
			continue
		sequencing_time = item['sequencing_time']
		sample_ID = item['sample_ID']
		project_type = item['project_type'].replace('基因检测(含胚系配对)','基因配对检测')

		###根据panel过滤
		if mode == 'leukemia':
			if project_panel[project_type] != 'leukaemia180_spikein':
				continue

		###打开本地文件
		file = f"/mnt/analysis6/NextSeq-NS500/{sequencing_time}_ngs/XueYeBing/backup/{sample_ID}/{sample_ID}.all.txt"
		try:
			fa = open(file, 'r')
		except:
			print (f'{file} is not existed.')
			continue
		###找到hgvs和物理位置对应关系
		falines = fa.readlines()
		fa.close()
		headline = falines[0].strip('\n\r').split('\t')
		hgvs_2_pos = {}
		for i in range(1,len(falines)):
			if falines[i].startswith('filter'):
				continue
			falines[i] = falines[i].strip('\n\r').split('\t')
			line = {}
			for o in range(len(headline)):
				line[headline[o]] = falines[i][o]

			key = f"{line['Gene']}:{line['HGVS_c']}"
			value = f"{line['chrom']}:{line['POS']}:{line['REF']}:{line['ALT']}"
			hgvs_2_pos[key] = value

		short_variant = item['short_variant']
		###判断该样本是否存在变异合并
		merge_value = 0
		for variant in short_variant:
			key = f"{variant['Gene_name']}:{variant['HGVS_c']}"
			if key not in hgvs_2_pos:
				merge_value = 1
				break
		if merge_value == 1:
			continue

		for variant  in short_variant:
			key = f"{variant['Gene_name']}:{variant['HGVS_c']}"
			pos = hgvs_2_pos[key]
			hgvs = f"{variant['Gene_name']}:{variant['Refseq_mRNA_Id']}:{variant['HGVS_c']}"
			Tumor_mutant_frequency = float(variant['Tumor_mutant_frequency'].strip('%'))
			status = variant['status']

			if pos not in variant_report:
				variant_report[pos] = {"samples": 1, "hgvs": hgvs, "vaf": {status: [Tumor_mutant_frequency]}}
			else:
				variant_report[pos]['samples'] = variant_report[pos]['samples'] + 1
				if status not in variant_report[pos]['vaf']:
					variant_report[pos]['vaf'][status] = [Tumor_mutant_frequency]
				else:
					variant_report[pos]['vaf'][status].append(Tumor_mutant_frequency)

	for item in variant_report:
		if variant_report[item]['samples'] < 3:
			continue
		if 'unknown' not in variant_report[item]['vaf']:
			continue
		cutoff = max(variant_report[item]['vaf']['unknown'])
		if item not in data:
			data[item] = {"hgvs": variant_report[item]['hgvs'], "cutoff": cutoff}
		else:
			data[item]["cutoff"] = cutoff
		print (f"{item}: {data[item]}")

	out = open(outfile, 'w')
	out.write(json.dumps(data))
	out.close()


