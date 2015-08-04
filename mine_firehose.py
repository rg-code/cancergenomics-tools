import os
import csv
import re

cancer_type_list = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'COADREAD', 'DLBC', 'ESCA', 'FPPP', 'GBM', 'GBMLGG', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']

primary_solid_tumor_BCR = []
metastatic_BCR = []
normal_BCR = []
primary_solid_tumor_Clinical = []
metastatic_Clinical = []
normal_Clinical = []
primary_solid_tumor_CN = []
metastatic_CN = []
normal_CN = []
primary_solid_tumor_LowP = []
metastatic_LowP = []
normal_LowP = []
primary_solid_tumor_Methylation = []
metastatic_Methylation = []
normal_Methylation = []
primary_solid_tumor_mRNA = []
metastatic_mRNA = []
normal_mRNA = []
primary_solid_tumor_mRNASeq = []
metastatic_mRNASeq = []
normal_mRNASeq = []
primary_solid_tumor_miR = []
metastatic_miR = []
normal_miR = []
primary_solid_tumor_miRSeq = []
metastatic_miRSeq = []
normal_miRSeq = []
primary_solid_tumor_RPPA = []
metastatic_RPPA = []
normal_RPPA = []
primary_solid_tumor_MAF = []
metastatic_MAF = []
normal_MAF = []

for cancer in cancer_type_list:
	file_name_reference = open(cancer + '.html')
	reference = file_name_reference.readlines()
	for o in range(0, len(reference)):
		section_name = re.search('div\sclass="sectionheader"',reference[o])
		if (section_name is None):
			continue
		else:
			for p in range(o, len(reference)):
				section_end = re.search('/table',reference[p])
				if section_end is None:
					continue
				else:
					end_number = p
			linesjoin = reference[o:o+end_number]
			linesjoined = ''.join(linesjoin)
			patients = re.findall('(?<=TCGA-)(.*?)(?=\n)',linesjoined, re.DOTALL)
			primary_solid_tumor_search = re.search('Primary\sSolid\sTumor',reference[o+1])
			metastatic_search = re.search('Metastatic',reference[o+1])
			blood_derived_normal_search = re.search('Blood\sDerived\sNormal',reference[o+1])
			solid_tissue_normal_search = re.search('Solid\sTissue\sNormal',reference[o+1])
			BCR_search = re.search('BCR',reference[o+1])
			Clinical_search = re.search('Clinical',reference[o+1])
			CN_search = re.search('CN',reference[o+1])
			LowP_search = re.search('LowP',reference[o+1])
			Methylation_search = re.search('Methylation',reference[o+1])
			mRNA_search = re.search('mRNA',reference[o+1])
			mRNASeq_search = re.search('mRNASeq',reference[o+1])
			miR_search = re.search('miR',reference[o+1])
			miRSeq_search = re.search('miRSeq',reference[o+1])
			RPPA_search = re.search('RPPA',reference[o+1])
			MAF_search = re.search('MAF',reference[o+1])
			if primary_solid_tumor_search is not None:
				if BCR_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_BCR:
							primary_solid_tumor_BCR.append(patients[r])
				if Clinical_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_Clinical:
							primary_solid_tumor_Clinical.append(patients[r])
				if CN_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_CN:
							primary_solid_tumor_CN.append(patients[r])
				if LowP_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_LowP:
							primary_solid_tumor_LowP.append(patients[r])
				if Methylation_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_Methylation:
							primary_solid_tumor_Methylation.append(patients[r])
				if mRNA_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_mRNA:
							primary_solid_tumor_mRNA.append(patients[r])
				if mRNASeq_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_mRNASeq:
							primary_solid_tumor_mRNASeq.append(patients[r])
				if miR_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_miR:
							primary_solid_tumor_miR.append(patients[r])
				if miRSeq_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_miRSeq:
							primary_solid_tumor_miRSeq.append(patients[r])
				if RPPA_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_RPPA:
							primary_solid_tumor_RPPA.append(patients[r])
				if MAF_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in primary_solid_tumor_MAF:
							primary_solid_tumor_MAF.append(patients[r])
			if metastatic_search is not None:
				if BCR_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_BCR:
							metastatic_BCR.append(patients[r])
				if Clinical_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_Clinical:
							metastatic_Clinical.append(patients[r])
				if CN_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_CN:
							metastatic_CN.append(patients[r])
				if LowP_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_LowP:
							metastatic_LowP.append(patients[r])
				if Methylation_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_Methylation:
							metastatic_Methylation.append(patients[r])
				if mRNA_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_mRNA:
							metastatic_mRNA.append(patients[r])
				if mRNASeq_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_mRNASeq:
							metastatic_mRNASeq.append(patients[r])
				if miR_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_miR:
							metastatic_miR.append(patients[r])
				if miRSeq_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_miRSeq:
							metastatic_miRSeq.append(patients[r])
				if RPPA_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_RPPA:
							metastatic_RPPA.append(patients[r])
				if MAF_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in metastatic_MAF:
							metastatic_MAF.append(patients[r])
			if blood_derived_normal_search is not None:
				if BCR_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_BCR:
							normal_BCR.append(patients[r])
				if Clinical_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_Clinical:
							normal_Clinical.append(patients[r])
				if CN_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_CN:
							normal_CN.append(patients[r])
				if LowP_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_LowP:
							normal_LowP.append(patients[r])
				if Methylation_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_Methylation:
							normal_Methylation.append(patients[r])
				if mRNA_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_mRNA:
							normal_mRNA.append(patients[r])
				if mRNASeq_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_mRNASeq:
							normal_mRNASeq.append(patients[r])
				if miR_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_miR:
							normal_miR.append(patients[r])
				if miRSeq_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_miRSeq:
							normal_miRSeq.append(patients[r])
				if RPPA_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_RPPA:
							normal_RPPA.append(patients[r])
				if MAF_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_MAF:
							normal_MAF.append(patients[r])
			if solid_tissue_normal_search is not None:
				if BCR_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_BCR:
							normal_BCR.append(patients[r])
				if Clinical_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_Clinical:
							normal_Clinical.append(patients[r])
				if CN_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_CN:
							normal_CN.append(patients[r])
				if LowP_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_LowP:
							normal_LowP.append(patients[r])
				if Methylation_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_Methylation:
							normal_Methylation.append(patients[r])
				if mRNA_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_mRNA:
							normal_mRNA.append(patients[r])
				if mRNASeq_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_mRNASeq:
							normal_mRNASeq.append(patients[r])
				if miR_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_miR:
							normal_miR.append(patients[r])
				if miRSeq_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_miRSeq:
							normal_miRSeq.append(patients[r])
				if RPPA_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_RPPA:
							normal_RPPA.append(patients[r])
				if MAF_search is not None:
					for r in range(0, len(patients)):
						if patients[r] not in normal_MAF:
							normal_MAF.append(patients[r])

f = open('primary_solid_tumor_BCR' + '.txt', 'w')
for i in primary_solid_tumor_BCR:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_BCR' + '.txt', 'w')
for i in metastatic_BCR:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_BCR' + '.txt', 'w')
for i in normal_BCR:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_Clinical' + '.txt', 'w')
for i in primary_solid_tumor_Clinical:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_Clincial' + '.txt', 'w')
for i in metastatic_Clinical:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_Clinical' + '.txt', 'w')
for i in normal_Clinical:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_CN' + '.txt', 'w')
for i in primary_solid_tumor_CN:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_CN' + '.txt', 'w')
for i in metastatic_CN:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_CN' + '.txt', 'w')
for i in normal_CN:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_LowP' + '.txt', 'w')
for i in primary_solid_tumor_LowP:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_LowP' + '.txt', 'w')
for i in metastatic_LowP:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_LowP' + '.txt', 'w')
for i in normal_LowP:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_Methylation' + '.txt', 'w')
for i in primary_solid_tumor_Methylation:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_Methylation' + '.txt', 'w')
for i in metastatic_Methylation:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_Methylation' + '.txt', 'w')
for i in normal_Methylation:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_mRNA' + '.txt', 'w')
for i in primary_solid_tumor_mRNA:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_mRNA' + '.txt', 'w')
for i in metastatic_mRNA:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_mRNA' + '.txt', 'w')
for i in normal_mRNA:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_mRNASeq' + '.txt', 'w')
for i in primary_solid_tumor_mRNASeq:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_mRNASeq' + '.txt', 'w')
for i in metastatic_mRNASeq:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_mRNASeq' + '.txt', 'w')
for i in normal_mRNASeq:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_miR' + '.txt', 'w')
for i in primary_solid_tumor_miR:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_miR' + '.txt', 'w')
for i in metastatic_miR:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_miR' + '.txt', 'w')
for i in normal_miR:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_miRSeq' + '.txt', 'w')
for i in primary_solid_tumor_miRSeq:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_miRSeq' + '.txt', 'w')
for i in metastatic_miRSeq:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_miRSeq' + '.txt', 'w')
for i in normal_miRSeq:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_RPPA' + '.txt', 'w')
for i in primary_solid_tumor_RPPA:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_RPPA' + '.txt', 'w')
for i in metastatic_RPPA:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_RPPA' + '.txt', 'w')
for i in normal_RPPA:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('primary_solid_tumor_MAF' + '.txt', 'w')
for i in primary_solid_tumor_MAF:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('metastatic_MAF' + '.txt', 'w')
for i in metastatic_MAF:
	f.write("TCGA-%s\n" % i)
f.close()
f = open('normal_MAF' + '.txt', 'w')
for i in normal_MAF:
	f.write("TCGA-%s\n" % i)
f.close()
