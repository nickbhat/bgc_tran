import sys

f = open('final_fun_families.tsv')
domains = []
for line in f.readlines():
	domains.append(line.split()[0])
f.close()
f = open('mibig_cathdb.crh')
from collections import defaultdict
bgcs = defaultdict(lambda: defaultdict(int))
for line in f.readlines():
	if not line.startswith("#") and line.strip() != '':
		prot = line.split()[0]
		prot = "_".join(prot.split("_")[:-1])
		domain = line.split()[1]
		bgcs[prot][domain] += 1
f.close()


print('BGC', end='\t')
for d in domains:
	print(d, end="\t")
print("")
for bgc in bgcs:
	print(bgc, end="\t")
	for d in domains:
		print(str(bgcs[bgc][d]), end="\t")
	print("")
