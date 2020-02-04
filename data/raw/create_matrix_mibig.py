import sys

f = open('../pfams/pfam_list.txt')
domains = []
for line in f.readlines():
	domains.append(line.split()[0])
f.close()

f = open('mibig_hits_filtered.crh')
from collections import defaultdict
bgcs = defaultdict(lambda: defaultdict(int))
ordered = defaultdict(list)
for line in f.readlines():
	if not line.startswith("#") and line.strip() != '':
		prot = line.split()[0]
		bgc = "_".join(prot.split("_")[:-1])
		domain = line.split()[1]
		bgcs[bgc][domain] += 1
		ordered[bgc].append(domain)
f.close()


print('BGC', end='\t')
for d in domains:
	print(d, end="\t")
print("Order")
print("")
for bgc in bgcs:
	print(bgc, end="\t")
	for d in domains:
		print(str(bgcs[bgc][d]), end="\t")

	print(",".join(ordered[bgc]))
	print("")
