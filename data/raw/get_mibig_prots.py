from Bio import SeqIO
import glob

for fn in glob.glob('./mibig_gbk_2.0/*.gbk'):
        name = fn.split("/")[-1].split(".")[0]
        i = 0
        try:
                for record in SeqIO.parse(fn, "genbank"):
                        for feature in record.features:
                                if feature.type == 'CDS':
                                        i += 1
                                        print(">" + name + "_" + str(i))
                                        print(feature.qualifiers['translation'][0])
        except:
                pass

