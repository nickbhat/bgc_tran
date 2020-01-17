def mibig_viewer(mibig_dir, mibig_transporters, bgc, to_label=['ABC_tran', 'BPD_transp_1','TonB_dep_Rec','ABC_membrane','ACR_tran','FecCD', 'ABC2_membrane','MatE', 'OEP', 'FtsX', 'MFS_3', 'MFS_1', 'ABC2_membrane_3','MacB_PCD', 'ABC2_membrane_4','MMPL', 'BPD_transp_2', 'Peripla_BP_2','SBP_bac_1','Peripla_BP_4','SBP_bac_5','SBP_bac_8']):
    from Bio import SeqIO
    from dna_features_viewer import GraphicFeature, GraphicRecord

    transporters = {}
    f = open(mibig_transporters)
    for line in f.readlines():
        if not line.startswith("#"):
            if line.split()[0] not in transporters:
                transporters[line.split()[0]] = line.split()[1]
            else:
                transporters[line.split()[0]] =  transporters[line.split()[0]] + " " + line.split()[1]
    f.close()
    
    
    features = []
    colors = {'biosynthetic': "#850000",
     'biosynthetic-additional': "#ea8686",
     'other': "#dbdbdb",
     'regulatory': "#7cd369",
     'resistance': "#307321",
     'transport': "#3c85cd"}

    genes = []
    length = 0
    last_end = 0
    i = 0
    for record in SeqIO.parse(mibig_dir.rstrip("/") + "/" + bgc + ".gbk", "genbank"):
        for feature in record.features:
            
            
            if feature.type == 'CDS':
                i += 1
                feature_name = bgc + "_" + str(i)
                
                try:
                    color = colors[feature.qualifiers['gene_kind'][0]]
                except:
                    color = colors['other']
                    
                if feature.location.start < last_end:
                    start = last_end + 1  
                else:
                    start = feature.location.start
                    
                if feature_name in transporters:
                    found = False
                    for transporter in to_label:
                        if transporter in transporters[feature_name]:
                            genes.append(GraphicFeature(start=start, end=feature.location.end, strand=feature.location.strand, color=color, label=transporters[feature_name]))
                            found = True
                            break
                    if not found:
                        genes.append(GraphicFeature(start=start, end=feature.location.end, strand=feature.location.strand, color=color))
                else:
                    genes.append(GraphicFeature(start=start, end=feature.location.end, strand=feature.location.strand, color=color))

                        
                if feature.location.end > length:
                    length = feature.location.end
                last_end = feature.location.end
    record = GraphicRecord(sequence_length=length, features=genes)
    record.plot(figure_width=15)
   
    #                 break
