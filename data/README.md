Data files:

```
./raw/                              - MiBIG files in various stages of processing
mibig_pfam20.txt                    - cath-resolve-hits filtered PFAM20 transporter hits for mibig2
mibig_all_pfams.txt                 - cath-resolve-hits filtered PFAM20, SBP, and Biosynthetic hits for mibig2
```

```
./cathdb/
antismash_cathdb.tsv                - Counts of transporter CATHDB HMMs in antiSMASH BGCs
gcfs_cathdb.tsv                     - Counts of transporter CATHDB HMMs in Gene Cluster Family BGCs
mibig_cathdb.tsv                    - Counts of transporter CATHDB HMMs in Mibig2.
```

```
./functional_groups/fg_table.tsv    - Chemical functional groups for all mibig2 BGCs with chemical structures. 
The other files in this directory create the above table.
```

```
./labels/
extra_antibiotics                   - Hand curated antibiotic labels for more mibig2 BGCs
extra_antifungal                    - Hand curated antifungal labels for more mibig2 BGCs
extra_siderophore                   - Hand curated siderophore labels for more mibig2 BGCs
specific_antibiotics                - A list of specific well studied antibiotic examples
```

```
./metadata/
antismash_biosynthetic_matrix.tsv   - Counts of common Biosynthetic PFAMs in each antiSMASH BGC.
biosynthetic_mibig2.tsv             - Coutns of common Biosynthetic PFAMs in each mibig2 BGC
antismash_metadata.tsv              - Metadata for antiSMASH BGCs.
final_fun_families.tsv              - A filtered version of original_fun_families.tsv that hand removes some HMMs that are not transporter specific.
mibig2_metadata.tsv                 - Metadata for mibig2 BGCs.
original_fun_families.tsv           - CathDB functional family HMMs and the TCDB families they correspond to / hit
tcdb_families.tsv                   - Names and Types (import vs export) for all Transporter Classificaiton Database (TCDB) labels
```

```
./pfam/
antismash_pfam20.tsv                - Counts of 20 important transporter PFAMs in antiSMASH BGCs
gcf_pfams.tsv                       - Counts of important transporter PFAMs in GCF BGCs
mibig2_sbp_hits.tsv                 - Counts of Substrate Binding Protein PFAMs in Mibig2
mibig_pfam20.tsv                    - Counts of 20 important transporter PFAMs in Mibig2
