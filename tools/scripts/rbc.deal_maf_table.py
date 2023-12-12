import sys, os
import pandas as pd
import argparse
import re


def _argparse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-m', dest='maf_file', action="store",help="maf file")
    parser.add_argument('-o', dest='out_file', action="store",help="out file")
    return parser.parse_args()

def select_value(value):
    a = str(value).split(",")
    if len(a) == 2:
        b = a[0].strip('[')
        return b
    else:
        return value

def main():
    _COLUMNS = [
        'Hugo_Symbol', 'Variant_Classification',"cDNA_Change",
        'Entrez_Gene_Id','NCBI_Build','Chromosome',
        'Start_Position','End_Position','Strand',
        # 'Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2',
        # "Genome_Change","Annotation_Transcript","Transcript_Strand","Transcript_Exon",
        # "Codon_Change","Protein_Change","Refseq_mRNA_Id","Refseq_prot_Id",
        # "ref_context", "t_alt_count", "t_ref_count",
        # "HGNC_Previous_Symbols","HGNC_Previous_Name","HGNC_Synonyms","HGNC_Name_Synonyms",
        # "HGNC_Chromosome","HGNC_Ensembl_Gene_ID","HGNC_Pubmed_IDs","HGNC_RefSeq_IDs",
        # "HGNC_Gene_Family_Name","HGNC_UCSC_ID(supplied_by_UCSC)",
        "dbSNP_ID","AF",
    ]

    # https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format_TCGAv2/#current-version-changes
    Variant_Classification = [
        'Missense_Mutation',
        'RNA',
        'Silent', # ABO: 297G>A
        'Splice_Site',
        'Nonsense_Mutation',  # ABCG2: c.376C>T
        'Frame_Shift_Del',  # FUT1: c.551_552delAG, c.881_882delTT
        'Frame_Shift_Ins', 
        'In_Frame_Del', 
        'In_Frame_Ins', 
        'Translation_Start_Site', 
        'Nonstop_Mutation', 
        'Targeted_Region',
    ]
    ##########################################
    parser = _argparse()

    maf_file = parser.maf_file
    out_file = parser.out_file

    # Check if maf_file exists. and is not empty
    if not os.path.exists(maf_file) or os.path.getsize(maf_file) == 0:
        print(f"maf file {maf_file} does not exist or is empty!")
        sys.exit(1)

    df = pd.read_csv(maf_file, sep="\t")
    # select data by column and variat classification.
    df = df[_COLUMNS]

    # filter Variant Classification
    if len(Variant_Classification) > 0:
        df = df[df['Variant_Classification'].isin(Variant_Classification)]

    # former dbSNP subsitution:
    former_dbsnp = {
        'rs1022328332': 'rs688976' # blood
    }
    df['dbSNP_ID'] = df['dbSNP_ID'].replace(former_dbsnp)
    # output

    ##############################
    ## fix AF column [0.5,0.5]
    df['AF'] = df['AF'].apply(select_value)
    df.to_csv(out_file, sep="\t", index = False, header=True)


if __name__ == "__main__":
    main()
