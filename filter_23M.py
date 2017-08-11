import pandas as pd
import numpy as np
import pysam 
import sys
import ponytools as pc

from subprocess import Popen, PIPE
from ponytools import Fasta, MNEc2MAnnot


header_text = '''\
##fileformat=VCFv4.1
##fileDate=20170228
##dbSNP_meta_start
## Submitter Contact Information
##TYPE: CONT
##HANDLE: EGGL
##NAME: Robert Schaefer, Molly McCue
##EMAIL: schae234@umn.edu,mccu0173@umn.edu
##LAB: Equine Genetics and Genomics Lab
##INST: University of Minnesota
##||
##TYPE: PUB
##HANDLE: EGGL
##TITLE: Development of a high-density, 2M SNP genotyping array and 670k SNP imputation array for the domestic horse
##AUTHORS: Robert J. Schaefer, Mikkel Schubert, Ernest Bailey, Danika L. Bannasch, Eric Barrey, Gila Kahila Bar-Gal, Gottfried Brem, Samantha A. Brooks, Ottmar Distl, Ruedi Fries, Carrie J. Finno, Vinzenz Gerber, Bianca Haase, Vidhya Jagannathan, Ted Kalbfleisch, Tosso Leeb, Gabriella Lindgren, Maria Susana Lopes, Nuria Mach, Artur Câmara Machado, James N. MacLeod, Annette McCoy, Julia Metzger, Cecilia Penedo, Sagi Polani, Stefan Rieder, Imke Tammen, Jens Tetens, Georg Thaller, Andrea Verini-Supplizi, Claire M. Wade, Barbara Wallner, Ludovic Orlando,, James R. Mickelson, Molly E. McCue
##YEAR: 2017
##STATUS: 2
##||
##TYPE: METHOD
##HANDLE: EGGL
##ID: MNEc_SNP_Chip
##METHOD_CLASS: WGS
##SEQ_BOTH_STRANDS: YES
##TEMPLATE_TYPE: DIPLOID
##MULT_PCR_AMPLIFICATION: UNKNOWN
##MULT_CLONES_TESTED: NO
##METHOD: Whole genome, 100 base pair (bp), paired-end Illumina HiSeq reads were generated for 153 horses (including Twilight), representing 24 distinct breeds, at a depth between 1.7X and 64X, with a median depth of 13X. Read mapping was performed using the PALEOMIX genome mapping protocol to efficiently process samples in parallel and to assess individual sample quality. Each sample was mapped to the EquCab2.0 reference genome to produce a total of nearly 48 billion unique reads being aligned to the nuclear genome. Variants were identified by extending the PALEOMIX framework to identify SNPs using GATK and samtools. To maximize efficiency of variant calling in individual breeds, and to minimize bias due to variable sequencing depth of coverage, individuals were broken up into 16 variant calling groups by estimated depth of coverage and breed. Variants were called using permissive parameters in both the GATK UnifiedGenotyper as well as SAMtools ‘mpileup’ utilities. Approximately 23 million potential SNPs were called by both GATK and SAMtools. These 23 million SNPs kept for further analysis and validation.
##||
##TYPE: POPULATION
##HANDLE: EGGL
##ID: DISCOVERY
##POPULATION: Discovery Population for the MNEc2M Array 
##||
##TYPE:SNPASSAY
##HANDLE:EGGL
##BATCH: MNEc_23M
##MOLTYPE: Genomic
##METHOD:WGS
##COMMENT: Genotypes were discovered from Whole Genome Sequence as described in the METHOD section. Methods are further descibed in the aforemention Publication.
##PRIVATE:
##Note: 
###dbSNP_end_meta
##filedate=20170224
##handle=eggl
##batch=MNEc23M
##reference=GCF_000002305.2
##CONTIGS
##INFO=<ID=VRT,Number=1,Type=Integer,Description="Variation type,1 - SNV: single nucleotide variation,2 - DIV: deletion/insertion variation,3 -HETEROZYGOUS: variable, but undefined at nucleotide level,4 - STR: short tandem repeat (microsatellite) variation, 5 - NAMED: insertion/deletion variation of named repetitive element,6 - NO VARIATON: sequence scanned for variation, but none observed,7 - MIXED: cluster contains submissions from 2 or more allelic classes (not used),8 - MNV: multiple nucleotide variation with alleles of common length greater than 1,9 - Exception">
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="The Recalibrated Variant Quality Score described in Schaefer et al. (2017)">
##INFO=<ID=FLANK-5,Number=1,Type=String,Description="The 5 prime flanking sequence of the variant">
##INFO=<ID=FLANK-3,Number=1,Type=String,Description="The 3 prime flanking sequence of the variant">
##INFO=<ID=MNEc2M-ID,Number=1,Type=String,Description="The MNEc2M ID from Schaefer et. al (2017)">
##INFO=<ID=MNEc670k,Number=.,Type=Flag,Description="A Flag indicating that the SNP is present on the MNEc670k SNP Array.">
##FORMAT=<ID=AC,Number=1,Type=Integer,Description="The ALT allele count in the MNEc2M Discovery Cohort">
##FORMAT=<ID=FRQ,Number=.,Type=Float,Description="Frequency of the alternate allele in the MNEc2M Discovery Cohort">
##population_id=DISCOVERY
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDISCOVERY'''

def main():
    # Open Reference Files
    vcf23M = pysam.VariantFile('/project/MNEc2M/VCFs/MNEc/23M_WGS.vcf.gz')
    print('Populating MNEc Annotation File',file=sys.stderr)
    annot = pc.MNEcAnnot.MNEc2MAnnot()
    print('Populating FASTA files',file=sys.stderr)
    refseq = Fasta.from_file("GCF_000002305.2_EquCab2.0_genomic.fna",nickname=(r'.*chromosome ([\dX]+).*',r'chr\1'))
    fasta = Fasta.from_file("/project/Data/Fasta/EquCab2_wChrun1_2/Equus_cab_nucl_wChrUn1_2.fasta")
    print('Reading in blast results',file=sys.stderr)
    blast = pd.read_table("23M_WGS.blast.txt",names=['id','chr','perc','start','end'])
    # We only want blast sites that are to chrunk sites
    blast = blast[np.core.defchararray.startswith(blast.chr.values.astype(str),'NW')]
    blast = blast.groupby('id').first()
    blast = dict(zip(blast.index,zip(blast.chr,blast.start+35)))
    
    a = annot._annot.reset_index()
    MNEc2M_ids = dict(zip(zip(a.chrom,a.pos),a.MNEcID))
    MNEc670k_ids = set(a.query('in670==True').MNEcID)
    del a

    # generate a list of the contigs
    contigs = []
    for id,chrom in refseq.chroms.items():
        if id in refseq.nicknames:
            contigs.append('##contig=<ID={},length={},ALT={}>'.format(id,len(chrom),refseq.nicknames[id]))
        else:
            contigs.append('##contig=<ID={},length={}>'.format(id,len(chrom)))
    with open('23M_dbSNP.vcf.header','w') as HEADER:
        print(header_text.replace('##CONTIGS',"\n".join(contigs)),file=HEADER)
    with open('23M_dbSNP.vcf.body','w') as BODY:
        for i,var in enumerate(vcf23M):
            # Get the coordinates from blast
            flank5 = fasta[var.chrom][var.pos-35:var.pos-1].upper()
            flank3 = fasta[var.chrom][var.pos+1:var.pos+35].upper()
            sequence =  flank5 + var.ref + flank3
            chrom = var.chrom
            pos = var.pos

            id = "{}.{}{}>{}".format(chrom.replace('chr',''),pos,var.ref,var.alts[0])
            ref = var.ref
            alts = ','.join(var.alts)
            qual = '{:.2f}'.format(var.qual).rstrip('0').rstrip('.')
            passing = 'PASS'

            # Fix the chun
            if chrom == "chrUn":
                # Get the blast information
                try:
                    chrom2,pos2 = blast[chrom+'.'+str(pos)] 
                except KeyError as e:
                    continue 
                chrom = chrom2
                pos = pos2

            # Get the flank information
            info = ';'.join([
                'VQSLOD={:.4}'.format(float(var.info['VQSLOD'])),
                'VRT=1',
                'FLANK-5={}'.format(flank5),
                'FLANK-3={}'.format(flank3),
            ])
            fmt = '{}:{:.4}'.format(
                var.info['AC'][0],
                var.info['AF'][0],
            )
            if (var.chrom,var.pos) in MNEc2M_ids:
                mnid = MNEc2M_ids[(var.chrom,var.pos)] 
                info += ";MNEc2M-ID={}".format(mnid)
                if mnid in MNEc670k_ids:
                    info += ";MNEc670k"
            print('\t'.join(map(str,[
                chrom,
                pos,
                id,
                ref,
                alts,
                qual,
                passing,
                info,
                'AC:FRQ',
                fmt
            ])),file=BODY)
            if i % 100000 == 0 and i > 0:
                print('Processed {} records'.format(i),file=sys.stderr)

if __name__ == '__main__':
    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(
        mode='Verbose', color_scheme='Linux', call_pdb=1
    )       
    main()



