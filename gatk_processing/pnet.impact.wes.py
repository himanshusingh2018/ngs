import pandas as pd, os, pysam, sys

#GENERATE PICARD CLEAN AND MARK DUPLICATES
def picard_process(bam_base_dir, bam):
    input = f'{bam_base_dir}{bam[0]}/{bam[1]}/{bam}.bam'
    output_picard_clean = f'{picard}{bam}.clean.bam'
    cmd_picard_clean = f'picard -Xmx64g CleanSam INPUT={input} OUTPUT={output_picard_clean}'
    print(cmd_picard_clean)
    os.system(cmd_picard_clean)
    print(f'\nSuccessfully generated:\tpicard/{bam}.clean.bam\n')

    #GENERATE PICARD MARKDUPLICATES
    cmd_picard_mkdup = f"picard -Xmx64g MarkDuplicates INPUT={output_picard_clean} OUTPUT={picard}{bam}.mkdup.bam METRICS_FILE={picard}{bam}.mkdup.txt VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"
    os.system(cmd_picard_mkdup)
    print(cmd_picard_mkdup)
    print(f'\nSuccessfully generated:\t:picard/{bam}.mkdup.bam\n')

def read_keys(keyfile):
    header = ['DMP ID', 'Mirror/Anonymized ID', 'Patient_Group_ID', 'PartC Constent Status', 'OncoTree TumorType code',
              'SampleType', 'PrimarySite', 'MetastasisSite', 'TissueType', 'Anonymized ProjectName',
              'SampleCoverage', 'SomaticCalling Status', 'Major Allele Contamination', 'Minor Allele Contamination']
    df = pd.read_csv(keyfile, sep=",",header=None, names=header, usecols=['DMP ID', 'Mirror/Anonymized ID'])
    return(df)

def read_sample_id(pnet_samples):
    df = pd.read_csv(pnet_samples, sep="\t", header=0, usecols=['DMP ID'])
    df = df['DMP ID'].str.split('-', n=2, expand=True)[0] + '-' + df['DMP ID'].str.split('-', n=2, expand=True)[1]
    return(df)

def sample_name(bam_base_dir, bam):#extract bam header
    sample_name = pysam.AlignmentFile(f'{bam_base_dir}{bam[0]}/{bam[1]}/{bam}.bam','rb').header['RG'][0].get('SM')
    return(sample_name)

jvm_memory = '60G'
hg19_ref = '/home/singhh5/hg19/Homo_sapiens_assembly19.fasta'
keyfile = '/juno/res/dmpcollab/dmprequest/12-245/key.txt'
header = '/juno/res/dmpcollab/dmprequest/12-245/header.txt'
bam_base_dir = '/juno/res/dmpcollab/dmpshare/share/irb12_245/'#BAM file directory
keys = read_keys(keyfile)#read all keys
picard = '/home/singhh5/juno_project/picard/'
vcf = '/home/singhh5/juno_project/vcf/'
#samp_id = 'P-0000850'
samp_id = sys.argv[1]
sample = keys.loc[keys['DMP ID'].str.contains(samp_id)]#extract sample from key
id_list = sample['Mirror/Anonymized ID'].tolist() if not sample.empty else []

normal_sample_header = []
for element in id_list:
     if os.path.exists(f'{bam_base_dir}{element[0]}/{element[1]}/{element}.bam'):
        picard_process(bam_base_dir,element)
        if element.endswith('-N'):
            print(element)
            normal_sample_header.append(sample_name(bam_base_dir, element))
print(normal_sample_header)
cmd_gatk_tumor_vs_normal = f'''
gatk Mutect2 \
--java-options -Xmx{jvm_memory} \
-R {hg19_ref} \
-I {" -I ".join([f"{picard}{s}.mkdup.bam" for s in id_list])} \
-normal {" -normal ".join([s for s in normal_sample_header])} \
--native-pair-hmm-threads 10 \
-O {vcf}{samp_id}.vcf.gz
'''
os.system(cmd_gatk_tumor_vs_normal)
print(cmd_gatk_tumor_vs_normal)
print(f'Successfully generated:\t{samp_id}.vcf.gz')

os.system('rm '+' '.join([f"{picard}{s}*" for s in id_list]))
