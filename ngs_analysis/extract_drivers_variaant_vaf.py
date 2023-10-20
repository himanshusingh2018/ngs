import os, pandas as pd, sys, gzip
def read_vcf_gz(PID, fvcf):
    header = [line.strip().split('\t') for line in gzip.open(f'{PID}/vcf/{fvcf}', 'rt') if line.startswith('#CHROM')]
    vcf = pd.read_csv(f'{PID}/vcf/{fvcf}', header=None, compression='gzip', sep="\t", comment='#', names=header[0])
    vcf.rename(columns={'#CHROM': 'Chromosome', 'POS': 'Position'}, inplace=True)
    return(vcf)

driver_gene = pd.DataFrame([{"Gene": "ARID1A", "chr": "chr1", "start": 27022522, "end": 27108601},
                            {"Gene": "SETD2", "chr": "chr3", "start": 47057898, "end": 47205467},
                            {"Gene": "DAXX", "chr": "chr6", "start": 33286335, "end": 33290793},
                            {"Gene": "BRAF", "chr": "chr7", "start": 140433813, "end": 140624564},
                            {"Gene": "PTEN", "chr": "chr10", "start": 89623195, "end": 89728532},
                            {"Gene": "MEN1", "chr": "chr11", "start": 64570986, "end": 64578188},
                            {"Gene": "KMT2D", "chr": "chr12", "start": 49412758, "end": 49449107},
                            {"Gene": "TSC2", "chr": "chr16", "start": 2097990, "end": 2138713},
                            {"Gene": "TP53", "chr": "chr17", "start": 7571720, "end": 7590868},
                            {"Gene": "ATRX", "chr": "chrX", "start": 76760356, "end": 77041719}])

def extract_driver_variant(list_pos):
    a = vcf.loc[vcf['Position'].isin(list_pos)]
    a = a[['Chromosome', 'Position', 'REF', 'ALT'] + [col for col in a.columns if col.startswith('T')]]

    for col in [c for c in a.columns if c.startswith('T')]:
        a[col] = a[col].str.split(':').str[1] + '|' + a[col].str.split(':').str[2]
    print(a.reset_index(drop=True))

os.chdir('/Volumes/lilac_data_ziv/transciptome/paired_pnet/wes_analysis/')#set input dir path
vcf = read_vcf_gz(PID=sys.argv[1], fvcf=sys.argv[2])#Read vcf file
vcf['Chromosome'] = 'chr' + vcf['Chromosome'].astype(str)

#extracting driver mutations
driver_mutation = pd.DataFrame()
for i, row in driver_gene.iterrows():
    print(row['Gene'])
    df = vcf[(vcf.FILTER == 'PASS') & (row['chr'] == vcf['Chromosome']) & (vcf['Position'] >= row['start']) & (vcf['Position'] <= row['end'])]
    df.insert(0, 'Gene', row['Gene'])# Insert 'Gene' column at the beginning
    # Apply the split operation to columns from index 11 to the end
    cols_to_process = df.iloc[:, 10:]
    result = cols_to_process.apply(lambda x: x.str.split(':').str[1:3])
    df.iloc[:, 10:] = result
    driver_mutation = pd.concat([driver_mutation, df], ignore_index=True)# Append the DataFrame to driver_mutation

driver_mutation.reset_index(drop=True, inplace=True)# Reset the index
driver_mutation.to_csv(f'{sys.argv[1]}/{sys.argv[1]}.driver_mutation.txt',sep="\t",header=True,index=False)
#extracting variant allele depth
variant_ad = pd.DataFrame()
for i, row in driver_gene.iterrows():
    print(row['Gene'])
    df = vcf[(vcf.FILTER == 'PASS') & (row['chr'] == vcf['Chromosome']) & (vcf['Position'] >= row['start']) & (vcf['Position'] <= row['end'])]
    df.insert(0, 'Gene', row['Gene'])# Insert 'Gene' column at the beginning
    # Apply the split operation to columns from index 11 to the end
    cols_to_process = df.iloc[:, 10:]
    result = cols_to_process.apply(lambda x: x.str.split(':').str[1]).apply(lambda x: x.str.split(',').str[1])
    df.iloc[:, 11:] = result
    variant_ad = pd.concat([variant_ad, df], ignore_index=True)# Append the DataFrame to driver_mutation
variant_ad = variant_ad[['Gene'] + variant_ad.columns[11:].tolist()].drop_duplicates(keep='first')

variant_ad.to_csv(f'/Users/singhh5/Desktop/mutation/result/{sys.argv[1]}.variant_ad.txt',sep=",",header=True,index=False)
#variant_ad[['Gene'] + variant_ad.columns[11:].tolist()].to_csv(f'/Users/singhh5/Desktop/mutation/result/{sys.argv[1]}.variant_ad.txt',
#                                                               sep=",",header=True,index=False)

'''
set input dir path
os.chdir('/Volumes/lilac_data_ziv/transciptome/paired_pnet/wes_analysis/')#set input dir path
Run the script
python test.extract_drivers_variaant_vaf.py PID_145 pid145.filter.vcf.gz
output dir: /Volumes/lilac_data_ziv/transciptome/paired_pnet/wes_analysis/
'''
