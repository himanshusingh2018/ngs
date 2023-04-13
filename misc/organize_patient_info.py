import numpy as np
import pandas as pd

def organize_data(file, sheet, pid):
    '''
    RUN THE FUNCTION
    file = '/Volumes/lilac_data_ziv/transciptome/paired_pnet/Project_12502_F/ziv_batch1_sample_data_key.xlsx'
    sheet = 'Sample_Medication_IR'
    med, ir = organize_data(file, sheet,pid='35103404')

    INPUT FILE<Excel Sheet>
    Pt id	35103404	35103404	35250343	35250343
    67889	Start Date	Proce Date	Start Date	Proce Date
    35103	Medication	IP/OP		Medication	IP/OP
    35250	3/9/18	    OPTime/CPT	4/2/23	    OPTime/CP
    35257	ACETAMINO	Proce Des	HEPARIN 	Procedure Des
    35263	3/9/18	    Physician	            Physician
    35352	NORMAL SALI	1/3/18			        1/27/23
    35415	OP			OP
    35436	62305		S2095
    35447	SPINAL CANA	THREE DOSE
    38018	TISNADO,J	ZIV,ETAY
    38038	12/7/17		1/27/23
		    OP			OP
		    47000		79445
		    NBB NEEDL	THREE DOSE
		    GONZALEZ-A	ZIV,ETAY
		    5/24/16			
		    OP			
		    36245			
		    SIRSPHERE TR			
		    BOAS,FR						
    '''
    
    #Read the excel file
    df = pd.read_excel(file, sheet_name=sheet)
    df.columns = df.columns.astype(str)
    df=df.loc[:, ~df.columns.str.startswith('Unnamed')]

    #Extract and reshape medication data
    med = df[pid].dropna().values
    med = med.reshape(-1,2)
    med_final = pd.DataFrame(med, columns=['Medicine Start', 'Medicine'])
    
    #Extract and restructure IR data
    ir = df[f'{pid}.1'].dropna()
    header = ir.iloc[:5].tolist()
    ir_final = pd.DataFrame(columns=header)

    for i in range(5, len(ir),5):
        data = ir.iloc[i:i+5].tolist()
        row_df = pd.DataFrame([data], columns=header)
        ir_final = pd.concat([ir_final, row_df],ignore_index=True)

    return med_final, ir_final


file = '/Volumes/lilac_data_ziv/transciptome/paired_pnet/Project_12502_F/ziv_batch1_sample_data_key.xlsx'
sheet = 'Sample_Medication_IR'
#med, ir = organize_data(file, sheet,pid='35103404')
#print("Medication DataFrame:\n", med)
#print("IR DataFrame:\n", ir)
#print(med.shape)
#print(ir.shape)

# Assuming you have four PID values
pid_list = ['35103404','35250343','35257401',
            '35267185','35352845','35415284',
            '35436484','35447553','38017807',
            '38038217']

# Create an empty Excel file with a writer
with pd.ExcelWriter('out.xlsx') as writer:
    # Loop through the PID values
    for pid in pid_list:
        # Call the function to organize data and get the resulting DataFrames
        med, ir = organize_data(file, sheet, pid)
        # Write the DataFrames to respective sheets with sheet names as per the PID values
        med.to_excel(writer, sheet_name='med'+pid, index=False)
        ir.to_excel(writer, sheet_name='ir'+pid, index=False)

print("Data written to Excel file successfully!")
