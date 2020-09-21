import pandas as pd

# generate TAD results
def splitTAD(inputFile, outputFile):
    df = pd.read_csv(inputFile, sep = '\t', encoding='utf-8')
    tad = df[~df['boundary_strength_200000'].isnull()]
    tad = tad.loc[:, ['chrom', 'start', 'end', 'boundary_strength_200000']]
    res = pd.DataFrame()
    for i in range(1, 23):
        chromTad = tad[tad['chrom']=='chr'+str(i)]
        for j in range(1, len(chromTad)):
            chromTad.iloc[j, 1] = chromTad.iloc[j-1, 2]
        res = pd.concat([res, chromTad[1:]])
    res.to_csv(outputFile, sep='\t', index = False)
    del(res)


# generate compartment results
def splitCompartment(inputFile, outputFile, chrom_num):
    df = pd.read_csv(inputFile, sep = '\t', encoding='utf-8')
    com = df.loc[:, ['chrom', 'start', 'end', 'E1']]
    res = pd.DataFrame()
    for i in range(1, chrom_num + 1):
        chromCom = com[com['chrom']=='chr'+ str(i)]
        tmp = 0
        chromCom = chromCom.reset_index(drop=True)
        comIndex = []
        for index, row in chromCom.iterrows():
            if chromCom.iloc[index, 3]*tmp < 0:
                comIndex.append(index)
            tmp = chromCom.iloc[index, 3]
        chromCom = chromCom.iloc[comIndex,:]
        for k in range(1, len(chromCom)):
            chromCom.iloc[k, 1] = chromCom.iloc[k-1, 2]
        res = pd.concat([res, chromCom[1:]])
    res.to_csv(outputFile, sep='\t', index = False)
