import cooler
import json
import pandas as pd

mcoolFile = '/store/qzhong/test/hanhan/Rao2014-GM12878-MboI-allreps-filtered.mcool::resolutions/10000'
chromosom = 'chr2'
inputFile = '~/zxhTest/DemoData/result/Rao_compartment.bed'
target = pd.read_csv(inputFile, sep='\t')
#targetRange = [target[target['chrom']==chromosom].iloc[0,1]] + target[target['chrom']==chromosom].iloc[:,2].tolist()
TADRange = [46890000, 47150000, 47630000,48020000 ,48340000,48760000,48980000,49140000,49830000,50390000,50740000]
targetRange = [46890000, 48770000, 50740000]
chromRange = '%s:%s-%s'%(chromosom, str(targetRange[0]), str(targetRange[-1]))
quantileForEdge = 0.96
genre = 'fragment_compartment1'

# 读取mcool文件
c = cooler.Cooler(mcoolFile)
df = c.matrix(balance=False, as_pixels=True, join=True).fetch(chromRange)
# 将mcool文件写为json文件
# 生成边的json文件
dfEdge = pd.DataFrame()
dfEdge['source'] = df['start1'].apply(str).str.cat(df['end1'].apply(str), sep='-')  #将片段位置作为ID
dfEdge['target'] = df['start2'].apply(str).str.cat(df['end2'].apply(str), sep='-')
dfEdge['weight'] = (df['count'] - df['count'].min())/(df['count'].max() - df['count'].min())
dfEdge = dfEdge[dfEdge['weight'] > dfEdge['weight'].quantile([quantileForEdge])[quantileForEdge]]

# 生成点的json文件
dfNode = pd.DataFrame()
dfNode = pd.concat([dfEdge['source'],dfEdge['target']]).value_counts().rename_axis('name').reset_index(name='value')

# 判断点的类别
dfNode['category']=0
for i in range(0, len(dfNode)):
    node = int(dfNode.iloc[i, 0].split('-')[0])
    for j in range(0, len(targetRange)): 
        if node >= targetRange[j] and node<=targetRange[j+1]:
            dfNode.iloc[i, 2] = j 
dfNode['chr'] = chromosom

# 判断边的类别（是类内的边还是类间的边）
dfEdge['edgeType'] = ['inner' if dfNode[dfNode['name']==row['source']].iloc[0, 2] == dfNode[dfNode['name']==row['target']].iloc[0, 2] else 'inter' for index, row in dfEdge.iterrows()]

edges = dfEdge.to_json(orient='records')
nodes = dfNode.to_json(orient='records')
edgeParsed = json.loads(edges)
nodeParsed = json.loads(nodes)
print(dfEdge)
print(dfNode)
res = '{\n"nodes":' + json.dumps(nodeParsed, indent=4) + ',\n"edges":' + json.dumps(edgeParsed, indent=4) + '\n}'

with open('../result/fragment/%s_%s_%s.json'%(chromosom, str(quantileForEdge), genre), 'w') as f:
    f.write(res)

