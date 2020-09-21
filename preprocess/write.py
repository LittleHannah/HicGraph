import cooler
import json
import pandas as pd

project = 'Rao'
compartmentFile = '/store/qzhong/test/hanhan/Rao2014-GM12878-MboI-allreps-filtered.compartments.bed.cis.    vecs.tsv'
TADFile = '/store/qzhong/test/hanhan/Rao2014-GM12878-MboI-allreps-filtered.insulation_boundaries.bed'
mcoolFile = '/store/qzhong/test/hanhan/Rao2014-GM12878-MboI-allreps-filtered.mcool::resolutions/10000'
chromosom = 'chr2'
start = 12800000
end = 50690000

# 写compartment和TAD文件
def writeTarget(project, compartmentFile, TADFile, mcoolFile, chromosom, start, end, genre, quantileForEdge):
    # 读取分割文件的内容，并且把分割点写入一个数组中
    inputFile = 'result/split/%s_%s.bed'%(project, genre) 
    target = pd.read_csv(inputFile, sep='\t')
    targetRange = [target[target['chrom']==chromosom].iloc[0,1]] + target[target['chrom']==chromosom].iloc[:,2].tolist()
    print(targetRange)
    targetRange = [i for i in targetRange if i>=start and i<=end]
    print(targetRange)   
    # 初始化 
    dfEdge = pd.DataFrame()
    dfNode = pd.DataFrame()
    source = []
    target = []
    weight = []
    node = []
    value = []
    name = []
    category = []
    tad1 = ''
    tad2 = ''
    
    # 读取mcool文件
    c = cooler.Cooler(mcoolFile)
    for i in range(0, len(targetRange)-1):
        for j in range(0, len(targetRange)-1):
            t1 = '-'.join([str(targetRange[i]), str(targetRange[i+1])])
            t2 = '-'.join([str(targetRange[j]), str(targetRange[j+1])])
            df = c.matrix(balance=False, as_pixels=True, join=True).fetch('%s:%s' % (chromosom, t1), '%s:%s' % (chromosom, t2))
            if j<i:
                pass
            elif j==i:
                name.append(t1)
                value.append(df['count'].mean())
                #category.append(i)
            else:
                source.append(t1)
                target.append(t2)
                weight.append(df['count'].mean())
    dfEdge['source'] = source 
    dfEdge['target'] = target
    dfNode['name'] = name
    if len(weight) < 2:
        dfEdge['weight'] = weight
        dfNode['value'] = value
    else:
        dfEdge['weight'] = [(i-min(weight))/(max(weight)-min(weight)) for i in weight]
        dfNode['value'] = [(i-min(value))/(max(value)-min(value)) for i in value]
    dfNode['category']=0
    for i in range(0, len(dfNode)):
        node = int(dfNode.iloc[i, 0].split('-')[0])
        for j in range(0, len(targetRange)-1):
            if node >= targetRange[j] and node<targetRange[j+1]:
                dfNode.iloc[i, 2] = j
    dfEdge = dfEdge[dfEdge['weight'] > dfEdge['weight'].quantile([quantileForEdge])[quantileForEdge]]
    dfNode['chr'] = chromosom
    print(dfEdge)
    print(dfNode)

    # 将数据框写成json文件
    edges = dfEdge.to_json(orient='records')
    nodes = dfNode.to_json(orient='records')
    edgeParsed = json.loads(edges)
    nodeParsed = json.loads(nodes)
    res = '{\n"nodes":' + json.dumps(nodeParsed, indent=4) + ',\n"edges":' + json.dumps(edgeParsed, indent=4) + '\n    }'
    with open('result/%s/%s_%s_%s.json'%(genre, chromosom, quantileForEdge, genre), 'w') as f:
        f.write(res)

# 写fragment文件
def writeFragment(project, compartmentFile, TADFile, mcoolFile, chromosom, start, end, genre, quantileForEdge):
    # 读取分割文件并将分割点写入一个数组
    inputFile = 'result/split/%s_%s.bed'%(project, genre)
    target = pd.read_csv(inputFile, sep='\t')
    targetRange = [target[target['chrom']==chromosom].iloc[0,1]] + target[target['chrom']==chromosom].iloc[:,2].tolist()
    targetRange = [i for i in targetRange if i>=start and i<=end]

    chromRange = '%s:%s-%s'%(chromosom, str(targetRange[0]), str(targetRange[-1]))
    #genre = 'fragment_compartment'
    
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
    print('haha')
    
    with open('result/%s/fragment_%s_%s.json'%(genre, chromosom, str(quantileForEdge)), 'w') as f:
        f.write(res)
