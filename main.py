import argparse
from preprocess import split 
from preprocess import write
import warnings
warnings.filterwarnings('ignore')


parser = argparse.ArgumentParser(description = '-----------------Using HicGraph-----------------------')
parser.add_argument('-o', '--output', required = True, help='specify the output directory name')
parser.add_argument('--compartment', required=True, help='specify the directory of compartment file')
parser.add_argument('--tad', required=True, help='specify the directory of TAD file')
parser.add_argument('--mcool', required=True, help='specify the directory of mcool file')
parser.add_argument('-c', '--chr', required=True, help='the chromosome you want to visulize')
parser.add_argument('--start', required=True, help='start position of the region you want to visulize')
parser.add_argument('--end', required=True, help='end position of the region you want to visulize')
parser.add_argument('--genre',default = 'fragment_tad',  help='sepcify the type of graph you want to generate')
parser.add_argument('-q',  '--quantile', default=0, help='sepcify the quantile of edges you want to remove')
parser.add_argument('--raw', action='store_true', help='if the compartment file and TAD file are raw files obtained from cooltools, please choose this parameter')

args = parser.parse_args()
project = args.output
compartmentFile = args.compartment
TADFile = args.tad
mcoolFile = args.mcool
chromosom = args.chr
start = int(args.start)
end = int(args.end)
genre = args.genre
quantileForEdge = float(args.quantile)
raw = args.raw
print(args)

exit()

project = 'Rao' 
compartmentFile = '/store/qzhong/test/hanhan/Rao2014-GM12878-MboI-allreps-filtered.compartments.bed.cis.vecs.tsv'
TADFile = '/store/qzhong/test/hanhan/Rao2014-GM12878-MboI-allreps-filtered.insulation_boundaries.bed'
mcoolFile = '/store/qzhong/test/hanhan/Rao2014-GM12878-MboI-allreps-filtered.mcool::resolutions/10000'
chromosom = 'chr2'
start = 12800000
end = 50690000
quantileForEdge = 0.99
genre = 'compartment'

if raw == True:
    split.splitCompartment(compartmentFile, 'result/split/%s_compartment.bed'%project, 22)
    split.splitTAD(TADFile, 'result/split/%s_TAD.bed'%project)
write.writeTarget(project, compartmentFile, TADFile, mcoolFile, chromosom, start, end, genre, quantileForEdge)
write.writeFragment(project, compartmentFile, TADFile, mcoolFile, chromosom, start, end, genre, quantileForEdge)
