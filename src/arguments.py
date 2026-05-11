import argparse
import sys
from src import __version__

def get_args(arglist):
    parser = argparse.ArgumentParser(description='PreciseBridge: hg38 <-> yao mapping')
    
    # 1. 转换方向 (从哪到哪)
    parser.add_argument('-s', choices=['hg38', 'yao'], required=True, help='source Genome')
    parser.add_argument('-t', choices=['hg38', 'yao'], required=True, help='target Genome')
    
    # 批量位置文件：每行格式预计为 "chr pos"
    parser.add_argument('-file', type=str, required=True, help='#chrom #pos #---')

    # 输出文件位置
    parser.add_argument('-out', type=str, required=True, help='output file path')
    parser.add_argument('-V', '--version', help='show program version', action='version', version=__version__)
    args = parser.parse_args(arglist)
    if args.s == args.t:
        sys.exit("Error: source and target cannot be the same")
    # return parser.parse_args()
    # 关键修改：如果传入了 arglist，就解析它；否则解析 sys.argv
    return parser.parse_args(arglist)
