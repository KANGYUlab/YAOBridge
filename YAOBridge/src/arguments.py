# import argparse
# import sys
# from src import __version__

# def get_args(arglist):
#     parser = argparse.ArgumentParser(description='PreciseBridge: hg38 <-> yao mapping')
    
#     # 1. 转换方向 (从哪到哪)
#     parser.add_argument('-s', choices=['hg38', 'yao'], required=True, help='source Genome')
#     parser.add_argument('-t', choices=['hg38', 'yao'], required=True, help='target Genome')
    
#     # 批量位置文件：每行格式预计为 "chr pos"
#     parser.add_argument('-file', type=str, required=True, help='#chrom #pos #---')

#     # 输出文件位置
#     parser.add_argument('-out', type=str, required=True, help='output file path')
#     parser.add_argument('-V', '--version', help='show program version', action='version', version=__version__)
#     args = parser.parse_args(arglist)
#     if args.s == args.t:
#         sys.exit("Error: source and target cannot be the same")
#     # return parser.parse_args()
#     # 关键修改：如果传入了 arglist，就解析它；否则解析 sys.argv
#     return parser.parse_args(arglist)


import argparse
import sys
# 假设你的 __version__ 正常导入
__version__ = "2.0.0" 

def get_args(arglist):
    # 1. 创建主解析器
    parser = argparse.ArgumentParser(description='PreciseBridge: hg38 <-> yao coordinate mapping tool')
    parser.add_argument('-V', '--version', help='show program version', action='version', version=__version__)
    
    # 2. 创建子命令管理器 (dest='command' 可以用来判断用户到底用了哪一个功能)
    subparsers = parser.add_subparsers(dest='command', required=True, help='Available commands')

    # ==========================================
    # 功能一：liftpos 的参数
    # ==========================================
    parser_pos = subparsers.add_parser('liftpos', help='Mapping for positions file (chr pos(1-base))')
    parser_pos.add_argument('-s', choices=['hg38', 'yao'], required=True, help='source Genome')
    parser_pos.add_argument('-t', choices=['hg38', 'yao'], required=True, help='target Genome')
    parser_pos.add_argument('-file', type=str, required=True, help='input file path with #chrom #pos')
    parser_pos.add_argument('-out', type=str, required=True, help='output file path')

    # ==========================================
    # 功能二：liftbed 的参数（这里为你定制新的 BED 参数）
    # ==========================================
    parser_bed = subparsers.add_parser('liftbed', help='Mapping for BED files (chr start end(1-base))')
    parser_bed.add_argument('-s', choices=['hg38', 'yao'], required=True, help='source Genome')
    parser_bed.add_argument('-t', choices=['hg38', 'yao'], required=True, help='target Genome')
    # 比如 BED 独有的参数：输入 bed，输出 bed
    parser_bed.add_argument('-file', type=str, required=True, help='input file path with #chrom #start #end')
    parser_bed.add_argument('-out', type=str, required=True, help='output file path')
    # 甚至可以加一些 BED 专属的质控过滤参数
    parser_bed.add_argument('--min-overlap', type=float, default=0.95, help='')

    # 3. 解析参数
    args = parser.parse_args(arglist)
    
    # 4. 统一的源和目标防呆检查
    if args.s == args.t:
        sys.exit("Error: source and target cannot be the same")
        
    return args

# ==========================================
# 下游如何根据参数分流调用？（主函数伪代码示例）
# ==========================================
# if __name__ == '__main__':
#     # 接收命令行参数
#     args = get_args(sys.argv[1:])
    
#     # 根据 args.command 来判断用户到底敲的是哪个功能
#     if args.command == 'liftpos':
#         print(f"正在运行单点转换：从 {args.s} 到 {args.t}，读取文件 {args.file}")
#         # run_liftpos(args)
        
#     elif args.command == 'liftbed':
#         print(f"正在运行 BED 转换：从 {args.s} 到 {args.t}，读取 BED {args.bed}，最小重叠度 {args.min_overlap}")
#         # run_liftbed(args)