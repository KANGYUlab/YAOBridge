from .single_chainmap import *
from .mul_chainmap import *
from .treebuild import *
from .arguments import get_args
# from single_chainmap import *
# from mul_chainmap import *
# from treebuild import *
# from arguments import get_args
import sys
import os
def dedup(items):
    seen = set()
    uniq = []

    for item in items:
        key = tuple(item)

        if key not in seen:
            seen.add(key)
            uniq.append(item)

    return uniq

def main(arglist=None):
    # 1. 获取当前脚本文件的绝对路径（比如：/home/.../yaolinkhg38pcc/yaolinkhg38chain/hg38toyaochainlift.py）
    # # 获取当前文件所在目录的上一级（即项目根目录）
    # BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # 获取当前 PB.py 所在的目录，即 .../site-packages/src/
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    # 建立区间树
    single_yao_2_hg38_trees,single_hg38_2_yao_trees = build_accurate_bidirectional_trees(BASE_DIR+"/data/map1v1.pb.chain") 
    mul_yao_2_hg38_trees, mul_hg38_2_yao_trees = build_accurate_bidirectional_trees(BASE_DIR+"/data/map1vn.pb.chain")
    # 预设的染色体映射（建议放在脚本顶部）
    args = get_args(arglist)
    #判断是否输入文件，输出文件
    if not args.file:
        sys.exit("Error: please provide input file #chrom #pos ******")
    if not args.out:
        sys.exit("Error: please provide output file path")
    # 1. 获取 args.out 的绝对路径并提取其所在的目录
    output_dir = os.path.dirname(os.path.abspath(args.out))

    # 2. 拼接新的文件名
    failed_file_path = os.path.join(output_dir, "failed.txt")
    if args.s=="hg38" and args.t=="yao":
        with open(args.file, 'r') as f, open(args.out, 'w') as f2, open(failed_file_path, 'w') as f3:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith("#"): continue
                data = line.split()
                if len(data) < 2: continue
                ans1v1 = None
                ans1vn = None
                # 使用 try-except 处理映射中可能不存在的染色体键
                try:
                    # chrom_code = hg38_chrom_rev.get(data[0])
                    # if not chrom_code: continue
                    chrom_code = data[0]
                    pos = int(data[1])-1
                    ans1v1 = single_ref_to_query_pos(single_hg38_2_yao_trees[chrom_code], chrom_code, pos)
                    if ans1v1:
                        if ans1v1[0][-3]=="+":
                            try:
                                ans1v1[0][1] = str(int(ans1v1[0][1])+1)
                            except:
                                pass
                        else:
                            try:
                                ans1v1[0][1] = str(int(ans1v1[0][1])-1)
                            except:
                                pass
                        # ans1v1 = dedup(ans1v1) if ans1v1 else None
                    if ans1v1 and ans1v1[0][-1] == "1":
                        pass
                    else:
                        ans1vn = mul_ref_to_query_pos(mul_hg38_2_yao_trees[chrom_code], chrom_code, pos)
                        for item in ans1vn:
                            if item[-3] == "+":
                                try:
                                    item[1] = str(int(item[1])+1)
                                except:
                                    pass
                            else:
                                try:
                                    item[1] = str(int(item[1])-1)
                                except:
                                    pass
                        # ans1vn = dedup(ans1vn) if ans1vn else None
                    # str_1v1 = ",".join(map(str, ans1v1[0])) if ans1v1 else "None"
                    # str_1vn = ";".join([",".join(map(str, item)) for item in ans1vn]) if ans1vn else "None"
                    # 在处理完 ans1v1 和 ans1vn 的偏移逻辑后
                    combined_ans = []

                    if ans1v1:
                        combined_ans.extend(ans1v1)
                    if ans1vn:
                        combined_ans.extend(ans1vn)

                    # 统一进行去重处理
                    final_ans = dedup(combined_ans) if combined_ans else None
                    if final_ans:
                        for item in final_ans:
                            datalist=[str(chrom_code), str(pos+1), *item, *data[2:]]
                            # 在写入之前进行转换
                            f2.write('\t'.join(map(str, datalist)) + '\n')
                    else:
                        f3.write(line + '\n')
                except Exception:
                    continue
    if args.s=="yao" and args.t=="hg38":
        with open(args.file, 'r') as f, open(args.out, 'w') as f2, open(failed_file_path, 'w') as f3:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith("#"): continue
                data = line.split()
                if len(data) < 2: continue
                ans1v1 = None
                ans1vn = None
                # 使用 try-except 处理映射中可能不存在的染色体键
                try:
                    # chrom_code = hg38_chrom_rev.get(data[0])
                    # if not chrom_code: continue
                    chrom_code = data[0]
                    pos = int(data[1])-1
                    ans1v1 = single_query_to_ref_pos(single_yao_2_hg38_trees[chrom_code], chrom_code, pos)
                    if ans1v1:
                        if ans1v1[0][-3]=="+":
                            try:
                                ans1v1[0][1] = str(int(ans1v1[0][1])+1)
                            except:
                                pass
                        else:
                            try:
                                ans1v1[0][1] = str(int(ans1v1[0][1])-1)
                            except:
                                pass
                        # ans1v1 = dedup(ans1v1) if ans1v1 else None
                    if ans1v1 and ans1v1[0][-1] == "1":
                        pass
                    else:
                        ans1vn = mul_query_to_ref_pos(mul_yao_2_hg38_trees[chrom_code], chrom_code, pos)
                        for item in ans1vn:
                            if item[-3] == "+":
                                try:
                                    item[1] = str(int(item[1])+1)
                                except:
                                    pass
                            else:
                                try:
                                    item[1] = str(int(item[1])-1)
                                except:
                                    pass
                    #     ans1vn = dedup(ans1vn) if ans1vn else None
                    # str_1v1 = ",".join(map(str, ans1v1[0])) if ans1v1 else "None"
                    # str_1vn = ";".join([",".join(map(str, item)) for item in ans1vn]) if ans1vn else "None"
                    
                    # datalist = [str(chrom_code), str(pos+1), 'Found' if (ans1v1 or ans1vn) else 'NotFound', str_1v1, str_1vn, *data[2:]]
                    # f2.write('\t'.join(datalist) + '\n')
                    combined_ans = []

                    if ans1v1:
                        combined_ans.extend(ans1v1)
                    if ans1vn:
                        combined_ans.extend(ans1vn)

                    # 统一进行去重处理
                    final_ans = dedup(combined_ans) if combined_ans else None
                    if final_ans:
                        for item in final_ans:
                            datalist=[str(chrom_code), str(pos+1), *item, *data[2:]]
                            # 在写入之前进行转换
                            f2.write('\t'.join(map(str, datalist)) + '\n')
                    else:
                        f3.write(line + '\n')
                except Exception:
                    continue
if __name__ == '__main__':
    main(arglist = [
    "-s", "yao",
    "-t", "hg38",
    "-file", "/home/lfszxyy/old/annotation/gffConfidence-Total/teachersrc/githubPB/src/data/yao2hg38.test",
    "-out", "/home/lfszxyy/old/annotation/gffConfidence-Total/teachersrc/githubPB/src/data/yao2hg38answer.txt"
])
    # main()
#def chrtrans(,chr):
    
