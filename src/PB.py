from .single_chainmap import *
from .mul_chainmap import *
from .treebuild import *
from .single_bedmap import *
from .mul_bedmap import *
from .arguments import get_args
# from single_chainmap import *
# from mul_chainmap import *
# from treebuild import *
# from arguments import get_args
import sys
import os
import gzip

def dedup(items):
    seen = set()
    uniq = []

    for item in items:
        key = tuple(item)

        if key not in seen:
            seen.add(key)
            uniq.append(item)

    return uniq

def clean_pos(pos):
    if isinstance(pos, int):
        return pos
    # 如果是字符串，去掉 Gap 后转成 int
    return int(pos.replace("Gap", ""))

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
    if args.command == 'liftpos':
        #判断是否输入文件，输出文件
        if not args.file:
            sys.exit("Error: please provide input file #chrom #pos ******")
        if not args.out:
            sys.exit("Error: please provide output file path")
        # 1. 强制自动创建 args.out 里面包含的所有父文件夹（防止 FileNotFoundError）
        output_dir = os.path.dirname(os.path.abspath(args.out))
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        # 2. 安全切分后缀，不管它是 .vcf 还是 .tsv，都能在文件名后面加上 _failed.txt
        file_base, _ = os.path.splitext(args.out)
        failed_file_path = f"{file_base}_failed.txt"
        # 3. 动态判断并打开输入文件（兼容 .gz 和普通文本）
        is_gzip = args.file.endswith('.gz')
        open_func = gzip.open if is_gzip else open
        # 如果是 gz 文件，需要用 'rt' (read text) 模式读取字符串，普通文件用 'r'
        open_mode = 'rt' if is_gzip else 'r'

        if args.s == "hg38" and args.t == "yao":
            with open_func(args.file, open_mode) as f, open(args.out, 'w') as f2, open(failed_file_path, 'w') as f3:
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
                            try:
                                ans1v1[0][1] = str(int(ans1v1[0][1])+1)
                            except:
                                pass
                        if ans1v1 and ans1v1[0][-1] == "1":
                            pass
                        else:
                            ans1vn = mul_ref_to_query_pos(mul_hg38_2_yao_trees[chrom_code], chrom_code, pos)
                            for item in ans1vn:
                                try:
                                    item[1] = str(int(item[1])+1)
                                except:
                                    pass
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
                            try:
                                ans1v1[0][1] = str(int(ans1v1[0][1])+1)
                            except:
                                pass
                        if ans1v1 and ans1v1[0][-1] == "1":
                            pass
                        else:
                            ans1vn = mul_query_to_ref_pos(mul_yao_2_hg38_trees[chrom_code], chrom_code, pos)
                            for item in ans1vn:
                                try:
                                    item[1] = str(int(item[1])+1)
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
    if args.command == 'liftbed':
        #判断是否输入文件，输出文件
        if not args.file:
            sys.exit("Error: please provide input file #chrom #start #end ******")
        if not args.out:
            sys.exit("Error: please provide output file path")
        # 1. 强制自动创建 args.out 里面包含的所有父文件夹（防止 FileNotFoundError）
        output_dir = os.path.dirname(os.path.abspath(args.out))
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        # 2. 安全切分后缀，不管它是 .vcf 还是 .tsv，都能在文件名后面加上 _failed.txt
        file_base, _ = os.path.splitext(args.out)
        failed_file_path = f"{file_base}_failed.txt"
        # 3. 动态判断并打开输入文件（兼容 .gz 和普通文本）
        is_gzip = args.file.endswith('.gz')
        open_func = gzip.open if is_gzip else open
        # 如果是 gz 文件，需要用 'rt' (read text) 模式读取字符串，普通文件用 'r'
        open_mode = 'rt' if is_gzip else 'r'

        if args.s == "hg38" and args.t == "yao":
            with open_func(args.file, open_mode) as f, open(args.out, 'w') as f2, open(failed_file_path, 'w') as f3:
                lines=f.readlines()
                for line in lines:
                    line=line.strip()
                    if line.startswith("#"):
                        continue
                    cols=line.split("\t")
                    chrom=cols[0]
                    start=int(cols[1])-1
                    end=int(cols[2])-1
                    bed=[start,end]
                    datalist=cols[3:]
                    single_hits=single_hg38_2_yao_trees[chrom].overlap(*bed)
                    # for hit in hits:
                    #     print(hit)
                    mul_hits=mul_hg38_2_yao_trees[chrom].overlap(*bed)
                    if not single_hits and not mul_hits:
                        f3.write(line + '\n')
                        continue
                    pos_list=[]
                    for hit in single_hits: 
                        yao_pos=single_bed_ref_2_query(bed,hit)
                        pos_list.extend(yao_pos)
                    if len(pos_list)==2 and pos_list[0][0]!=pos_list[1][0]:
                        try:
                            pos1=pos_list[0][2]+1
                        except:
                            pos1=pos_list[0][2]
                        try:
                            pos2=pos_list[1][2]+1
                        except:
                            pos2=pos_list[1][2]
                        bedstrand=pos_list[0][3]
                        newchrom=pos_list[0][1]
                        score=str(pos_list[0][4])
                        level=pos_list[0][5]
                        if clean_pos(pos1)<clean_pos(pos2):
                            f2.write("\t".join([chrom,str(start+1),str(end+1),newchrom,str(pos1),str(pos2),bedstrand,score,level,*datalist])+"\n")  
                        else:
                            f2.write("\t".join([chrom,str(start+1),str(end+1),newchrom,str(pos2),str(pos1),bedstrand,score,level,*datalist])+"\n")  
                    else:
                        f3.write("\t".join([chrom,str(start+1),str(end+1),newchrom,str(yao_pos),score,level,*datalist])+"\n")

                    for hit in mul_hits:
                        yao_pos=mul_bed_ref_2_query(bed,hit)
                        if len(yao_pos)==2 and yao_pos[0][0]!=yao_pos[1][0]:
                            newchrom=yao_pos[0][1]
                            try:
                                pos1=yao_pos[0][2]+1
                            except:
                                pos1=yao_pos[0][2]
                            try:
                                pos2=yao_pos[1][2]+1
                            except:
                                pos2=yao_pos[1][2]
                            bedstrand=yao_pos[0][3]
                            score=str(yao_pos[0][4])
                            level=yao_pos[0][5]
                            if type(pos1) == str or type(pos2) == str:
                                f3.write("\t".join([chrom,str(start+1),str(end+1),newchrom,str(pos1),str(pos2),bedstrand,score,level,*datalist])+"\n")  
                            else:
                                if pos1 < pos2:
                                    f2.write("\t".join([chrom,str(start+1),str(end+1),newchrom, str(pos1),str(pos2),bedstrand,score,level,*datalist])+"\n")  
                                else:
                                    f2.write("\t".join([chrom,str(start+1),str(end+1),newchrom, str(pos2),str(pos1),bedstrand,score,level,*datalist])+"\n")  
        if args.s=="yao" and args.t=="hg38":
            with open(args.file, 'r') as f, open(args.out, 'w') as f2, open(failed_file_path, 'w') as f3:
                lines=f.readlines()
                for line in lines:
                    line=line.strip()
                    if line.startswith("#"):
                        continue
                    cols=line.split("\t")
                    chrom=cols[0]
                    start=int(cols[1])-1
                    end=int(cols[2])-1
                    datalist=cols[3:]
                    bed=[start,end]
                    single_hits=single_yao_2_hg38_trees[chrom].overlap(*bed)
                    # for hit in hits:
                    #     print(hit)
                    mul_hits=mul_yao_2_hg38_trees[chrom].overlap(*bed)
                    if not single_hits and not mul_hits:
                        f3.write(line + '\n')
                        continue
                    pos_list=[]
                    for hit in single_hits: 
                        yao_pos=single_bed_query_2_ref(bed,hit)
                        # print(yao_pos)
                        pos_list.extend(yao_pos)
                    if len(pos_list)==2 and pos_list[0][0]!=pos_list[1][0]:
                        try:
                            pos1=yao_pos[0][2]+1
                        except:
                            pos1=yao_pos[0][2]
                        try:
                            pos2=yao_pos[1][2]+1
                        except:
                            pos2=yao_pos[1][2]
                        bedstrand=pos_list[0][3]
                        newchrom=pos_list[0][1]
                        score=str(pos_list[0][4])
                        level=pos_list[0][5]
                        if clean_pos(pos1)<clean_pos(pos2):
                            f2.write("\t".join([chrom,str(start+1),str(end+1),newchrom,str(pos1),str(pos2),bedstrand,score,level,*datalist])+"\n")  
                        else:
                            f2.write("\t".join([chrom,str(start+1),str(end+1),newchrom,str(pos2),str(pos1),bedstrand,score,level,*datalist])+"\n")  
                    else:
                        f3.write("\t".join([chrom,str(start+1),str(end+1),*datalist])+"\n")

                    for hit in mul_hits:
                        yao_pos=mul_bed_query_2_ref(bed,hit)
                        if len(yao_pos)==2 and yao_pos[0][0]!=yao_pos[1][0]:
                            newchrom=yao_pos[0][1]
                            try:
                                pos1=yao_pos[0][2]+1
                            except:
                                pos1=yao_pos[0][2]
                            try:
                                pos2=yao_pos[1][2]+1
                            except:
                                pos2=yao_pos[1][2]
                            bedstrand=yao_pos[0][3]
                            score=str(yao_pos[0][4])
                            level=yao_pos[0][5]
                            # if type(pos1) == str or type(pos2) == str:
                            if len(yao_pos)==2:
                                if pos1 < pos2:
                                    f2.write("\t".join([chrom,str(start+1),str(end+1),newchrom, str(pos1),str(pos2),bedstrand,score,level,*datalist])+"\n")  
                                else:
                                    f2.write("\t".join([chrom,str(start+1),str(end+1),newchrom, str(pos2),str(pos1),bedstrand,score,level,*datalist])+"\n")  
                            else:
                                f3.write("\t".join([chrom,str(start+1),str(end+1),*datalist])+"\n") 



if __name__ == '__main__':
    main()
    # main(arglist = [
    # "liftbed",
    # "-s", "hg38",
    # "-t", "yao",
    # "-file", "/home/lfszxyy/old/annotation/gffConfidence-Total/teachersrc/speedhmPB/YAOBridge/src/data/hg38bed.test",
    # "-out", "/home/lfszxyy/old/annotation/gffConfidence-Total/teachersrc/speedhmPB/YAOBridge/src/data/hg382yaobed.test"
    # ])
    # main(arglist = [
    # "liftbed",
    # "-s", "hg38",
    # "-t", "yao",
    # "-file", "/home/lfszxyy/old/annotation/gffConfidence-Total/teachersrc/speedhmPB/YAOBridge/src/beddata/hg38bed.test",
    # "-out", "/home/lfszxyy/old/annotation/gffConfidence-Total/teachersrc/speedhmPB/YAOBridge/src/beddata/hg382yaobed.test"
    # ])

