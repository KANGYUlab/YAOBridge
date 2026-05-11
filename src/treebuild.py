import collections
from intervaltree import Interval, IntervalTree
def build_accurate_bidirectional_trees(chain_file):
    yao_to_hg38_trees = collections.defaultdict(IntervalTree)
    hg38_to_yao_trees = collections.defaultdict(IntervalTree)

    current_chain = None
    current_data_blocks = []

    def push_to_tree():
        """结算并存入树的内部函数"""
        if current_chain and current_data_blocks:
            # 存入 YAO -> hg38 树 (YAO 坐标作为 Key)
            yao_to_hg38_trees[current_chain['yao_name']].add(
                Interval(current_chain['yao_start'], current_chain['yao_end'], {
                    'hg38_name': current_chain['hg38_name'],
                    'hg38_start': current_chain['hg38_start'],
                    'hg38_end': current_chain['hg38_end'],
                    'strand': current_chain['strand'],
                    'score': current_chain['score'],
                    'level': current_chain['level'],
                    'cigar':current_data_blocks
                })
            )
            # 存入 hg38 -> YAO 树 (hg38 坐标作为 Key)
            hg38_to_yao_trees[current_chain['hg38_name']].add(
                Interval(current_chain['hg38_start'], current_chain['hg38_end'], {
                    'yao_name': current_chain['yao_name'],
                    'yao_start': current_chain['yao_start'],
                    'yao_end': current_chain['yao_end'],
                    'strand': current_chain['strand'],
                    'score': current_chain['score'],
                    'level': current_chain['level'],
                    'cigar':current_data_blocks
                })
            )

    with open(chain_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue

            if line.startswith('chain'):
                # 1. 遇到新 chain，先结算上一个
                push_to_tree()
                
                # 2. 解析新 chain 的头部
                parts = line.strip().split('\t')
                # 你的格式：第3列是 YAO (T)，第8列是 hg38 (Q)
                current_chain = {
                    'yao_name': parts[1],
                    'strand': parts[3],
                    'yao_start': int(parts[4]),
                    'yao_end': int(parts[5]),
                    'hg38_name': parts[6],
                    'hg38_strand': '+',
                    'hg38_start': int(parts[9]),
                    'hg38_end': int(parts[10]),
                    'score': float(parts[11]),
                    'level': parts[12]
                }
                current_data_blocks = []
            else:
                # 3. 它是数据行，直接追加
                current_data_blocks.append(line.split('\t'))

        # 4. 别忘了结算最后一个 chain
        push_to_tree()

    return yao_to_hg38_trees, hg38_to_yao_trees