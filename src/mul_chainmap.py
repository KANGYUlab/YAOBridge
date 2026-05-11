def mul_ref_to_query_pos(tree,chrom,pos):
    yao_pos=[]
    # overlap 返回的是一个 set，包含所有覆盖该位点的 Interval 对象
    # for tree in trees.values():
    #     hits = tree.at(pos) # 或者用 .overlap(pos, pos+1)
    hits = tree.at(pos)
    if not hits:
        print(f"No chains found at position {pos+1}")
    # print(hits)
    for interval in hits:
        # interval.begin: 该 chain 的起始
        # interval.end: 该 chain 的结束
        # interval.data: 你存进去的字典 (包含 cigar, dest_start 等)
        # print(f"Found chain: {interval.begin}-{interval.end}\n{interval.data}")
        chain_info = interval.data
        hg38_start = interval.begin
        hg38_end = interval.end
        yao_name = interval.data['yao_name']
        yao_start = interval.data['yao_start']
        yao_end = interval.data['yao_end']
        strand = interval.data['strand']
        level = interval.data['level']
        cigar = interval.data['cigar']
# 初始化双指针
        # r_curr 指向 hg38 (Reference)，q_curr 指向 YAO (Query)
        r_curr = hg38_start
        q_curr = yao_start  # 始终按正向累加进度

        # 开始遍历 CIGAR 块
        for block in cigar:
            size = int(block[0])
            
            # 1. 命中 Match 块：检查输入的 pos 是否落在当前的 hg38 块内
            if r_curr <= pos < r_curr + size:
                offset = pos - r_curr
                # 计算正向进度下的 YAO 位置
                q_pos_forward = q_curr + offset
                
                if strand == '+':
                    final_yao_pos = q_pos_forward
                else:
                    # 负链翻转：用 YAO 区间的总长度减去进度偏移
                    # 公式：yao_end - (q_pos_forward - yao_start) - 1
                    final_yao_pos = yao_end - (q_pos_forward - yao_start) - 1
                
                # print(f"RESULT: hg38:{pos} -> YAO:{final_yao_pos} ({strand})")
                yao_pos.append([yao_name,final_yao_pos, strand, chain_info['score'],level])
                break # 找到即返回

            # 2. 未命中 Match，同步更新两个指针的进度
            r_curr += size
            q_curr += size

            # 3. 处理 Gap 逻辑
            if len(block) == 3:
                dq = int(block[1])  # YAO 的 gap (dq)
                dt = int(block[2])  # hg38 的 gap (dt)
                
                # 如果 pos 落在 hg38 的 gap (dt) 区域，说明 YAO 这边是插入，hg38 无法对应
                if r_curr <= pos < r_curr + dt:
                    offset = 0
                    # 计算正向进度下的 YAO 位置
                    q_pos_forward = q_curr + offset
                    
                    if strand == '+':
                        final_yao_pos = q_pos_forward
                    else:
                        # 负链翻转：用 YAO 区间的总长度减去进度偏移
                        # 公式：yao_end - (q_pos_forward - yao_start) - 1
                        final_yao_pos = yao_end - (q_pos_forward - yao_start) - 1
                    # print(f"Position {pos} falls in hg38 gap (Deletion in YAO).")
                    yao_pos.append([yao_name,str(final_yao_pos)+"Gap", strand, chain_info['score'],level])
                    break # 找到即返回
                
                # 跳过 Gap
                r_curr += dt
                q_curr += dq
            else:
                # 最后一个 block，没有 gap
                break
                
    return yao_pos
def mul_query_to_ref_pos(tree,chrom,pos):
    hg38_pos=[]
    # overlap 返回的是一个 set，包含所有覆盖该位点的 Interval 对象
    # for tree in trees.values():
    #     hits = tree.at(pos) # 或者用 .overlap(pos, pos+1)
    hits = tree.at(pos)
    if not hits:
        print(f"No chains found at position {pos+1}")
    for interval in hits:
        # interval.begin: 该 chain 的起始
        # interval.end: 该 chain 的结束
        # interval.data: 你存进去的字典 (包含 cigar, dest_start 等)
        # print(f"Found chain: {interval.begin}-{interval.end}\n{interval.data}")
        chain_info = interval.data
        hg38_name = interval.data['hg38_name']
        hg38_start = interval.data['hg38_start']
        hg38_end = interval.data['hg38_end']
        yao_start = interval.begin
        yao_end = interval.end
        strand = interval.data['strand']
        level = interval.data['level']
        cigar = interval.data['cigar']
        # if strand == '+':
        #     target_q = pos
        # else:
        #     target_q = yao_start+(yao_end-pos-1)
        # for block in cigar:
        #     qstart=yao_start
        #     rstart=hg38_start
        #     size,Indel,Del = int(block[0]),int(block[1]),int(block[2])
    # 1. 坐标归一化：把“倒着走”的负链坐标换算成“从头开始数”的进度
        if strand == '+':
            target_q = pos
        else:
            # 负链公式：起点 + (终点 - 1 - 当前位点)
            target_q = yao_start + (yao_end - pos - 1)

        # 2. 初始化双指针（必须在循环外，记录当前匹配到了哪里）
        curr_q = yao_start
        curr_r = hg38_start

        # 3. 遍历 block
        for block in cigar:
            # 解压数据：size(匹配长度), Indel(YAO gap), Del(hg38 gap)
            size = int(block[0])
            
            # --- 判定命中：看 target_q 是否落在当前的匹配块内 ---
            if curr_q <= target_q < curr_q + size:
                offset = target_q - curr_q
                
                # 如果是正链，hg38 也是往后加；如果是负链，hg38 实际上是往前走
                # 但如果你在这里处理的是绝对坐标映射，通常返回 curr_r + offset
                final_pos = curr_r + offset
                # print(f"YAO:{pos} -> hg38:{final_pos}")
                hg38_pos.append([hg38_name,final_pos, strand, chain_info['score'],level])
                break # 找到即返回

            # --- 更新进度：跳过当前的 Match 块 ---
            curr_q += size
            curr_r += size

            # --- 处理 Gap ---
            if len(block) == 3:
                dq = int(block[1])  # YAO 的 gap
                dt = int(block[2])  # hg38 的 gap
                
                # 如果 target_q 落在 YAO 的 gap 区域，说明 hg38 这边是缺失
                if curr_q <= target_q < curr_q + dq:
                    offset = 0
                    
                    # 如果是正链，hg38 也是往后加；如果是负链，hg38 实际上是往前走
                    # 但如果你在这里处理的是绝对坐标映射，通常返回 curr_r + offset
                    final_pos = curr_r + offset
                    #print(f"Position {pos} falls in hg38 gap (Deletion in YAO).")
                    #print(f"YAO:{pos} 落在{curr_q} Gap 区域，无 hg38 对应坐标")
                    hg38_pos.append([hg38_name,str(final_pos)+"Gap", strand, chain_info['score'],level])
                    break # 找到即返回
                
                # 更新指针，跳过 Gap
                curr_q += dq
                curr_r += dt
            else:
                # 最后一个 block，没有 gap
                break
    return hg38_pos            