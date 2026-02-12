# Streamlined ETE Toolkit
# Copyright (c) 2026 ZJX

# Debug flag - set to True when running as main module
DEBUG = __name__ == '__main__'

def inspect_newick_string(tree_blocks):
    # inspect the original newick from newick_blocks
    newick_tree = ''
    last_block = -1 # 0: ( 1: taxon 2: )
    for i in tree_blocks:
        if i[0] == '(':
            if last_block in [1,2]:
                newick_tree += ',' # a sister clade with a single taxon (1) or a clade (2)
            newick_tree += i
            last_block = 0
        elif i[0] == ')':
            # No need to add a comma
            newick_tree += i
            last_block = 2
        else: # OTU
            if last_block in [1,2]:
                newick_tree += ','
            newick_tree += i
            last_block = 1
    newick_tree += ';'
    return newick_tree
def build_tree_block(newick_string):
    tree_blocks = ['']
    for i in newick_string: # split newick string into blocks to parse taxon names
        if i == ' ':
            continue
        # ( / ) indicates the increase / decrease of stack depth
        elif i == ')':
            tree_blocks.append(i)
        elif i == '(':
            tree_blocks[-1] += i
            tree_blocks.append('')
        elif i == ',':
            # tree_blocks.append(i)
            tree_blocks.append('')
            # for the next possible OTU or (
        elif i == ';':
            break
        else:
            tree_blocks[-1] += i 
    return tree_blocks

def find_mrca_position(tree_blocks, target_taxa_set):
    """Find the MRCA position of a set of taxa in a newick string
    
    :param newick_string: the newick tree
    :param target_taxa_set: a set of target taxa
    """
    # if isinstance(DEBUG, bool):
    #     if DEBUG == True:
    #         print("Query:", target_taxa_set)
    # NOTE: We transverse the tree from the ROOT! Therefore, the smaller the stack depth, the closer to the root
    depth = 0
    lower_bound = -1 # lower bound is the lowest stack depth of all clades derived from the target MRCA
    for index, i in enumerate(tree_blocks):
        if i[0] in '()': # Note that ')' may have extra annotations
            if i[0] == '(':
                depth += 1
            elif i[0] == ')':
                depth -= 1
            if target_taxa_set == []:
                if depth < lower_bound:
                    # This is the MRCA!
                    return index
            lower_bound = min(depth, lower_bound) if lower_bound >= 0 else lower_bound
        else:
            i = i.split(':')[0] if ':' in i else i # split the taxon name
            if i in target_taxa_set:
                target_taxa_set.remove(i)
                lower_bound = depth if lower_bound < 0 else lower_bound
        
        if isinstance(DEBUG, bool):
            if DEBUG == True:
                print(f"depth: {depth}\tblock: {i}\tlower_bound:{lower_bound}\tindex: {index}")

def annotate_mrca_position(tree_blocks, target_taxa_set):
    if isinstance(DEBUG, bool):
        if DEBUG == True:
            print("Query:", target_taxa_set)
    # split newick string into blocks
    # tree_blocks = build_tree_block(newick_string)
    
    if isinstance(DEBUG, bool):
        if DEBUG == True:
            print(tree_blocks)
    # NOTE: We transverse the tree from the ROOT! Therefore, the smaller the stack depth, the closer to the root

    mrca_index = find_mrca_position(tree_blocks, target_taxa_set)
    tree_blocks[mrca_index] += "[&MRCA]"

    return tree_blocks
    



    
def query_genetic_distance(tree_blocks, target_taxa_sets):
    """
    Stat genetic distances between multiple paired taxa in a tree
    
    :param tree_blocks: the newick tree block built by function build_tree_block
    :param target_taxa_sets: a list of target taxa sets, e.g. [['A', 'B'], ['C', 'D']]
    """
    # tree_blocks = build_tree_block(newick_string)

    records = [[0, -1, -1, 0] for _ in target_taxa_sets]

    # index 0: genetic distance
    # index 1: status (-1: unstarted, 0: started, 1: ended)
    # index 2: lowest depth (MRCA reconstruction, for branch length statistics of the next half of the tree)
    # index 3: whether we are now at the right branch

    depth = 0

    for index, i in enumerate(tree_blocks):
        if i[0] in '()':
            if i[0] == ')':
                depth -= 1
                for j in records:
                    if j[1] == 0:
                        # j[2] = min(j[2], depth) if j[2] > 0 else depth
                        if j[2] < 0:
                            j[2] = depth
                        elif j[2] > depth:
                            j[2] = depth
                            j[0] = j[0] + float(i.split(':')[1])
                    elif j[1] == 1:
                        if j[3] == 0:
                            if j[2] <= depth and j[2] >= 0: 
                                j[0] += float(i.split(':')[1])
                                if j[2] == depth:
                                    j[2] = -1 # We have reached MRCA. Stop adding genetic distance
                        else:
                            j[3] -= 1

            if i[0] == '(':
                depth += 1
                for j in records:
                    if j[1] == 1:
                        j[3] += 1 # unrelated branch detected. DO NOT ADD TO GENETIC DISTANCE
        else: # OTU Block
            name = i.split(':')[0]
            for jndex, j in enumerate(target_taxa_sets):
                if name in j:
                    # start or end a record
                    if records[jndex][1] == -1:
                        records[jndex][1] = 0
                        records[jndex][0] = float(i.split(':')[1])
                    if records[jndex][1] == 0:
                        records[jndex][1] = 1
                        records[jndex][0] = records[jndex][0] + float(i.split(':')[1])
    
    return [i[0] for i in records]

def extract_subtree(tree_blocks, target_taxa_set):
    """
    Extract a subtree from a newick string
    
    :param tree_blocks: the newick tree block built by function build_tree_block
    :param target_taxa_set: a set of target taxa
    """
    # tree_blocks = build_tree_block(newick_string)
    
    # 1. find MRCA node
    mrca_index = find_mrca_position(tree_blocks, target_taxa_set)

    # 2. reverse parse to extract subtree
    depth = 0
    reversed_tree_blocks = []
    for index in range(mrca_index, -1, -1):
        if tree_blocks[index][0] == ')':
            depth += 1 # NOTE: Reversed order!
            # reversed_tree_blocks.append(tree_blocks[index])
        elif tree_blocks[index][0] == '(':
            depth -= 1
            if depth == 0:
                break
            # reversed_tree_blocks.append(tree_blocks[index])
        reversed_tree_blocks.append(tree_blocks[index])
    
    # reverse the blocks to the original order
    reversed_tree_blocks = reversed_tree_blocks[::-1]
    return reversed_tree_blocks

def rotate_clade(tree_blocks, node_index):
    """
    Rotate a clade in a tree
    
    :param tree_blocks: the newick tree block built by function build_tree_block
    :param index: the index of the clade to be rotated
    """
    # print(tree_blocks, node_index)
    # tree_blocks = build_tree_block(newick_string)
    # 1. find the clade
    depth = 1
    rotated_tree_blocks = tree_blocks[:]
    splitted = False
    tmp_blocks_1, tmp_blocks_2 = [tree_blocks[node_index]], []
    for index in range(node_index-1, -1, -1):
        # Why iterate from index-1:
        # tb[index] is the right parentheses of MRCA of the clade
        # We need to find out the split point of the clade, where relative depth is 1
        # By iterating from index-1, we can avoid preserving an extra state machine
        if tree_blocks[index][0] == ')':
            depth += 1
            if splitted:
                tmp_blocks_2.append(tree_blocks[index])
            else:
                tmp_blocks_1.append(tree_blocks[index])
        elif tree_blocks[index][0] == '(':
            depth -= 1
            if splitted:
                tmp_blocks_2.append(tree_blocks[index])
            else:
                tmp_blocks_1.append(tree_blocks[index])
            if depth == 1:
                splitted = True
            if depth == 0:
                print("ROTATING CLADE - ", tree_blocks[index: node_index+1])
                rotated_clade = ["("] + tmp_blocks_1[1::][::-1] + tmp_blocks_2[-2::-1] + [tmp_blocks_1[0]]
                print("ROTATED CLADE - ",rotated_clade)
                rotated_tree_blocks = tree_blocks[:index] + rotated_clade + tree_blocks[-1:node_index:-1][::-1]
                return rotated_tree_blocks
        else:
            if splitted:
                tmp_blocks_2.append(tree_blocks[index])
            else:
                tmp_blocks_1.append(tree_blocks[index])

def sort_topology(tree_blocks, ascending=True):
    """
    Ascending / Descending sort the topology of a tree
    NOTE: This function is not applicable to trees with unsolved nodes or may generate incorrect results
    
    :param tree_blocks: the newick tree block built by function build_tree_block
    :param ascending: True for ascending, False for descending

    NOTE: There's one critical bug remaining in this function
    """

    depth = 0
    otu_nums = [0]
    temp = 0 # store the number of OTU in the right of a comma
    sorted_tree_blocks = tree_blocks[:]
    print(tree_blocks)
    for index, i in enumerate(tree_blocks):
        if i[0] in '()':
            if i[0] == ')':
                # depth -= 1
                swap = False
                if len(otu_nums) > depth + 1:
                    swap = ((otu_nums[depth +1] !=0 and otu_nums[depth +1] < temp and ascending) or (otu_nums[depth +1] !=0 and otu_nums[depth +1] > temp and not ascending))
                if swap:
                    print("SWAP - index:", index, "OTU_NUM:", otu_nums[depth], "temp:", temp)
                    sorted_tree_blocks = rotate_clade(sorted_tree_blocks, index) 
                    if isinstance(DEBUG, bool) and DEBUG == True:
                        print(inspect_newick_string(sorted_tree_blocks))
                    # # NOTE: The index of remaining nodes won't change, so no extra work
                    temp += otu_nums[depth]
                    otu_nums[depth] = 0
                depth -= 1
                # if otu_nums[depth] < temp:
            elif i[0] == '(':
                depth += 1
                if len(otu_nums) < depth + 1:
                    otu_nums.append(temp)
                    temp = 0
                else:
                    otu_nums[depth] = temp
                    temp = 0
        else:
            temp += 1
    
    return sorted_tree_blocks

if __name__ == "__main__":
    nwk = "((A,(B,C):NODE1):NODE2,(D,E):NODE3);"
    tree_blocks = build_tree_block(nwk)
    print(tree_blocks)
    # print(annotate_mrca_position(tree_blocks, ['A', 'B', 'C']))
    # print(extract_subtree(tree_blocks, ['A', 'B', 'C']))
    # print(query_genetic_distance(tree_blocks, [['A', 'B'], ['C', 'D']]))  
    # print(inspect_newick_string(rotate_clade(tree_blocks, 7)))
    print(inspect_newick_string(sort_topology(tree_blocks, ascending=True)))
