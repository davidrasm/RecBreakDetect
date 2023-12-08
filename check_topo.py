# -*- coding: utf-8 -*-
"""
-------------------------------------------------
   File Name：     check_topo
   Description :
   Author :       shicen
   date：          7/11/22
-------------------------------------------------
   Change Activity:
                   7/11/22:
-------------------------------------------------
"""


# def print_topology(t_newick):
#     t1_t = []
#     pairs = []
#     node_list = [""]
#     optr = []
#     t1 = t_newick.strip(';')
#     first_node = []
#
#     i = 0
#     while i < len(t1):
#         op = t1[i]
#         if op not in ['(', ')', ',']:
#             node_list[-1] += op
#             i += 1
#         # print(t1[i])
#         elif op == '(':
#             optr.append(op)
#             i += 1
#             continue
#         elif op == ')':
#             if optr[-1] == ',':
#                 node1 = node_list.pop()
#                 node2 = node_list.pop()
#                 if not first_node:
#                     first_node.append(node1)
#                     first_node.append(node2)
#                 pairs.append(sorted({node1, node2}))
#                 pairs = sorted(pairs)
#                 node_list.append("".join(sorted([node1, node2])))
#                 optr.pop()
#                 continue
#             elif optr[-1] == '(':
#                 optr.pop()
#                 i += 1
#                 continue
#             else:
#                 break
#         elif op == ',':
#             node_list.append("")
#             optr.append(op)
#             i += 1
#             continue
#
#     # print(pairs)
#     return pairs, first_node
import math
import re


# @staticmethod
def delete_branch_length(tree_newick):
    branch_length_pattern = re.compile(r':\d+\.\d*')
    new_tree = re.sub(branch_length_pattern, "", tree_newick)

    return new_tree


def print_topology(t_newick):

    t1 = delete_branch_length(t_newick)
    t1 = t1.replace("n", "").strip().strip(';')
    # print(t1)
    pairs = []
    node_list = []
    optr = []
    first_node = []

    i = 0
    while i < len(t1):
        op = t1[i]
        if op not in ['(', ')', ',']:
            node_list.append(2 ** int(op))
            i += 1
        # print(t1[i])
        elif op == '(':
            optr.append(op)
            i += 1
            continue
        elif op == ')':
            if optr[-1] == ',':
                node1 = node_list.pop()
                node2 = node_list.pop()
                if not first_node:
                    first_node.append(math.log2(node1))
                    first_node.append(math.log2(node2))
                pairs.append(node1 + node2)
                node_list.append(node1 + node2)
                optr.pop()
                continue
            elif optr[-1] == '(':
                optr.pop()
                i += 1
                continue
            else:
                break
        elif op == ',':
            # node_list.append("")
            optr.append(op)
            i += 1
            continue

    # print(pairs)
    return pairs, first_node


if __name__ == "__main__":

    t1 = '((n0,n2),((n1,n3),n4));'
    t1 = '(n1,((n0),(n2)));'
    t2 = '(n1,((n2,(n0))));'

    topo_1, _ = print_topology(t1)
    topo_2, _ = print_topology(t2)

    print(topo_1)
    print(topo_2)

    print(topo_1[0] == topo_2[0])
