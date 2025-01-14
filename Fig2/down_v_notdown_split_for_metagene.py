#take sig_down and not_down genes from Sam's file (tss_regions_chr6_ring1b_ezh2_peak_overlap_expr_BY_geneID_11_1_24.csv) and get the TSS coordinates from gencode.vM25.comprehensive.tss.
#for each genotype/treatment output a down file and a not_down file. Each file has 2 columns: geneid and tss coordinate
import sys

input_genes = open(sys.argv[1])
input_tss = open(sys.argv[2])


matrix_genes = []
matrix_tss = []

down_list = []
not_down_list = []

def split_down_not_down(col): #split into two lists down or not_down
    for item in matrix_genes:
        gene = item.split(",")[1]
        if gene != 'original_geneID':
            up_down = item.split(",")[col]
            if up_down == 'not_down':
                not_down_list.append(gene)
            elif up_down == 'sig_down':
                down_list.append(gene)
            elif up_down == 'NA':
                continue
            else:
                print(f"You haven't accounted for {up_down}")
def assign_tss(filename1, filename2):
    o_down = open(filename1, 'w')
    o_not_down = open(filename2, 'w')
    for item in matrix_tss:
        geneid = item.split("\t")[0]
        tss_coord = item.split("\t")[1]
        if geneid in down_list:
            o_down.write(f"{geneid}\t{tss_coord}\n")
        if geneid in not_down_list:
            o_not_down.write(f"{geneid}\t{tss_coord}\n")
    o_down.close()
    o_not_down.close()
    down_list.clear()
    not_down_list.clear()
'''
def write_out_files(filename1, filename2):
    o_down = open(filename1, 'w')
    o_not_down = open(filename2, 'w')
    for item in geneid_tss_down:
        o_down.write(f"{item}\t{geneid_tss_down[item]}\n")
    for item in geneid_tss_not_down:
        o_not_down.write(f"{item}\t{geneid_tss_not_down[item]}\n")
    o_down.close()
    o_not_down.close()
'''


for line in input_genes:
    line = line.strip()
    matrix_genes.append(line)

for line in input_tss:
    line = line.strip()
    matrix_tss.append(line)

split_down_not_down(8)
assign_tss('x1000ng_down.txt', 'x1000ng_not_down.txt')

split_down_not_down(10)
assign_tss('x10ng_down.txt', 'x10ng_not_down.txt')

split_down_not_down(12)
assign_tss('x872a_down.txt', 'x872a_not_down.txt')

input_genes.close()
input_tss.close()
