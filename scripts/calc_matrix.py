import argparse
from Bio import AlignIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help="Aligned FASTA file")
    parser.add_argument('--output', required=True, help="Output Matrix TSV")
    args = parser.parse_args()

    alignment = AlignIO.read(args.input, "fasta")
    num_seqs = len(alignment)
    ids = [rec.id for rec in alignment]
    seqs = [str(rec.seq).upper() for rec in alignment]
    
    matrix = [[0] * num_seqs for _ in range(num_seqs)]
    
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            s1 = seqs[i]
            s2 = seqs[j]
            dist = 0
            length = len(s1)
            
            for k in range(length):
                c1 = s1[k]
                c2 = s2[k]
                
                if c1 in '-N' or c2 in '-N':
                    continue
                
                if c1 != c2:
                    dist += 1
            
            matrix[i][j] = dist
            matrix[j][i] = dist

    with open(args.output, "w") as f:
        f.write("snp-dists\t" + "\t".join(ids) + "\n")
        for i, row in enumerate(matrix):
            f.write(ids[i] + "\t" + "\t".join(map(str, row)) + "\n")

if __name__ == "__main__":
    main()