import argparse
import csv
from Bio import SeqIO


def get_serotype(conclusion):
    if not conclusion or conclusion == "NA":
        return "Unknown"
    conclusion = str(conclusion).upper()

    if "DENV-1" in conclusion or "DENV1" in conclusion: return "DENV-1"
    if "DENV-2" in conclusion or "DENV2" in conclusion: return "DENV-2"
    if "DENV-3" in conclusion or "DENV3" in conclusion: return "DENV-3"
    if "DENV-4" in conclusion or "DENV4" in conclusion: return "DENV-4"

    if "NC_001477" in conclusion: return "DENV-1"
    if "NC_001474" in conclusion: return "DENV-2"
    if "NC_001475" in conclusion: return "DENV-3"
    if "NC_002640" in conclusion: return "DENV-4"

    return "Unknown"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequences', required=True)
    parser.add_argument('--metadata', required=True)
    args = parser.parse_args()

    sample_serotype = {}
    with open(args.metadata, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample_serotype[row['sample_id']] = get_serotype(row.get('conclusion', 'NA'))

    serotype_seqs = {}
    for record in SeqIO.parse(args.sequences, "fasta"):
        serotype = sample_serotype.get(record.id, "Unknown")
        if serotype == "Unknown":
            continue
        serotype_seqs.setdefault(serotype, []).append(record)

    for serotype, records in serotype_seqs.items():
        if len(records) < 3:
            print(f"Skipping {serotype}: only {len(records)} sequences")
            continue
        filename = f"{serotype.replace('-', '')}.fasta"
        SeqIO.write(records, filename, "fasta")
        print(f"Wrote {len(records)} sequences for {serotype} to {filename}")


if __name__ == "__main__":
    main()