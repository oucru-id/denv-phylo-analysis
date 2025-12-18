import json
import argparse
import os
import csv
import re  
from Bio import SeqIO

def load_references(ref_paths):
    refs = {}
    for path in ref_paths:
        ref_id = os.path.basename(path).replace('.fasta', '')
        record = SeqIO.read(path, "fasta")
        refs[ref_id] = str(record.seq).upper()
    return refs

def parse_fhir(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    sample_id = os.path.basename(file_path).replace('.fhir.json', '').replace('.merged', '')
    consensus_seq = None
    
    metadata = {
        "sample_id": sample_id,
        "patient_id": "NA",
        "latitude": "NA",
        "longitude": "NA",
        "conclusion": "NA"
    }
    
    if 'entry' in data:
        for entry in data['entry']:
            res = entry.get('resource', {})
            r_type = res.get('resourceType')

            if r_type == 'Patient':
                metadata["patient_id"] = res.get('id', 'NA')

            elif r_type == 'DiagnosticReport':
                conclusions = []
                if res.get('conclusion'):
                    conclusions.append(res.get('conclusion'))
                if conclusions:
                    metadata["conclusion"] = "; ".join(conclusions)

            elif r_type == 'Observation':
                is_consensus = False
                coding = res.get('code', {}).get('coding', [])
                if any(c.get('code') == '86206-0' for c in coding):
                    is_consensus = True
                elif res.get('code', {}).get('text') == 'Viral Consensus Genome Sequence':
                    is_consensus = True
                
                if is_consensus and 'valueString' in res:
                    raw_seq = res['valueString']
                    consensus_seq = re.sub(r'[^ATCGN]', '', raw_seq.upper())

    return sample_id, consensus_seq, metadata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True, help="List of FHIR JSON files")
    parser.add_argument('--references', nargs='+', required=True, help="List of Reference FASTA files")
    args = parser.parse_args()

    refs = load_references(args.references)
    
    sample_seqs = {}
    all_metadata = []

    for f in args.inputs:
        sid, seq, meta = parse_fhir(f)
        if seq:
            sample_seqs[sid] = seq
            all_metadata.append(meta)
        else:
            print(f"No consensus sequence found for {sid}")

    with open("metadata.tsv", "w", newline='') as f:
        fieldnames = ["sample_id", "patient_id", "latitude", "longitude", "conclusion"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(all_metadata)
        for rid in refs:
            writer.writerow({
                "sample_id": rid,
                "patient_id": "Reference",
                "latitude": "NA",
                "longitude": "NA",
                "conclusion": f"Reference {rid}"
            })

    with open("unaligned_sequences.fasta", "w") as f:
        for rid, seq in refs.items():
            f.write(f">{rid}\n{seq}\n")
        for sid, seq in sample_seqs.items():
            f.write(f">{sid}\n{seq}\n")

if __name__ == "__main__":
    main()