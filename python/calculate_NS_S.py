import gzip
import json
import re
from typing import Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.codonalign.codonseq import *
from Bio.codonalign import codonseq
from Bio.Data import CodonTable


def calculate_ns_s(sequence: Seq) -> Tuple[int, int]:

    # sanity checks
    assert len(sequence) % 3 == 0, "The sequence must be a multiple of 3"
    for b in sequence:
        assert b in ["A", "C", "G", "T"], "The sequence contains other bases than A, C, G, T: {}".format(b)

    # translate the sequence
    protein = sequence.translate()
    ns = 0
    s = 0

    # introduce all possible point mutations
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        sequence_start = sequence[0:i]
        sequence_end = sequence[i + 3:]
        for j in range(3):
            base = codon[j]
            for k in {"A", "C", "G", "T"}.difference(base):
                # mutation k in codon position j in ith codon within sequence
                new_codon = codon[0:j] + k + codon[j+1:]
                new_sequence = sequence_start + new_codon + sequence_end
                new_protein = new_sequence.translate()
                if new_protein == protein:
                    # counts a synonymous mutation if the translated sequence does not change
                    s += 1
                else:
                    # counts a non synonymous mutation if the translated sequence changes
                    ns += 1

    return ns, s


def calculate_ns_s_2(sequence: Seq) -> Tuple[int, int]:
    codon_sequence = CodonSeq(str(sequence))
    s, ns = codonseq._count_site_NG86(
        codon_lst=codonseq._get_codon_list(codon_sequence)[0:-1],
        codon_table=CodonTable.standard_dna_table
    )
    return ns, s


# calculates NS and S over genes in SARS-CoV-2
def calculate_ns_s_over_genes():
    data = {"transcript": [], "gene": [], "NS": [], "S": []}
    with gzip.open("../references/Sars_cov_2.ASM985889v3.cds.all.fa.gz", "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            transcript = re.sub("\\..*", "", record.id)
            gene = re.compile(r'.*gene_symbol:(\w+)').match(record.description).group(1)
            ns, s = calculate_ns_s_2(sequence=record.seq)
            data.get("transcript").append(transcript)
            data.get("gene").append(gene)
            data.get("NS").append(ns)
            data.get("S").append(s)

    df = pd.DataFrame(data=data)
    df.to_csv("../references/genes_NS_S_2.csv", index=False)


calculate_ns_s_over_genes()


# calculates NS and S over domains in SARS-CoV-2
def calculate_ns_s_over_domains():

    gene_sequences = {}
    with gzip.open("../references/Sars_cov_2.ASM985889v3.cds.all.fa.gz", "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            gene = re.compile(r'.*gene_symbol:(\w+)').match(record.description).group(1)
            gene_sequences[re.sub("\\..*", "", record.id)] = (gene, record.seq)

    print(gene_sequences)

    data = json.load(open("../references/sars_cov_2.json"))
    domains = {"transcript": [], "gene": [], "domain": [], "start": [], "end": [], "NS": [], 'S': []}
    for g in data.get('genes'):
        transcript = g.get("transcripts", [])[0]
        transcript_id = transcript.get('id')
        transcript_sequence = gene_sequences.get(transcript_id)[1]
        protein_features = transcript.get("translations", [])[0].get("protein_features")
        pfam_protein_features = [f for f in protein_features if f.get("dbname") == "Pfam"]
        gene = g.get("name")
        for d in pfam_protein_features:
            domains.get("transcript").append(transcript_id)
            domains.get("gene").append(gene)
            domains.get("domain").append(d.get('description'))
            # convert here 1-based protein coordinate to 0-based DNA coordinates
            start = (int(d.get('start')) - 1) * 3
            end = (int(d.get('end')) - 1) * 3
            ns, s = calculate_ns_s_2(transcript_sequence[start:end])
            domains.get("NS").append(ns)
            domains.get("S").append(s)
            domains.get("start").append(start)
            domains.get("end").append(end)
    df = pd.DataFrame(data=domains)
    df.to_csv("../references/domains_NS_S_2.csv", index=False)


calculate_ns_s_over_domains()
