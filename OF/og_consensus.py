#!/usr/bin/env python3
import os
import re
import glob
import argparse
from collections import Counter, defaultdict

import pandas as pd

# ---------- Helpers ----------

def read_emapper_annotations(path):
    """
    Robust reader: finds the last header line starting with '#'
    and uses it as column names; then reads data rows.
    """
    header = None
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('#'):
                if line.lower().startswith('#query') or line.lower().startswith('#query_name'):
                    header = line[1:].strip().split('\t')
    if header is None:
        # Fallback to a common set of columns; your version may differ.
        header = [
            'query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_level',
            'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC',
            'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'BRITE',
            'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'
        ]
    df = pd.read_csv(path, sep='\t', comment='#', header=None, names=header, dtype=str)
    # Coerce score/evalue if present
    if 'score' in df.columns:
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
    if 'evalue' in df.columns:
        # Smaller is better
        df['evalue'] = pd.to_numeric(df['evalue'], errors='coerce')
    return df

def split_list_field(val, sep=','):
    """
    Split a list-like field into normalized tokens.
    Handles empty/NaN safely.
    """
    if pd.isna(val) or val == '':
        return []
    # Some fields may use commas; some contain decorations like PF12345(1e-10)
    parts = [p.strip() for p in str(val).replace(';', ',').split(sep)]
    parts = [p for p in parts if p]
    return parts

def normalize_pfam_token(tok):
    # PFAMs often look like PF00001.20(1.2e-05)
    # Keep the PFxxxxx accession only
    m = re.match(r'(PF\d{5}(?:\.\d+)?)', tok)
    return m.group(1) if m else tok

def freq_terms(series, normalize_func=lambda x: x, threshold=0.3, topn=5):
    """
    For a series of list-like strings, compute term frequencies and
    return a list of "term|support_fraction" strings for terms meeting threshold.
    """
    all_terms = []
    for v in series.fillna(''):
        for t in split_list_field(v):
            tnorm = normalize_func(t)
            if tnorm:
                all_terms.append(tnorm)
    if not all_terms:
        return []
    total = series.shape[0]
    cnt = Counter(all_terms)
    ranked = [(term, n/total, n) for term, n in cnt.most_common()]
    selected = [(t, f) for (t, f, n) in ranked if f >= threshold]
    if topn is not None:
        selected = selected[:topn]
    return [f"{t}|{f:.2f}" for t, f in selected]

def majority_text(series, tie_break_by=None):
    """
    Pick most frequent non-empty text. If tie, choose the row with best 'score' (max),
    otherwise shortest/cleanest string.
    """
    values = [s for s in series.fillna('') if s.strip()]
    if not values:
        return ''
    cnt = Counter(values)
    top_freq = cnt.most_common()
    best = [v for v, n in top_freq if n == top_freq[0][1]]
    if len(best) == 1 or tie_break_by is None:
        # pick the "cleanest" shortest string to avoid overly verbose descriptions
        return sorted(best, key=lambda s: (len(s), s))[0]
    # tie-break using score: pick the description of the row with max score
    idx_max = tie_break_by['score'].idxmax() if 'score' in tie_break_by else None
    if idx_max is not None:
        val = tie_break_by.loc[idx_max, series.name]
        return val if isinstance(val, str) else best[0]
    return sorted(best, key=lambda s: (len(s), s))[0]

def summarize_orthogroup(og_id, df, threshold=0.3, topn=5):
    """
    Build a dict of consensus annotation for one orthogroup.
    """
    n_seqs = df.shape[0]
    # Count non-empty per field
    def nonempty(col): 
        return (df[col].astype(str).str.strip() != '').sum() if col in df.columns else 0

    consensus = {
        'orthogroup': og_id,
        'n_sequences': n_seqs,
        'n_with_description': nonempty('Description'),
        'n_with_preferred_name': nonempty('Preferred_name'),
        'n_with_GO': nonempty('GOs'),
        'n_with_KO': nonempty('KEGG_ko'),
        'n_with_EC': nonempty('EC'),
        'n_with_PFAM': nonempty('PFAMs'),
        'COG_category_consensus': '',
        'Description_consensus': '',
        'Preferred_name_consensus': '',
        'GO_BP_consensus': '',
        'GO_MF_consensus': '',
        'GO_CC_consensus': '',
        'KO_consensus': '',
        'KEGG_Pathway_consensus': '',
        'KEGG_Module_consensus': '',
        'EC_consensus': '',
        'PFAM_core_domains': '',
    }

    # Text consensus
    if 'Description' in df.columns:
        consensus['Description_consensus'] = majority_text(df['Description'], tie_break_by=df)
    if 'Preferred_name' in df.columns:
        consensus['Preferred_name_consensus'] = majority_text(df['Preferred_name'], tie_break_by=df)

    # COG categories (single-letter codes, can be multiple per gene)
    if 'COG_category' in df.columns:
        letters = []
        for v in df['COG_category'].fillna(''):
            letters += [c for c in v if c.isalpha()]
        if letters:
            cnt = Counter(letters)
            top = [f"{k}|{cnt[k]/n_seqs:.2f}" for k in [x for x,_ in cnt.most_common()]]
            # report those ≥ threshold, or at least top 1 if none meet threshold
            top_keep = [t for t in top if float(t.split('|')[1]) >= threshold]
            if not top_keep and top:
                top_keep = top[:1]
            consensus['COG_category_consensus'] = ','.join(top_keep)

    # GO terms → split by namespace (BP/MF/CC) if provided in IDs (GO:xxxx; emapper does not tag namespace)
    # Here, we report overall top terms; optionally you can separate namespaces later via GOATOOLS.
    if 'GOs' in df.columns:
        go_terms = freq_terms(df['GOs'], normalize_func=lambda x: x, threshold=threshold, topn=topn)
        # As emapper doesn't include namespace, keep a single list; or split later if you have a GO mapping.
        consensus['GO_BP_consensus'] = ','.join(go_terms)  # keep in one field; rename later if needed
        consensus['GO_MF_consensus'] = ''
        consensus['GO_CC_consensus'] = ''

    # KEGG
    if 'KEGG_ko' in df.columns:
        ko = freq_terms(df['KEGG_ko'], threshold=threshold, topn=topn)
        consensus['KO_consensus'] = ','.join(ko)
    if 'KEGG_Pathway' in df.columns:
        kp = freq_terms(df['KEGG_Pathway'], threshold=threshold, topn=topn)
        consensus['KEGG_Pathway_consensus'] = ','.join(kp)
    if 'KEGG_Module' in df.columns:
        km = freq_terms(df['KEGG_Module'], threshold=threshold, topn=topn)
        consensus['KEGG_Module_consensus'] = ','.join(km)

    # EC numbers
    if 'EC' in df.columns:
        ec = freq_terms(df['EC'], threshold=threshold, topn=topn)
        consensus['EC_consensus'] = ','.join(ec)

    # PFAM “core” domains (normalize PFxxxxx accessions)
    if 'PFAMs' in df.columns:
        pf = freq_terms(df['PFAMs'], normalize_func=normalize_pfam_token, threshold=threshold, topn=topn)
        consensus['PFAM_core_domains'] = ','.join(pf)

    return consensus

# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser(description="Summarize EggNOG-mapper annotations into per-orthogroup consensus.")
    ap.add_argument('--ann_dir', required=True, help="Directory with *.emapper.annotations files.")
    ap.add_argument('--out', default='orthogroup_consensus.tsv', help="Output TSV for consensus.")
    ap.add_argument('--members_out', default='orthogroup_membership.tsv', help="Optional long table of per-gene annotations.")
    ap.add_argument('--threshold', type=float, default=0.3, help="Support fraction threshold (e.g., 0.3 means ≥30% of members).")
    ap.add_argument('--topn', type=int, default=5, help="Max terms to report per field.")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.ann_dir, '*.emapper.annotations')))
    if not files:
        raise SystemExit(f"No .emapper.annotations found in {args.ann_dir}")

    consensus_rows = []
    members = []

    for fp in files:
        og_id = os.path.basename(fp).replace('.emapper.annotations', '')
        df = read_emapper_annotations(fp)
        df['orthogroup'] = og_id
        members.append(df)

        cons = summarize_orthogroup(og_id, df, threshold=args.threshold, topn=args.topn)
        consensus_rows.append(cons)

    cons_df = pd.DataFrame(consensus_rows).sort_values('orthogroup')
    cons_df.to_csv(args.out, sep='\t', index=False)

    # Also write a tidy per-gene table (useful for QA/QC and custom plots)
    members_df = pd.concat(members, ignore_index=True)
    # Keep a subset of informative columns if present
    keep_cols = [c for c in [
        'orthogroup','query','score','evalue','COG_category','Description','Preferred_name',
        'GOs','EC','KEGG_ko','KEGG_Pathway','KEGG_Module','PFAMs'
    ] if c in members_df.columns]
    members_df[keep_cols].to_csv(args.members_out, sep='\t', index=False)

    print(f"✅ Wrote {args.out} (consensus) and {args.members_out} (per-gene details).")

if __name__ == '__main__':
    main()
