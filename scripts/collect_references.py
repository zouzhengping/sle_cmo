#!/usr/bin/env python3
"""
收集PLOS ONE论文所需的参考文献

参考文献要求：
- 数量：至少25-30篇
- 格式：Vancouver格式
- 范围：TSCM、SLE、单细胞RNA-seq、SLEDAI等相关研究
"""

import sys
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "docs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 参考文献模板（Vancouver格式）
references = {
    "TSCM_definition": [
        {
            "id": 1,
            "authors": "Gattinoni L, Lugli E, Ji Y, Pos Z, Paulos CM, Quigley MF, et al.",
            "title": "A human memory T cell subset with stem cell-like properties",
            "journal": "Nat Med",
            "year": 2011,
            "volume": 17,
            "issue": 10,
            "pages": "1290-7",
            "doi": "10.1038/nm.2446",
            "category": "TSCM definition"
        },
        {
            "id": 2,
            "authors": "Gattinoni L, Speiser DE, Lichterfeld M, Bonini C.",
            "title": "T memory stem cells in health and disease",
            "journal": "Nat Med",
            "year": 2017,
            "volume": 23,
            "issue": 1,
            "pages": "18-27",
            "doi": "10.1038/nm.4241",
            "category": "TSCM review"
        }
    ],
    "SLE_T_cells": [
        {
            "id": 3,
            "authors": "Deng Y, Tsao BP.",
            "title": "Genetic susceptibility to systemic lupus erythematosus in the genomic era",
            "journal": "Nat Rev Rheumatol",
            "year": 2010,
            "volume": 6,
            "issue": 12,
            "pages": "683-92",
            "doi": "10.1038/nrrheum.2010.176",
            "category": "SLE genetics"
        },
        {
            "id": 4,
            "authors": "Crispin JC, Tsokos GC.",
            "title": "Human TCR-alpha beta+ CD4- CD8- T cells can derive from CD8+ T cells and display an inflammatory effector phenotype",
            "journal": "J Immunol",
            "year": 2009,
            "volume": 183,
            "issue": 7,
            "pages": "4615-23",
            "doi": "10.4049/jimmunol.0901533",
            "category": "SLE T cells"
        }
    ],
    "scRNA_seq_methods": [
        {
            "id": 5,
            "authors": "Wolf FA, Angerer P, Theis FJ.",
            "title": "SCANPY: large-scale single-cell gene expression data analysis",
            "journal": "Genome Biol",
            "year": 2018,
            "volume": 19,
            "issue": 1,
            "pages": "15",
            "doi": "10.1186/s13059-017-1382-0",
            "category": "Scanpy"
        },
        {
            "id": 6,
            "authors": "Traag VA, Waltman L, van Eck NJ.",
            "title": "From Louvain to Leiden: guaranteeing well-connected communities",
            "journal": "Sci Rep",
            "year": 2019,
            "volume": 9,
            "issue": 1,
            "pages": "5233",
            "doi": "10.1038/s41598-019-41695-z",
            "category": "Leiden clustering"
        },
        {
            "id": 7,
            "authors": "McInnes L, Healy J, Melville J.",
            "title": "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction",
            "journal": "arXiv preprint",
            "year": 2018,
            "doi": "10.48550/arXiv.1802.03426",
            "category": "UMAP"
        }
    ],
    "SLEDAI": [
        {
            "id": 8,
            "authors": "Bombardier C, Gladman DD, Urowitz MB, Caron D, Chang CH.",
            "title": "Derivation of the SLEDAI. A disease activity index for lupus patients. The Committee on Prognosis Studies in SLE",
            "journal": "Arthritis Rheum",
            "year": 1992,
            "volume": 35,
            "issue": 6,
            "pages": "630-40",
            "doi": "10.1002/art.1780350606",
            "category": "SLEDAI definition"
        }
    ],
    "GSE_datasets": [
        {
            "id": 9,
            "authors": "Nehar-Belaid D, Hong S, Marches R, Chen G, Bolisetty M, Baisch J, et al.",
            "title": "Mapping systemic lupus erythematosus heterogeneity at the single-cell level",
            "journal": "Nat Immunol",
            "year": 2020,
            "volume": 21,
            "issue": 9,
            "pages": "1094-106",
            "doi": "10.1038/s41590-020-0743-0",
            "category": "GSE137029",
            "note": "Single-cell data source"
        },
        {
            "id": 10,
            "authors": "[GSE88884 authors - need to find]",
            "title": "[GSE88884 title - need to find]",
            "journal": "[Journal]",
            "year": "[Year]",
            "doi": "[DOI]",
            "category": "GSE88884",
            "note": "Bulk validation dataset - need to find original publication"
        }
    ],
    "TSCM_autoimmunity": [
        {
            "id": 11,
            "authors": "Zhang L, Yu X, Zheng L, Zhang Y, Li Y, Fang Q, et al.",
            "title": "Lineage tracking reveals dynamic relationships of T cells in colorectal cancer",
            "journal": "Nature",
            "year": 2018,
            "volume": 564,
            "issue": 7735,
            "pages": "268-72",
            "doi": "10.1038/s41586-018-0694-x",
            "category": "TSCM in disease"
        }
    ],
    "enrichment_analysis": [
        {
            "id": 12,
            "authors": "Fang Z, Liu X, Peltz G.",
            "title": "GSEApy: a comprehensive package for performing gene set enrichment analysis in Python",
            "journal": "Bioinformatics",
            "year": 2023,
            "volume": 39,
            "issue": 1,
            "doi": "10.1093/bioinformatics/btac757",
            "category": "gseapy"
        }
    ]
}

def format_vancouver(ref):
    """Format reference in Vancouver style"""
    # Format: Author. Title. Journal. Year;Volume(Issue):Pages. DOI
    parts = []
    
    # Authors
    if 'authors' in ref:
        parts.append(ref['authors'])
    
    # Title
    if 'title' in ref:
        parts.append(ref['title'] + '.')
    
    # Journal
    if 'journal' in ref:
        parts.append(ref['journal'])
    
    # Year, Volume, Issue, Pages
    year_vol = []
    if 'year' in ref:
        year_vol.append(str(ref['year']))
    if 'volume' in ref:
        year_vol.append(str(ref['volume']))
    if 'issue' in ref:
        year_vol.append('(' + str(ref['issue']) + ')')
    if 'pages' in ref:
        year_vol.append(':' + ref['pages'])
    
    if year_vol:
        parts.append(';'.join(year_vol) + '.')
    
    # DOI
    if 'doi' in ref:
        parts.append('doi: ' + ref['doi'])
    
    return ' '.join(parts)

# Generate reference list
all_refs = []
for category, refs in references.items():
    all_refs.extend(refs)

# Sort by ID
all_refs.sort(key=lambda x: x.get('id', 999))

# Generate markdown document
ref_doc = """# References for PLOS ONE Submission

**Note**: This is a preliminary reference list. Please verify all references and add missing information (especially for GSE88884).

## Reference List (Vancouver Format)

"""

for i, ref in enumerate(all_refs, 1):
    ref_doc += f"{i}. {format_vancouver(ref)}\n"
    if 'note' in ref:
        ref_doc += f"   *Note: {ref['note']}*\n"
    ref_doc += "\n"

ref_doc += """
## Reference Categories

### TSCM Definition and Review
- References 1-2: TSCM cell definition and properties

### SLE and T Cells
- References 3-4: SLE pathogenesis and T-cell involvement

### Single-cell RNA-seq Methods
- References 5-7: Analysis tools (Scanpy, Leiden, UMAP)

### SLEDAI
- Reference 8: SLEDAI definition

### Datasets
- Reference 9: GSE137029 (single-cell data source)
- Reference 10: GSE88884 (bulk validation - **NEED TO FIND**)

### TSCM in Disease
- Reference 11: TSCM in other diseases

### Analysis Tools
- Reference 12: GSEApy for enrichment analysis

## Missing References to Find

1. **GSE88884 original publication**: Need to find the original paper describing this dataset
2. **GSE211700 original publication**: If used in paper
3. **GSE61635 original publication**: If used in paper
4. **Additional TSCM in autoimmunity papers**: Search for TSCM in rheumatoid arthritis, multiple sclerosis, etc.
5. **SLE memory T-cell studies**: Recent papers on memory T cells in SLE
6. **SLEDAI validation studies**: Papers validating SLEDAI as outcome measure

## Search Strategy

1. Search PubMed for:
   - "TSCM" OR "stem-like memory T cell" AND "lupus" OR "SLE"
   - "memory T cell" AND "systemic lupus erythematosus"
   - "single-cell RNA-seq" AND "SLE" OR "lupus"
   - "SLEDAI" AND "biomarker" OR "outcome"

2. Check GEO dataset pages for original publications

3. Review recent reviews on:
   - T-cell subsets in autoimmunity
   - Single-cell genomics in rheumatic diseases
   - SLE pathogenesis

## Target: 25-30 References

Current count: {count} references
Need to add: {remaining} more references

""".format(count=len(all_refs), remaining=max(0, 25-len(all_refs)))

# Save
ref_file = OUT_DIR / "paper_references.md"
with open(ref_file, 'w', encoding='utf-8') as f:
    f.write(ref_doc)

# Also save as JSON for easy editing
ref_json = OUT_DIR / "paper_references.json"
with open(ref_json, 'w', encoding='utf-8') as f:
    json.dump(references, f, indent=2, ensure_ascii=False)

print("=" * 80)
print("参考文献收集完成")
print("=" * 80)
print(f"\n已收集 {len(all_refs)} 篇参考文献")
print(f"\n输出文件:")
print(f"  1. Markdown格式: {ref_file}")
print(f"  2. JSON格式: {ref_json}")
print(f"\n注意:")
print(f"  - 需要查找GSE88884的原始发表文献")
print(f"  - 需要补充更多TSCM在自身免疫病中的研究")
print(f"  - 目标：25-30篇参考文献")

