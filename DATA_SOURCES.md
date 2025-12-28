# TCGS-SEQUENTION: Complete Database Sources
## Where to Find Data for Synchronous Parallel Emergence Testing

---

## 1. E. COLI LTEE (LONG-TERM EVOLUTION EXPERIMENT)

### Primary Data Sources

#### A. Good et al. 2017 - Metagenomics (60,000 generations)
**"The dynamics of molecular evolution over 60,000 generations"**
- **GitHub Repository**: https://github.com/benjaminhgood/LTEE-metagenomic
- **Nature Paper**: https://doi.org/10.1038/nature24287
- **Contents**:
  - Allele frequency trajectories for all 12 populations
  - ~60,000 generations of data
  - Time-resolved mutation frequencies
  - Population-level statistics

```bash
# Download command
git clone https://github.com/benjaminhgood/LTEE-metagenomic.git
```

#### B. Tenaillon et al. 2016 - Mutation Catalog (50,000 generations)
**"Tempo and mode of genome evolution in a 50,000-generation experiment"**
- **Nature Paper**: https://doi.org/10.1038/nature18959
- **Supplementary Data**: https://www.nature.com/articles/nature18959#Sec24
- **Contents**:
  - Complete mutation lists for all populations
  - Timing of mutations (generation of detection)
  - Gene annotations
  - Parallel mutation identification

#### C. Barrick & Lenski - breseq Mutation Calls
- **breseq Software**: https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing
- **LTEE Genomes**: https://barricklab.org/twiki/bin/view/Lab/ToolsLTEESequenceData
- **Contents**:
  - Raw sequencing data
  - Processed mutation calls
  - Ancestral reference (REL606)

#### D. Lenski Lab Main Resources
- **LTEE Homepage**: http://myxo.css.msu.edu/ecoli/
- **Fitness Data**: http://myxo.css.msu.edu/ecoli/summdata.html
- **Publications**: http://myxo.css.msu.edu/ecoli/papers.html

### NCBI Data Archives

#### SRA (Sequence Read Archive)
- **BioProject PRJNA295606**: LTEE metagenomics
  - https://www.ncbi.nlm.nih.gov/bioproject/PRJNA295606
  
- **BioProject PRJNA188481**: LTEE clones
  - https://www.ncbi.nlm.nih.gov/bioproject/PRJNA188481

```bash
# Download SRA data using SRA Toolkit
prefetch PRJNA295606
fastq-dump --split-files SRR*
```

#### GenBank Reference Genomes
- **REL606 (Ancestor)**: NC_012967.1
  - https://www.ncbi.nlm.nih.gov/nuccore/NC_012967.1

---

## 2. DROSOPHILA DATASETS

### DGRP (Drosophila Genetic Reference Panel)
- **Main Site**: http://dgrp2.gnets.ncsu.edu/
- **Data Download**: http://dgrp2.gnets.ncsu.edu/data.html
- **Contents**:
  - 205 inbred lines from Raleigh, NC
  - Whole genome sequences
  - Phenotype measurements
  - GWAS results

### DGN (Drosophila Genome Nexus)
- **Main Site**: https://www.johnpool.net/genomes.html
- **FTP**: ftp://ftp.hgsc.bcm.edu/DGN/
- **Contents**:
  - Population samples from multiple continents
  - VCF files with variants
  - Population structure data

### Evolve & Resequence (E&R) Studies
- **Dryad Repository**: https://datadryad.org/
  Search: "Drosophila evolve resequence"
- **Key datasets**:
  - Barghi et al. 2019: Thermal adaptation
  - Graves et al. 2017: Hypoxia adaptation
  - Burke et al. 2010: Accelerated development

---

## 3. OTHER MICROBIAL EVOLUTION DATASETS

### Pseudomonas fluorescens
- **Rainey Lab**: https://www.raineylab.com/resources
- **Key papers**: Adaptive radiation studies

### Saccharomyces cerevisiae (Yeast)
- **1000 Genomes**: http://1002genomes.u-strasbg.fr/
- **SGD**: https://www.yeastgenome.org/
- **E&R Studies**: Search Dryad/NCBI for "yeast experimental evolution"

### LTEE Spinoff Experiments
- **STLE (Short-Term)**: Various publications
- **Citrate+ Evolution**: Blount et al. papers

---

## 4. DOWNLOAD SCRIPTS

### Automated LTEE Data Download

```bash
#!/bin/bash
# download_ltee_data.sh
# Downloads all primary LTEE data sources

mkdir -p ltee_data && cd ltee_data

echo "=== Downloading LTEE Data ==="

# 1. Good et al. 2017 metagenomics
echo "Downloading Good et al. 2017 metagenomics..."
git clone https://github.com/benjaminhgood/LTEE-metagenomic.git
echo "Done: LTEE-metagenomic/"

# 2. Tenaillon et al. 2016 supplementary
echo "Downloading Tenaillon et al. 2016 supplementary..."
# Note: May need manual download from Nature website
wget -O tenaillon2016_supp.zip \
  "https://static-content.springer.com/esm/art%3A10.1038%2Fnature18959/MediaObjects/41586_2016_BFnature18959_MOESM1_ESM.zip" \
  2>/dev/null || echo "Manual download needed from Nature website"

# 3. REL606 reference genome
echo "Downloading REL606 reference genome..."
wget -O REL606.fna.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/985/GCF_000017985.1_ASM1798v1/GCF_000017985.1_ASM1798v1_genomic.fna.gz"
gunzip REL606.fna.gz

# 4. Processed mutation data (if available)
echo "Checking for processed data..."
# Add URLs as they become available

echo ""
echo "=== Download Complete ==="
echo "Contents:"
ls -la

echo ""
echo "Next steps:"
echo "1. Check LTEE-metagenomic/README for data format"
echo "2. Download Tenaillon supplementary from Nature if needed"
echo "3. Run Sequention protocol on allele frequency data"
```

### Python Data Loader for LTEE

```python
#!/usr/bin/env python3
"""
LTEE Data Loader for Sequention Protocol
"""

import pandas as pd
import numpy as np
from pathlib import Path

def load_good2017_allele_frequencies(data_dir: str) -> dict:
    """
    Load allele frequency data from Good et al. 2017.
    
    Parameters
    ----------
    data_dir : str
        Path to LTEE-metagenomic directory
    
    Returns
    -------
    dict
        Population ID -> DataFrame with columns [generation, position, frequency]
    """
    data_path = Path(data_dir)
    
    populations = {}
    
    # Population directories
    pop_ids = ['Ara-1', 'Ara-2', 'Ara-3', 'Ara-4', 'Ara-5', 'Ara-6',
               'Ara+1', 'Ara+2', 'Ara+3', 'Ara+4', 'Ara+5', 'Ara+6']
    
    for pop_id in pop_ids:
        # Look for timecourse files
        freq_file = data_path / f"{pop_id}_timecourse.txt"
        if freq_file.exists():
            df = pd.read_csv(freq_file, sep='\t')
            populations[pop_id] = df
            print(f"Loaded {pop_id}: {len(df)} mutations")
    
    return populations


def load_tenaillon2016_mutations(excel_file: str) -> pd.DataFrame:
    """
    Load mutation timing data from Tenaillon et al. 2016 supplementary.
    
    Parameters
    ----------
    excel_file : str
        Path to supplementary Excel file
    
    Returns
    -------
    pd.DataFrame
        Mutation catalog with timing information
    """
    df = pd.read_excel(excel_file, sheet_name='Mutations')
    return df


def extract_gene_timing(mutations_df: pd.DataFrame, gene: str) -> dict:
    """
    Extract timing of mutations in a specific gene across populations.
    
    Parameters
    ----------
    mutations_df : pd.DataFrame
        Full mutation catalog
    gene : str
        Gene name (e.g., 'spoT', 'topA')
    
    Returns
    -------
    dict
        Population ID -> generation of first detection
    """
    gene_muts = mutations_df[mutations_df['gene'] == gene]
    
    timing = {}
    for pop_id in gene_muts['population'].unique():
        pop_muts = gene_muts[gene_muts['population'] == pop_id]
        first_gen = pop_muts['generation'].min()
        timing[pop_id] = first_gen
    
    return timing


if __name__ == "__main__":
    # Example usage
    print("LTEE Data Loader")
    print("================")
    print()
    print("To use:")
    print("1. Clone: git clone https://github.com/benjaminhgood/LTEE-metagenomic.git")
    print("2. Run: populations = load_good2017_allele_frequencies('LTEE-metagenomic/')")
    print("3. Analyze with Sequention protocol")
```

---

## 5. QUICK START GUIDE

### Step 1: Download LTEE Metagenomics
```bash
git clone https://github.com/benjaminhgood/LTEE-metagenomic.git
cd LTEE-metagenomic
ls -la  # Examine structure
```

### Step 2: Examine Data Format
```bash
# The repository contains:
# - Allele frequency timecourses
# - Processed SFS data
# - Analysis scripts from the paper
head -20 README.md
```

### Step 3: Run Sequention Protocol
```python
from sequention import ProtocolConfig, SequentionProtocol

config = ProtocolConfig()
config.synchrony.delta_t = 2000  # 2000 generations window for LTEE

protocol = SequentionProtocol(config)
protocol.load_data(
    generic_file="LTEE-metagenomic/data/allele_freqs.csv",
    generic_format={
        'population_col': 'population',
        'time_col': 'generation',
        'value_cols': ['frequency']
    }
)

result = protocol.run()
print(f"Decision: {result.decision}")
```

---

## 6. KEY CONTACTS

### Lenski Lab (Michigan State University)
- **Richard Lenski**: lenski@msu.edu
- **Lab Website**: http://myxo.css.msu.edu/
- For: LTEE frozen stocks, unpublished data, collaboration

### Desai Lab (Harvard University)  
- **Michael Desai**: mdesai@oeb.harvard.edu
- For: Metagenomics analysis, population genetics methods

### Good (UC Berkeley)
- **Benjamin Good**: bgood@berkeley.edu
- For: LTEE-metagenomic data questions, analysis methods

---

## 7. DATA SUMMARY TABLE

| Dataset | Organisms | Generations | Populations | URL |
|---------|-----------|-------------|-------------|-----|
| LTEE Metagenomics | E. coli | 60,000 | 12 | github.com/benjaminhgood/LTEE-metagenomic |
| LTEE Mutations | E. coli | 50,000 | 12 | Nature Supp (Tenaillon 2016) |
| DGRP | Drosophila | N/A | 205 lines | dgrp2.gnets.ncsu.edu |
| DGN | Drosophila | N/A | ~200 | johnpool.net/genomes.html |
| 1002 Yeast | S. cerevisiae | N/A | 1,011 | 1002genomes.u-strasbg.fr |

---

## 8. CITATION REQUIREMENTS

When using these data, cite:

**LTEE Metagenomics:**
> Good BH, McDonald MJ, Barrick JE, Lenski RE, Desai MM. The dynamics of molecular evolution over 60,000 generations. Nature. 2017;551:45-50.

**LTEE Mutations:**
> Tenaillon O, Barrick JE, Ribeck N, et al. Tempo and mode of genome evolution in a 50,000-generation experiment. Nature. 2016;536:165-170.

**DGRP:**
> Mackay TF, et al. The Drosophila melanogaster Genetic Reference Panel. Nature. 2012;482:173-178.

---

*Document prepared: December 2025*
*For TCGS-SEQUENTION Protocol v1.0.0*
