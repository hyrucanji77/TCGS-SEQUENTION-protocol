# Data Sources

## LTEE Gene Timing Data

The mutation timing data in `ltee_genes.json` is compiled from:

### Primary Sources

1. **Good et al. 2017** - *Nature* 551:45-50
   - "The dynamics of molecular evolution over 60,000 generations"
   - DOI: [10.1038/nature24287](https://doi.org/10.1038/nature24287)
   - Metagenomics data repository: https://github.com/benjaminhgood/LTEE-metagenomic

2. **Tenaillon et al. 2016** - *Nature* 536:165-170
   - "Tempo and mode of genome evolution in a 50,000-generation experiment"
   - DOI: [10.1038/nature18959](https://doi.org/10.1038/nature18959)
   - Mutation data in Supplementary Tables

3. **Maddamsetti et al. 2020** - *Genome Biology and Evolution* 12:1591-1603
   - "Divergent evolution of mutation rates and biases in the Long-Term Evolution Experiment"
   - DOI: [10.1093/gbe/evaa178](https://doi.org/10.1093/gbe/evaa178)
   - Mutator status analysis

### LTEE Background

The Long-Term Evolution Experiment (LTEE) was started by Richard E. Lenski on February 24, 1988, at UC Irvine (later moved to Michigan State University). Key features:

- **12 populations** of *E. coli* B strain REL606
- **Strict isolation**: Each population maintained in separate flask
- **Daily transfers**: 1:100 dilution into fresh DM25 medium
- **~6.64 generations per day**
- **Frozen archive**: Samples frozen every 500 generations

### Population Metadata

| Population | Mutator Status | Onset (generations) |
|------------|---------------|---------------------|
| Ara-1 | MMR (mutL) | ~26,000 |
| Ara-2 | MMR (mutS) | ~7,000 |
| Ara-3 | Non-mutator | - |
| Ara-4 | Non-mutator | - |
| Ara-5 | Non-mutator | - |
| Ara-6 | MMR | ~36,000 |
| Ara+1 | MMR | ~2,500 |
| Ara+2 | Non-mutator | - |
| Ara+3 | MMR | ~1,000 |
| Ara+4 | Non-mutator | - |
| Ara+5 | Non-mutator | - |
| Ara+6 | MutT | ~30,000 |

### Data Format

The JSON file contains mutation timing data structured as:

```json
{
  "gene_name": {
    "population_id": first_detection_generation,
    ...
  },
  ...
}
```

Example:
```json
{
  "spoT": {
    "Ara-1": 1500,
    "Ara-2": 2000,
    ...
  }
}
```

### Citation

If using this data, please cite:

1. The original LTEE papers (Good et al. 2017, Tenaillon et al. 2016)
2. This protocol repository

## License

The compiled data is provided under CC BY 4.0, consistent with the original publications' data sharing policies.
