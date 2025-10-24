# Copilot Instructions for pafkit

## Project Overview
pafkit is a Python toolkit for processing and visualizing PAF (Pairwise mApping Format) alignment files, specifically designed for whole-genome synteny analysis between species. The toolkit focuses on comparative genomics workflows where genomes are aligned using minimap2 and then processed for visualization and structural variant detection.

## Core Architecture

### Module Structure
- **modules/paf.py**: Central `PAF` class that parses alignment files and manages collections of `Alignment` objects
- **modules/genome.py**: `Genome` class handles sequence metadata (chromosome names, lengths, cumulative positions)  
- **modules/aln.py**: `Alignment` class represents individual alignment records with query/target coordinates
- **modules/homology_blocks.py**: Groups alignments into syntenic blocks for fission/fusion detection
- **modules/pafplotter.py**: Matplotlib-based plotting utilities for synteny visualization

### Data Flow Pattern
1. Parse genome index files (.fai format) into `Genome` objects
2. Parse PAF alignment files into `PAF` objects containing `Alignment` collections
3. Filter alignments by length, identity, mapping quality thresholds
4. Apply coordinate transformations for visualization (cumulative positioning, strand reorientation)
5. Generate matplotlib plots or export filtered PAF data

## Key Conventions

### Coordinate Systems
- **Internal coordinates**: Always 0-based, half-open intervals throughout the codebase
- **User input**: 1-based, inclusive intervals (converted internally via `start - 1`)
- **Strand handling**: Negative strand alignments automatically reorient query coordinates during plotting

### File Format Expectations
- **Genome files**: Tab-delimited with sequence_name and length columns (samtools .fai format preferred)
- **PAF files**: Standard minimap2 output format, supports stdin input via `-i -`
- **Config files**: Custom format for multi-genome progressive synteny (see `config_hominins.txt`)

### Filtering Pipeline
All scripts use consistent filtering parameters: `--min-aln-length`, `--min-identity`, `--min-quality`, `--min-reference-length`, `--min-query-length`. Apply filters before visualization to avoid rendering performance issues with large datasets.

## Essential Workflows

### Basic Synteny Plotting
```bash
# Generate alignments with minimap2
minimap2 -x asm5 -c reference.fa query.fa > alignments.paf

# Create genome index files  
samtools faidx reference.fa
samtools faidx query.fa

# Plot synteny map
./plot_paf.py -i alignments.paf \
  -r reference.fa.fai -q query.fa.fai \
  -o output.png --min-aln-length 20000
```

### PAF Subsetting and Chaining
```bash
# Subset to specific chromosome and filter low-coverage sequences
./subset_paf.py -i alignments.paf -r chr1 \
  --aligned-fraction 0.2 \
  --query query.fa.fai | \
./plot_paf.py -i - -o chr1_plot.png \
  -r reference.fa.fai -q query.fa.fai
```

## Development Patterns

### Adding New Filtering Options
1. Add argument to script's `parse_args()` function
2. Implement filter logic in the main processing loop
3. Apply filter before calling `paf.filter_*()` methods
4. Test with various genome size datasets

### Visualization Customization
- Modify `synteny_plot()` function in plotting scripts
- Colors and styling defined as function parameters (ref_color, query_color, aln_color)
- Matplotlib patches used: `Rectangle` for chromosomes, `Polygon` for alignment ribbons
- Y-axis positioning: reference at y=2, query at y=1, ribbons between

### Memory Considerations
- Large PAF files (>1GB) should be processed in streaming mode
- Use `--min-aln-length` aggressively to reduce memory footprint
- Consider chromosome-by-chromosome processing for whole-genome comparisons

## Testing and Debugging

### Validate Alignments
```bash
# Check alignment parsing without plotting
python3 -c "import modules.paf as paf; p = paf.PAF(); p.parse_paf('alignments.paf'); print(p)"

# Test filtering thresholds
./subset_paf.py -i alignments.paf -l 50000 --min-identity 0.95 -o /dev/null
```

### Common Issues
- **Coordinate mismatches**: Verify 0-based vs 1-based coordinate handling in new functions
- **Strand orientation**: Check `target_sequence_strand` field processing in alignment visualization
- **Memory usage**: Profile with large PAF files, implement streaming if needed
- **Plot readability**: Adjust `gap-size` parameter for chromosome spacing in dense plots

## Dependencies and Environment
- **Core**: matplotlib, numpy (see requirements.txt)
- **External tools**: minimap2 for alignment generation, samtools for genome indexing
- **Python version**: Developed with 3.12.7, compatible with most 3.x versions
- **File permissions**: Make scripts executable with `chmod +x *.py`