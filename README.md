# cortex

Reference-free variant assembly

## Installation

```
git clone --recursive https://github.com/iqbal-lab/cortex.git
cd cortex
bash install.sh
make cortex_var
```

By default, `cortex` will compile a binary called `cortex_var_31_c1`.
This supports a maximum k-mer size of 31, and up to 1 colour in the graph.
To increase these, you need to re-compile `cortex` as follows:
```
make cortex_var MAXK=127 NUM_COLS=8
```
The `MAXK` can only take values of the form `32 x N - 1`, in the range 31 to 255.
The `NUM_COLS` parameter must be a positive integer.

## Dependencies

* `htslib` (bundled)
* `seq_file` (bundled)
* `string_buffer` (bundled)
* `zlib`

## Usage

Build a single-colour binary:
```
cortex_var --se_list <filename> --pe_list <filename> --format FASTQ \
  --quality_score_threshold 5 --remove_pcr_duplicates \
  --remove_low_coverage_supernodes 1 --dump_binary some_name.ctx
```

Build a multicolour graph from single-colour graphs and call variants between colours 1 and 2:
```
cortex_var --colour_list <filename> --detect_bubbles1 1/2 \
  --output_bubbles1 vars_between_cols1_and_2
```

Load a multicolour graph from single-colour graphs and call heterozygous variants in colour 0:
```
cortex_var --colour_list <filename> --detect_bubbles1 0/0 \
  --output_bubbles1 hets_in_colour_0
```

## Licence

[GPLv3](https://raw.githubusercontent.com/iqbal-lab/cortex/master/gpl.txt)

## Feedback

File questions or issues to the [Issue Tracker](https://github.com/iqbal-lab/cortex/issues)

## Citation

Zamin Iqbal, Isaac Turner, Gil McVean.
*High-throughput microbial population genomics using the Cortex variation assembler*
**Bioinformatics**  2013;2(15):275-6. 
[(link)](https://doi.org/10.1093/bioinformatics/bts673)
