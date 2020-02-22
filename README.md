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

Zamin Iqbal, Mario Caccamo, Isaac Turner, Paul Flicek, Gil McVean.
*De novo assembly and genotyping of variants using colored de Bruijn graphs*
**Nature Genetics**  44, pages226–232(2012) 
[(link)](https://www.nature.com/articles/ng.1028)

## For future reference
We are currently working on a new variant calling package which wraps cortex and samtools, which can be found here
[(link)](https://github.com/iqbal-lab-org/clockwork). We have been building singularity containers of these, and these should shortly be hosted externally and easily available



