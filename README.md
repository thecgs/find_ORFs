# find_ORFs

```

usage: find_ORFs -i [input_file] [-o [output_file]] [-outfmt [int]] [-phase [int]] [-strand [int]] [-min_len [int]] [-max_len [int]] [-translate] [-codontable [int]]
                 [-stop_codons [str] [[str] ...]] [-start_codons [str] [[str] ...]] [-start_codon_model [int]] [-remove_stop_codon] [-remove_nested] [-h] [-v]

Find open Reading Frames (ORFs)

required arguments:
  -i [input_file], --input [input_file]
                        A file of fasta format

optional arguments:
  -o [output_file], --output [output_file]
                        A file of output. defualt: stdout.
  -outfmt [int], --outfmt [int]
                        Output file format. 0: fasta; 1: tsv; 2: gff. defualt: 0.
  -phase [int], --phase [int]
                        Start address of sequence. 0:all; 1: first base; 2: second base; 3: third base. defualt: 0.
  -strand [int], --strand [int]
                        Search strand, 0: both; 1: +; 2: -. defualt: 0.
  -min_len [int], --min_len [int]
                        Min Length. defualt: 0.
  -max_len [int], --max_len [int]
                        Max Length. defualt: INF.
  -translate, --translate
                        Translate DNA to protein sequence. defualt: False.
  -codontable [int], --codontable [int]
                        Genetic code to use. choices: 1-6, 9-14, 16, 21-31, 33. defualt: 1.
  -stop_codons [str] [[str] ...], --stop_codons [str] [[str] ...]
                        Stop codons. defualt: None.
  -start_codons [str] [[str] ...], --start_codons [str] [[str] ...]
                        Start codons. defualt: None.
  -start_codon_model [int], --start_codon_model [int]
                        ORF start codon to use. 0: only "ATG"; 1: "ATG" and alternative initiation codons. defualt: 0.
  -remove_stop_codon, --remove_stop_codon
                        Remove portein seqence stop codon. defualt: False.
  -remove_nested, --remove_nested
                        Ignore nested ORFs:. defualt: False.
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

date:2023/02/05 author:guisen chen email:thecgs001@foxmail.com

```
