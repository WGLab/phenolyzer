# Phenolyzer
Phenolyzer is a software that implements phenotype-based prioritization of candidate  genes for human diseases.

## Introduction
Prior biological knowledge and phenotype information may help pinpoint disease contributory genes in whole genome/exome sequencing studies on human diseases. We developed a computational tool called Phenolyzer, which follows a biologist's natural thought processes through four steps: term interpretation, seed gene generation, seed gene growth and data integration. Compared to competing approaches, Phenolyzer has superior performance on finding known disease genes, and on prioritizing recently published novel disease genes.

## Releases
Phenolyzer is also available as a web-sever at [here](http://phenolyzer.wglab.org).

## Mannual
The more detailed manual is available at [here](http://phenolyzer.wglab.org/download/Phenolyzer_manual.pdf) 

## Dependency

[Bioperl](http://www.bioperl.org/wiki/Main_Page)

Graph::Directed

## Installation

Please clone the repository into your computer:
```
git clone https://github.com/WGLab/phenolyzer
```
Then enter phenolyzer directory:
```
cd phenolyzer
```

## Use your own Addon databases

Put your own addon databases into lib/compiled_database.
For details of how to use your own Addon databases, please refer to [FAQ](http://phenolyzer.wglab.org/FAQ.php#collapse-14)

Use 
```
-addon [Gene-disease database file name, delimited by comma] -addon_weight [default:1] -addon_gg [Gene-gene database file name, delimited by comma] -addon_gg_weight [default:1]
```

## Pre-requisites
- Download the databases for CNV annotation (No need to do this for other functions)
```
perl bin/annotate.pl -downdb -buildver hg19 -webfrom annovar refGene lib/humandb
perl bin/annotate.pl -downdb -buildver hg18 -webfrom annovar refGene lib/humandb
```

## Synopsis

- Print help message
```
perl disease_annotation.pl --help
```

- Prioritize 'sleep' genes: 
```
perl disease_annotation.pl sleep -p -ph -logistic -out out/sleep/out
```

- Use the terms in 'disease' file:
```
perl disease_annotation.pl disease -f -p -ph -logistic -out out/disease/out
```

- Use the cnv.bed region:
```
perl disease_annotation.pl alzheimer -bedfile cnv.bed -p -ph -logistic -out out/alzheimer/out
```

- Use the Mentha gene-gene interaction database as Addon
```
perl disease_annotation.pl alzheimer -p -ph -logistic -out out/alzheimer_addon/out -addon_gg DB_MENTHA_GENE_GENE_INTERACTION -addon_gg_weight 0.05
```

- To generate exactly the same result as Phenolyzer web server default settings
```
perl disease_annotation.pl alzheimer -p -ph -logistic -out out/alzheimer_addon_all/out -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
```

- Integrate with wANNOVAR output to prioritize variant
```
perl calculate_score.pl <phenolyzer_gene_list> <wannovar_genome_summary.txt>
```

- Input multiple diseases (alzheimer and brain)
```
perl disease_annotation.pl "alzheimer;brain" -p -ph -logistic -out out/sd
```

## License Agreement
By using the software, you acknowledge that you agree to the terms below:

For academic and non-profit use, you are free to fork, download, modify, distribute and use the software without restriction.

For commercial use, you are required to contact [Stevens Institute of Innovation](https://stevens.usc.edu/contact-us/) at USC directly to discuss licensing options.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Contact
- Hui Yang (yanghui@usc.edu)
- Kai Wang (kaichop@gmail.com)


