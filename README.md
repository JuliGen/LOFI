# 👩🏻‍💻OperonMapper <img src ='https://papik.pro/uploads/posts/2021-12/1639240390_33-papik-pro-p-dinozavr-klipart-33.png' width =400 align="right">
> *This repository is for a project to create a tool for searching for operons*

## Instruction: 
### `calculate_dist_score`🧑🏻‍🔬 
This function allows you to calculate the intergenic distance, as well as calculate the metric value of a gene’s belonging to an operon.
The function takes 2 arguments as input. 

- `input_gff` - path to gff file (after the bakta annotation)
- `output_file` - output file name.

As a result, the `operon_map_result` folder is created, which contains the annotation file. The file includes additional columns such as intergenic_distance_next, intergenic_distance_prev, score_operon

#### Example:
```python
from Code.operon_score import calculate_dist_score

calculate_dist_score(input_gff = 'Data/annot_bakta.gff3')

#OUTPUT: 'Created annot_bakta_calc_dist_score.txt'

calculate_dist_score(input_gff = 'Data/annot_bakta.gff3', output_file = 'operon_dist_score.txt')

#OUTPUT: 'Created operon_dist_score.txt'

```
