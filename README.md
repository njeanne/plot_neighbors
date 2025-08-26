# Plot neighbors

From a CSV of the amino acids neighborhood contacts file and a YAML parameters file, produced by the script 
trajectories_neighbors.py (https://github.com/njeanne/trajectories_neighbors), create a heatmap representing the 
residues neighborhood contacts.

A Region Of Interest (ROI) is defined with a range of amino acids selected in the protein, on the heatmap the contacts 
on the ordinate axis will be the ones belonging to this ROI.

If a domains CSV file is used with the option `--domains`, a plot and a CSV file of the neighborhood contacts by domains 
will be produced. An example of the domains CSV file is provided in `data/traj_test_domains.csv`
For this CSV, if some domains are embedded in other domains, they can be displayed in the outputs with the option 
`--embedded-domains`.

## Conda environment

A [conda](https://docs.conda.io/projects/conda/en/latest/index.html) YAML environment file is provided: 
`conda_env/python3_env.yml`. The file contains all the dependencies to run the script.
The conda environment is generated using the command:
```shell script
# create the environment
conda env create -f conda_env/python3_env.yml

# activate the environment
conda activate python3
```

## Usage

The script can be tested with the test data provided in the `data` directory, which contains the output CSV and the 
YAML parameters files of the [trajectories_neighbors](https://github.com/njeanne/trajectories_neighbors) script, 
respectively `AB248520-3e_WT_ORF1_neighbors.csv` and `AB248520-3e_WT_ORF1_analysis.yaml`. The commands are:

```shell script
conda activate python3

./plot_neighbors.py --out results --roi 682-840 --format svg\
--parameters data/AB248520-3e_WT_ORF1_analysis.yaml data/AB248520-3e_WT_ORF1_neighbors.csv

conda deactivate
```

An optional domains CSV file can also be provided with the `--domains` option, `data/AB248520-3e_WT_ORF1_domains.csv`. 
The commands are:

```shell script
conda activate python3

./plot_neighbors.py --out results --roi 682-840  --format svg \
--domains data/AB248520-3e_WT_ORF1_domains.csv --residues-distance 10 \
--parameters data/AB248520-3e_WT_ORF1_analysis.yaml  data/AB248520-3e_WT_ORF1_neighbors.csv

conda deactivate
```

The optional parameter used:
- `--domains`: to produce a number of neighborhood contacts plot of the whole protein or the region of interest if used 
on the whole protein.
- `--residues-distance`: the minimal number of residues between two residues to validate a contact.

## Outputs

The script outputs are:
- a heatmap of the neighborhood contacts:

![neighbors heatmap](doc/_static/heatmap.svg)

And if the `--domains` option is used: 
- a CSV file of the neighborhood contacts with the minimal residues distance for the whole protein or the region of interest if used 
on the whole protein.
- a plot of the neighborhood contacts with the minimal residues distance for the whole protein or the region of interest if used on 
the whole protein:

![neighbors heatmap](doc/_static/outliers.svg)