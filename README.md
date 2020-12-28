# simple_adl
A new version of the simple binner dwarf galaxy search algorithm ([`simple`](https://github.com/DarkEnergySurvey/simple)) intended for use with data pulled directly from the Astro Data Lab.

## Installation

Clone this repository to the machine you wish to run the code from.

Run the following command (consider adding to your `.bashrc`/`.zshrc`/etc.):
```
export PYTHONPATH=<path to simple_adl>:$PYTHONPATH
```
where `<path to simple_adl>` is the path to the inner simple_adl repository (ie, the full string should end with `*/simple_adl:$PYTHONPATH`).

Create the conda environment with
```
conda env create -f conda.yaml
```
This will create a conda environment called `simple_adl` that includes all of the necessary software.

## Use

Run `init_simple.py` to generate a default `config.yaml` configuration file.
This can be edited as needed; in particular, you may want to change the `profile` entry.

To search for [DELVE 2](https://arxiv.org/abs/2009.08550), run the following:
```
python search.py --ra 28.77 --dec -68.25 
```
or
```
python search.py --nside 32 --ipix 11812
```

To plot DELVE 2, run
```
python plot_hotspot.py --ra 28.77 --dec -68.25 --mod 19.26
```
