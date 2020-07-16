# SVAFotate
Annotate Structural Variants with Population Allele Frequency Info 



## Installing SVAFotate


0) Installing Miniconda

- If Miniconda is not installed on your system, install it from [miniconda](https://conda.io/en/latest/miniconda.html)


1) Set up new conda environment 

```
$ conda create --name svafotate-env python=3
```

```
$ conda activate svafotate-env
```


2) Install package requirements 

```
$ conda install --file https://raw.githubusercontent.com/mikecormier/SVAFotate/master/requirements.txt?token=AIM3FBTOHCDMZJYYLOX3OYC7B7FGW
```


3) Install SVAFotate

```
$ pip install -U git@github.com:mikecormier/SVAFotate.git
```


4) Check that SVAFotate installed Correctly 

```
$ svafotate --version

svafotate 0.0.1
```

```
$ svafotate -h 

usage: svafotate [-h] [-v] {annotate,pickle-source,custom-annotation} ...

SVAFotate: Structural Variant VCF Annotation tools
==================================================

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Installed version

[sub-commands]:
  {annotate,pickle-source,custom-annotation}
    annotate            Annotate SV VCF File
    pickle-source       Pickle Source Bed
    custom-annotation   Add custom annotation(s) to source annotation file
```

