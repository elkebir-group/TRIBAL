# Installation 



### Dependencies
`tribal` has the following dependencies:
```
    python==3.9
    numpy=1.26.3
    networkx==3.1
    ete3==3.1.2
    glpk==5.0
    pyomo==6.7.3
    pygraphviz>=1.10
    pandas==2.1.1
    biopython==1.81
    mafft>=7.526
    phylip>=3.697
```

`tribal` can be installed from [Github](#installing-from-github) or from [bioconda](#installing-from-bioconda).  We recommend installing from [bioconda](#installing-from-bioconda).

### Installing from bioconda

```
conda create -n tribal -c bioconda tribal
conda activate tribal
```

### Installing from GitHub
```
git clone https://github.com/elkebir-group/TRIBAL.git

```



Dependencies can be installed using a package mangager such as `conda` `mamba` or `micromamba`, using the included `tribal.yml` file.

#### Using conda package manager to install dependencies


```
conda create -f tribal.yml 

```

*Note: `tribal` is not currently compatible arm64 architecture. Please using the following command instead if you are using an M1 Macbook.

```
CONDA_SUBDIR=osx-64 conda env create -f tribal.yml --solver libmamba

```

To build and install `tribal` into this environment, follow the instructions below.

```
conda activate tribal
pip install hatchling
hatchling build
pip install dist/tribal-0.1.0-py3-none-any.whl

```

#### Using micromamba (mamba) package manager to install dependencies


```
micromamba create -f tribal.yml 

```

*Note: `tribal` is not currently compatible arm64 architecture. Please using the following command instead if you are using an M1 Macbook.

```
micromamba create -f tribal.yml --platform osx-64 

```

#### Build and install the `tribal` package

To build and install `tribal` into this environment, follow the instructions below, replacing `conda` with `mamba` or `micromamba` as necessary. 

```
conda activate tribal
pip install hatchling
hatchling build
pip install dist/tribal-0.1.0-py3-none-any.whl 

```

`tribal` can be imported as a package or run via a  command line interface

```
import tribal
```

```
tribal preprocess --help
```

```
tribal fit --help
```

See [Package Overivew](package.md) for more detailed usage intstructions. 


