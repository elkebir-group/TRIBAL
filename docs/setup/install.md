# Installation 

`tribal` can be installed from [GitHub](https://github.com/elkebir-group/TRIBAL) or from `bioconda`.  We recommend installing from `bioconda`.

### Installing from Bioconda

```
conda create -n tribal -c bioconda tribal
```

### Installing from GitHub

```
git clone https://github.com/elkebir-group/TRIBAL.git

```

`tribal` has the following dependencies.
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

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages images and other files.
