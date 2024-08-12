# Installation 
`tribal` can be installed via `bioconda` (recommended) or `github` and can be run as a [python package](package.md) or via the [command line](cli.md). 


### Dependencies
`tribal` has the following dependencies:
```
    python>=3.9
    numpy>=1.26.3
    networkx>=3.1
    ete3>=3.1.2
    glpk>=5.0
    pyomo>=6.7.3
    pygraphviz>=1.10
    pandas>=2.1.1
    biopython>=1.81
    mafft>=7.526
    phylip>=3.697 (included with package)
```

`tribal` can be installed from [Github](#installing-from-github) or from [bioconda](#installing-from-bioconda).  We recommend installing from [bioconda](#installing-from-bioconda).

### Installing from bioconda

```bash
conda create -n tribal -c bioconda tribal
conda activate tribal
```

### Installing from GitHub
```
git clone https://github.com/elkebir-group/TRIBAL.git
cd TRIBAL

```

Dependencies can be installed using a package mangager such as `conda`, `mamba` or `micromamba`, using the included `tribal.yml` file.



```bash
conda create -f tribal.yml 
```

To build and install `tribal` into this environment, follow the instructions below.

```bash
conda activate tribal
pip install hatchling
hatchling build
pip install dist/tribal-0.1.0-py3-none-any.whl

```

### Verifying installation


`tribal` can be imported as a package or run via a  command line interface.  

To verify the package can be imported,  run the following in the terminal.

```bash
python -c "import tribal"
```

See [Package Overivew](package.md) for more detailed usage intstructions. 

To verify the CLI tool was properly installed, run the following in the terminal. 

```bash
tribal --help
```

See [Command Line Interface](cli.md) for more detailed usage intstructions.


