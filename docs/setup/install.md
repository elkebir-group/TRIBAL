# Installation 
`tribal` can be installed via `bioconda` (coming soon) or `github` and can be run as a [python package](package.md) or via the [command line](cli.md). 


### Dependencies
`tribal` has the following dependencies:
```
    python >=3.9,<3.11
    numpy >=1.26,<2.0
    pandas
    networkx >=3.1
    pygraphviz >=1.10
    ete3 >=3.1.2
    mafft ==7.526
    glpk >=5.0
    pyomo >=6.7
    biopython >=1.81
```



`tribal` can be installed from [Github](#installing-from-github) or from [bioconda](#installing-from-bioconda).  We recommend installing from [bioconda](#installing-from-bioconda).

### Installing from bioconda (recommended)


```bash
conda create -n tribal -c bioconda tribal
conda activate tribal
```
See [Veryifing installation](#verifying-installation) to make sure the package was installed correctly. 

### Installing from GitHub
```
git clone https://github.com/elkebir-group/TRIBAL.git
cd TRIBAL

```

Dependencies can be installed using a package mangager such as `conda`, `mamba` or `micromamba`, using the included `tribal.yml` file.



```bash
conda env create -f tribal.yml 
```

To build and install `tribal` into this environment, follow the instructions below:

```bash
conda activate tribal
python -m build 
pip install dist/tribal-0.1.0-py3-none-any.whl 

```

`dnapars` is a requirement to run the preprocessing tools. The source code is included in the package and can be installed using the `dnapars_install.sh` script. 

```bash
chmod +x dnapars_install.sh
./dnapars_install.sh

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


