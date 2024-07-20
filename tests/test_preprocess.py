# test_preprocess.py
import pytest
import pandas as pd
from tribal.preprocess import process_clonotype, preprocess

def test_process_clonotype():
    df = pd.read_csv("tribal/data/data.csv")
    roots = pd.read_csv("tribal/data/roots.csv")

    j = "Clonotype_1036"
    use_light_chain = True
    verbose = False

    j, linforest = process_clonotype(j, df, roots, use_light_chain, verbose)
    assert j == "Clonotype_1036"
    assert linforest is not None

def test_preprocess_single_core():
    df = pd.read_csv("tribal/data/data.csv")
    roots = pd.read_csv("tribal/data/roots.csv")

   
    df = df[df['clonotype'].isin(common_data["clonotype_list"])]

    clonotypes, dat_filt = preprocess(df, roots, common_data["isotypes"], cores=1, verbose=False)

    assert len(clonotypes) > 0
    assert not dat_filt.empty


def test_preprocess_multi_core():
    df = pd.read_csv("tribal/data/data.csv")
    roots = pd.read_csv("tribal/data/roots.csv")

    df = df[df['clonotype'].isin(common_data["clonotype_list"])]

    isotypes = ['IGHM', 'IGHA2', 'IGHG2', 'IGHG1', 'IGHA1', 'IGHG4', 'IGHG3', 'IGHE']
    clonotypes, dat_filt = preprocess(df, roots, common_data["isotypes"], cores=3, verbose=False)

    assert len(clonotypes) > 0
    assert not dat_filt.empty


if __name__ == "__main__":
    pytest.main()




