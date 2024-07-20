# tests/__init__.py

import pytest

# Example fixture that can be used across multiple test modules
@pytest.fixture(scope='module')
def common_data():
    return {
        'isotypes': ['IGHM', 'IGHA2', 'IGHG2', 'IGHG1', 'IGHA1', 'IGHG4', 'IGHG3', 'IGHE'],
         'clonotype_list' : ['Clonotype_1036', 'Clonotype_1050', 'Clonotype_10884']
    }