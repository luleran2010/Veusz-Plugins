import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import BSVasprun
from veusz.plugins import *


class ImportPluginBandStructure(ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "Band Structure plugin"
    author = "Leran Lu"
    description = "Reads the band structure from vasprun.xml"

    # Uncomment this line for the plugin to get its own tab
    #promote_tab='Example'

    file_extensions = set(['*'])

    def __init__(self):
        ImportPlugin.__init__(self)
        self.fields = [
            ImportFieldCheck("import_fermi", descr="Import Fermi energy", default=True),
            ImportFieldCheck("sub_fermi", descr="Substract Fermi energy"),
            ImportFieldCheck('details', descr='Detailed Information')
        ]

    def doImport(self, params: ImportPluginParams):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """
        datasets = []
        vr = BSVasprun(params.filename, parse_projected_eigen=True)
        bs = vr.get_band_structure()
        efermi = bs.efermi
        if params.field_results['import_fermi']:
            datasets.append(ImportDataset1D('efermi', [efermi]))
        if not params.field_results['sub_fermi']:
            efermi = 0

        nbands = bs.bands[Spin.up].shape[0]
        distances = np.array(bs.distance).reshape((1,-1)) 
        distances = np.repeat(distances, nbands, axis=0)
        distances = np.concatenate((distances, np.array([np.nan]*nbands).reshape((nbands,1))), axis=1).flatten()
        datasets.append(ImportDataset1D('distances', distances))

        for spin in bs.bands.keys():
            name = 'up' if spin == Spin.up else 'dw'
            dat = bs.bands[spin]-efermi
            dat = np.concatenate((dat, np.array([np.nan]*nbands).reshape((nbands,1))), axis=1).flatten()
            datasets.append(ImportDataset1D('energies_'+name, dat))

        tickd = []
        tickl = []
        for branch in bs.branches:
            if branch['start_index'] == 0:
                tickd.append(distances[branch['start_index']])
                tickl.append(branch['name'].split('-')[0])
                if tickl[-1] == 'GAMMA':
                    tickl[-1] = '\Gamma'
            tickd.append(distances[branch['end_index']])
            tickl.append(branch['name'].split('-')[1])
            if tickl[-1] == 'GAMMA':
                tickl[-1] = '\Gamma'
        datasets += [
            ImportDataset1D('tickd', tickd),
            ImportDatasetText('tickl', tickl)
        ]

        if params.field_results['details']:
            datasets += [
                ImportDataset1D('nbands', [nbands]),
                ImportDataset1D('distances1', bs.distance)
            ]

        return datasets

# add the class to the registry. An instance also works, but is deprecated
importpluginregistry.append(ImportPluginBandStructure)