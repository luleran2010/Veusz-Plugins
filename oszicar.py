import numpy as np
from pymatgen.io.vasp.outputs import Oszicar
from veusz.plugins import *


class ImportPluginOszicar(ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "OSZICAR plugin"
    author = "Leran Lu"
    description = "Reads energies of ionic steps from VASP OSZICAR"

    # Uncomment this line for the plugin to get its own tab
    #promote_tab='Example'

    file_extensions = set(['*'])

    def __init__(self):
        ImportPlugin.__init__(self)
        self.fields = [
            ImportFieldText('quantities', descr='Quantities (e.g. "E0 dE")'),
            ImportFieldCheck('indices', descr='Create Indices', default=True),
            ImportFieldCheck('sub_final', descr="Substract final energy"),
            ]

    def doImport(self, params: ImportPluginParams):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """
        datasets = []
        ionic_steps = Oszicar(params.filename).ionic_steps
        nsteps = len(ionic_steps)
        quantities = [i.strip() for i in str(params.field_results['quantities']).split()]
        for quantity in quantities:
            if quantity in ionic_steps[0]:
                dataset = np.array([i[quantity] for i in ionic_steps])
                if params.field_results['sub_final'] and quantity in ['E0', 'F']:
                    dataset -= dataset[-1]
                datasets.append(ImportDataset1D(quantity, dataset))

        if params.field_results['indices']:
            datasets = [ImportDataset1D('indices', np.arange(nsteps)+1)] + datasets

        return datasets

# add the class to the registry. An instance also works, but is deprecated
importpluginregistry.append(ImportPluginOszicar)