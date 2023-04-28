import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import BSVasprun, Vasprun
from pymatgen.io.vasp.outputs import Oszicar
import veusz.plugins as plugins

def branch2labels(branch: str):
    convert = lambda label: '\Gamma' if label == 'GAMMA' else label
    return [convert(label) for label in branch.split('-')]

class ImportPluginBandStructure(plugins.ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "Band Structure plugin"
    author = "Leran Lu"
    description = "Reads the band structure from vasprun.xml"

    # Uncomment this line for the plugin to get its own tab
    #promote_tab='Example'

    file_extensions = set(['.xml'])

    def __init__(self):
        plugins.ImportPlugin.__init__(self)
        self.fields = [
            plugins.ImportFieldCheck("import_fermi", descr="Import Fermi energy", default=True),
            plugins.ImportFieldCheck("sub_fermi", descr="Substract Fermi energy"),
            plugins.ImportFieldCheck('details', descr='Detailed information')
        ]

    def doImport(self, params: plugins.ImportPluginParams):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """
        datasets = []
        vr = BSVasprun(params.filename, parse_projected_eigen=True)
        bs = vr.get_band_structure()
        efermi = bs.efermi
        if params.field_results['import_fermi']:
            datasets.append(plugins.ImportDataset1D('efermi', [efermi]))
        if not params.field_results['sub_fermi']:
            efermi = 0

        nbands = bs.bands[Spin.up].shape[0]
        breaks = [branch['start_index'] for branch in bs.branches]
        distances = np.insert(bs.distance, breaks, np.nan, axis=0).reshape((1,-1))
        distances = np.repeat(distances, nbands, axis=0).flatten()
        datasets.append(plugins.ImportDataset1D('distances', distances))

        for spin in bs.bands.keys():
            name = 'up' if spin == Spin.up else 'dw'
            dat = bs.bands[spin]-efermi
            dat = np.insert(dat, breaks, np.nan, axis=1).flatten()
            datasets.append(plugins.ImportDataset1D('bands_'+name, dat))

        dist = bs.distance[bs.branches[0]['start_index']]
        last_left, last_right = branch2labels(bs.branches[0]['name'])
        tickd, tickl = [dist], [last_left]
        for i, branch in enumerate(bs.branches):
            dist = bs.distance[branch['end_index']]
            tickd.append(dist)
            left, right = branch2labels(branch['name'])
            if i > 0 and last_right != left:
                tickl[-1] = last_right + '|' + left
            tickl.append(right)
            last_left, last_right = left, right
        datasets += [
            plugins.ImportDataset1D('tickd', tickd),
            plugins.ImportDatasetText('tickl', tickl)
        ]

        if params.field_results['details']:
            datasets += [
                plugins.ImportDataset1D('nbands', [nbands]),
                plugins.ImportDataset1D('distances1', bs.distance)
            ]

        return datasets

class ImportPluginDOS(plugins.ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "DOS plugin"
    author = "Leran Lu"
    description = "Reads the DOS and PDOS from vasprun.xml"

    # Uncomment this line for the plugin to get its own tab
    #promote_tab='Example'

    file_extensions = set(['.xml'])

    def __init__(self):
        plugins.ImportPlugin.__init__(self)
        self.fields = [
            plugins.ImportFieldCheck("import_epdos", descr="Import Elementwise PDOS", default=True),
            plugins.ImportFieldCombo("efermi_style", descr="Fermi energy style", items=['Direct', 'Non-zero'], default='Non-zero', editable=False),
            plugins.ImportFieldCheck("import_fermi", descr="Import Fermi energy", default=True),
            plugins.ImportFieldCheck("sub_fermi", descr="Substract Fermi energy")
            ]

    def doImport(self, params: plugins.ImportPluginParams):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """
        datasets = []
        vr = Vasprun(params.filename, parse_projected_eigen=True)
        dos = vr.complete_dos
        efermi = dos.efermi
        if params.field_results['efermi_style'] == 'Non-zero':
            efermi = dos.energies[np.logical_and(dos.energies<dos.efermi, dos.densities[Spin.up]>0)][-1]
        if params.field_results['import_fermi']:
            datasets.append(plugins.ImportDataset1D('efermi', [efermi]))
        if not params.field_results['sub_fermi']:
            efermi = 0

        datasets.append(plugins.ImportDataset1D('energies', dos.energies-efermi))
        for spin in dos.densities.keys():
            name = 'tdos_up' if spin == Spin.up else 'tdos_dw'
            datasets.append(plugins.ImportDataset1D(name, dos.densities[spin]))

        edos = dos.get_element_dos()
        for elem in edos.keys():
            odos = dos.get_element_spd_dos(elem)
            for spin in edos[elem].densities.keys():
                sspin = 'up' if spin == Spin.up else 'dw'
                datasets.append(plugins.ImportDataset1D('pdos_'+str(elem)+'_'+sspin, edos[elem].densities[spin]))
                for k in odos.keys():
                    datasets.append(plugins.ImportDataset1D('pdos_'+str(elem)+'_'+str(k)+'_'+sspin, data=odos[k].densities[spin]))

        return datasets

class ImportPluginOszicar(plugins.ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "OSZICAR plugin"
    author = "Leran Lu"
    description = "Reads energies of ionic steps from VASP OSZICAR"

    # Uncomment this line for the plugin to get its own tab
    #promote_tab='Example'

    file_extensions = set(['*'])

    def __init__(self):
        plugins.ImportPlugin.__init__(self)
        self.fields = [
            plugins.ImportFieldText('quantities', descr='Quantities (e.g. "E0 dE")'),
            plugins.ImportFieldCheck('indices', descr='Create Indices', default=True),
            plugins.ImportFieldCheck('sub_final', descr="Substract final energy"),
            ]

    def doImport(self, params: plugins.ImportPluginParams):
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
                datasets.append(plugins.ImportDataset1D(quantity, dataset))

        if params.field_results['indices']:
            datasets = [plugins.ImportDataset1D('indices', np.arange(nsteps)+1)] + datasets

        return datasets

plugins.importpluginregistry += [
    ImportPluginBandStructure,
    ImportPluginDOS,
    ImportPluginOszicar
]