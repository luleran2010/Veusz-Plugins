import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun
from veusz.plugins import *


class ImportPluginDOS(ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "DOS plugin"
    author = "Leran Lu"
    description = "Reads the DOS and PDOS from vasprun.xml"

    # Uncomment this line for the plugin to get its own tab
    #promote_tab='Example'

    file_extensions = set(['*'])

    def __init__(self):
        ImportPlugin.__init__(self)
        self.fields = [
            ImportFieldCheck("import_epdos", descr="Import Elementwise PDOS", default=True),
            ImportFieldCombo("efermi_style", descr="Fermi energy style", items=['Direct', 'Non-zero'], default='Non-zero', editable=False),
            ImportFieldCheck("import_fermi", descr="Import Fermi energy", default=True),
            ImportFieldCheck("sub_fermi", descr="Substract Fermi energy")
            ]

    def doImport(self, params: ImportPluginParams):
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
            datasets.append(ImportDataset1D('efermi', [efermi]))
        if not params.field_results['sub_fermi']:
            efermi = 0

        datasets.append(ImportDataset1D('energies', dos.energies-efermi))
        for spin in dos.densities.keys():
            name = 'dos_up' if spin == Spin.up else 'dos_dw'
            datasets.append(ImportDataset1D(name, dos.densities[spin]))

        edos = dos.get_element_dos()
        for elem in edos.keys():
            odos = dos.get_element_spd_dos(elem)
            for spin in edos[elem].densities.keys():
                sspin = 'up' if spin == Spin.up else 'dw'
                datasets.append(ImportDataset1D('epdos_'+str(elem)+'_'+sspin, edos[elem].densities[spin]))
                for k in odos.keys():
                    datasets.append(ImportDataset1D('eopdos_'+str(elem)+'_'+sspin, data=odos[k].densities[spin]))

        return datasets

# add the class to the registry. An instance also works, but is deprecated
importpluginregistry.append(ImportPluginDOS)