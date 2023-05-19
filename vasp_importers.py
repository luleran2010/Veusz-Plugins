import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.outputs import BSVasprun, Vasprun
from pymatgen.io.vasp.outputs import Oszicar
import veusz.plugins as plugins

import os

class MyBandStructure:
    def __init__(self, filename: str, hybrid: bool, kpoints_fn: str='', kpath_fn: str='') -> None:
        if not hybrid:
            self._read_pbe(filename)
        else:
            self._read_hybrid(filename, kpoints_fn, kpath_fn)
        
    def _read_pbe(self, filename: str):
        vr = BSVasprun(filename, parse_projected_eigen=True)
        bs = vr.get_band_structure()
        self.efermi = bs.efermi
        self.nbands = bs.nb_bands
        self.breaks = [branch['start_index'] for branch in bs.branches]
        self.end_indices = [branch['end_index'] for branch in bs.branches]
        self.distances = bs.distance
        self.bands = {spin: bands for spin, bands in bs.bands.items()}

        convert = lambda label: '\Gamma' if label == 'GAMMA' else label
        self.branches = [[convert(i) for i in branch['name'].split('-')] for branch in bs.branches]

    def _read_hybrid(self, filename: str, kpoints_fn: str, kpath_fn: str):
        vr = BSVasprun(filename, parse_projected_eigen=True)
        bs = vr.get_band_structure()
        self.efermi = bs.efermi
        self.nbands = bs.nb_bands
        
        kpoints = Kpoints.from_file(kpoints_fn)
        terms = kpoints.comment.split()
        nkpts = int(terms[6])
        nbranches = int(terms[7])
        npts = [int(term) for term in terms[8:8+nbranches]]
        beg = int(terms[4])
        kpts = np.array([kp.cart_coords for kp in bs.kpoints[beg:]])
        self.breaks = np.array([np.sum(npts[0:i]) for i in np.arange(len(npts))], dtype=int)
        self.end_indices = np.array([np.sum(npts[0:i])-1 for i in np.arange(1, len(npts)+1)])

        previous_distance = 0.0
        distance = []
        head = 0
        for length in npts:
            diff = np.linalg.norm(np.diff(kpts[head:head+length], axis=0), axis=1)
            dist = np.array([np.sum(diff[0:i]) for i in np.arange(len(diff)+1)]) + previous_distance
            distance.append(dist)
            previous_distance = dist[-1]
            head += length
        self.distances = np.concatenate(distance, axis=0)

        self.bands = {spin: bands[:,beg:] for spin, bands in bs.bands.items()}

        kpath = Kpoints.from_file(kpath_fn)
        convert = lambda label: '\Gamma' if label == 'GAMMA' else label
        self.branches = [[convert(kpath.labels[2*i]), convert(kpath.labels[2*i+1])] for i in np.arange(len(kpath.labels)/2, dtype=int)]

    @staticmethod
    def parse_path(s: str) -> list:
        terms = s.split('-')
        branches = []
        last_term = terms[0]
        for cur_term in terms[1:]:
            if '|' in cur_term:
                l, r = cur_term.split('|')
                branches.append([last_term, l])
                last_term = r
            else:
                branches.append([last_term, cur_term])
                last_term = cur_term
        return branches

    def change_path(self, s: str) -> bool:
        branches = MyBandStructure.parse_path(s)
        indices = []
        flip = []
        for branch in branches:
            found = False
            for i in np.arange(len(self.branches)):
                if self.branches[i][0] == branch[0] and self.branches[i][1] == branch[1]:
                    indices.append(i)
                    flip.append(False)
                    found = True
                    break
                elif self.branches[i][0] == branch[1] and self.branches[i][1] == branch[0]:
                    indices.append(i)
                    flip.append(True)
                    found = True
                    break
            if not found:
                return False
        
        distances = []
        bands = {spin: [] for spin in self.bands.keys()}
        npts = []
        self.branches = branches
        last_distance = 0
        for i, f in zip(indices, flip):
            distance = self.distances[self.breaks[i]:self.end_indices[i]+1]
            if f:
                distance = np.flip(distance)
            distance = np.abs(distance - distance[0]) + last_distance
            last_distance = distance[-1]
            distances.append(distance)
            for key in bands.keys():
                band = self.bands[key][:,self.breaks[i]:self.end_indices[i]+1]
                bands[key].append(band if not f else np.fliplr(band))
            npts.append(len(distance))
        self.distances = np.concatenate(distances)
        for key in bands.keys():
            bands[key] = np.concatenate(bands[key], axis=1)
        self.bands = bands
        self.breaks = np.array([np.sum(npts[0:i]) for i in np.arange(len(npts))], dtype=int)
        self.end_indices = np.array([np.sum(npts[0:i])-1 for i in np.arange(1, len(npts)+1)], dtype=int)

        return True

class ImportPluginBandStructure(plugins.ImportPlugin):
    """Plugins to import band structure from vasprun.xml"""

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
            plugins.ImportFieldCheck('hybrid', descr='Hybrid functionals', default=False),
            plugins.ImportFieldText('kpath', descr='K-Path (blank for defualt)', default=''),
            plugins.ImportFieldCheck('details', descr='Detailed information'),
        ]

    def doImport(self, params: plugins.ImportPluginParams):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """
        datasets = []

        kpoints_fn = ''
        kpath_fn = ''
        if params.field_results['hybrid']:
            dirname = os.path.dirname(params.filename)
            kpoints_fn = os.path.join(dirname, 'KPOINTS')
            kpath_fn = os.path.join(dirname, 'KPATH.in')
            
        mbs = MyBandStructure(params.filename, hybrid=params.field_results['hybrid'], kpoints_fn=kpoints_fn, kpath_fn=kpath_fn)
        if params.field_results['kpath'] != '':
            mbs.change_path(params.field_results['kpath'])

        efermi = mbs.efermi
        if params.field_results['import_fermi']:
            datasets.append(plugins.ImportDataset1D('efermi', [efermi]))
        if not params.field_results['sub_fermi']:
            efermi = 0

        distances = np.insert(mbs.distances, mbs.breaks, np.nan, axis=0).reshape((1,-1))
        distances = np.repeat(distances, mbs.nbands, axis=0).flatten()
        datasets.append(plugins.ImportDataset1D('distances', distances))

        for spin in mbs.bands.keys():
            name = 'up' if spin == Spin.up else 'dw'
            dat = mbs.bands[spin]-efermi
            dat = np.insert(dat, mbs.breaks, np.nan, axis=1).flatten()
            datasets.append(plugins.ImportDataset1D('bands_'+name, dat))

        dist = mbs.distances[0]
        last_left, last_right = mbs.branches[0]
        tickd, tickl = [dist], [last_left]
        for i, branch in enumerate(mbs.branches):
            dist = mbs.distances[mbs.end_indices[i]]
            tickd.append(dist)
            left, right = branch
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
                plugins.ImportDataset1D('nbands', [mbs.nbands]),
                plugins.ImportDataset1D('distances1', mbs.distances)
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
            datasets.append(plugins.ImportDataset1D(name, dos.densities[spin]-efermi))

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