import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.outputs import BSVasprun, Vasprun
from pymatgen.io.vasp.outputs import Oszicar
import veusz.plugins as plugins

import os

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

        nbands = bs.nb_bands
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

class ImportPluginBandStructureHybrid(plugins.ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "Band Structure plugin (Hybrid)"
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

        dirname = os.path.dirname(params.filename)
        kpoints = Kpoints.from_file(os.path.join(dirname, 'KPOINTS'))
        terms = kpoints.comment.split()
        nkpts = int(terms[6])
        nbranches = int(terms[7])
        npts = [int(term) for term in terms[8:8+nbranches]]
        beg = int(terms[4])
        kpts = np.array([kp.cart_coords for kp in bs.kpoints[beg:]])
        breaks = np.array([np.sum(npts[0:i]) for i in np.arange(len(npts))], dtype=int)
        end_indices = np.array([np.sum(npts[0:i])-1 for i in np.arange(1, len(npts)+1)])

        nbands = bs.nb_bands

        previous_distance = 0.0
        distance = []
        head = 0
        for length in npts:
            diff = np.linalg.norm(np.diff(kpts[head:head+length], axis=0), axis=1)
            dist = np.array([np.sum(diff[0:i]) for i in np.arange(len(diff)+1)]) + previous_distance
            distance.append(dist)
            previous_distance = dist[-1]
            head += length
        distance = np.concatenate(distance, axis=0)

        distances = np.insert(distance, breaks, np.nan, axis=0).reshape((1,-1))
        distances = np.repeat(distances, nbands, axis=0).flatten()
        datasets.append(plugins.ImportDataset1D('distances', distances))

        for spin in bs.bands.keys():
            name = 'up' if spin == Spin.up else 'dw'
            dat = bs.bands[spin][:,beg:]-efermi
            dat = np.insert(dat, breaks, np.nan, axis=1).flatten()
            datasets.append(plugins.ImportDataset1D('bands_'+name, dat))

        kpath = Kpoints.from_file(os.path.join(dirname, 'KPATH.in'))
        convert = lambda label: '\Gamma' if label == 'GAMMA' else label
        branches = [(convert(kpath.labels[2*i]), convert(kpath.labels[2*i+1])) for i in np.arange(len(kpath.labels)/2, dtype=int)]

        dist = distance[0]
        last_left, last_right = branches[0]
        tickd, tickl = [dist], [last_left]
        for i, branch in enumerate(branches):
            dist = distance[end_indices[i]]
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
                plugins.ImportDataset1D('nbands', [nbands]),
                plugins.ImportDataset1D('distances1', distance)
            ]

        return datasets