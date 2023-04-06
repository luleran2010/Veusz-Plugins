import numpy as np
import h5py
from veusz.plugins import *


class ImportPluginPhononDispersion(ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "Phonon Dispersion plugin"
    author = "Leran Lu"
    description = "Reads the phonon dispersion from the Phonopy-generated HDF5 file"

    # Uncomment this line for the plugin to get its own tab
    #promote_tab='Example'

    file_extensions = set(['.h5', '.hdf5'])

    def __init__(self):
        ImportPlugin.__init__(self)
        self.fields = [
            ImportFieldCheck('details', descr='Detailed Information')
        ]

    def doImport(self, params: ImportPluginParams):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """

        phf = h5py.File(params.filename)
        distance = phf['distance'][:]
        frequency = phf['frequency'][:]
        label = phf['label'][:]
        nqpoint = phf['nqpoint'][0]
        path = phf['path'][:]
        segment_nqpoint = phf['segment_nqpoint'][:]

        nbands = frequency.shape[-1]
        dist = np.concatenate((distance, np.array([np.nan]*3).reshape((3,1))), axis=1).flatten()
        dist = dist.reshape((1,-1))
        dist = np.repeat(dist, nbands, axis=0)
        dist = np.concatenate((dist, np.array([np.nan]*nbands).reshape((nbands,1))), axis=1).flatten()

        bk = np.array([np.nan]*nbands).reshape((1,nbands))
        freq = np.concatenate((frequency[0,:,:], bk, frequency[1,:,:], bk, frequency[2,:,:], bk), axis=0).T
        freq = np.concatenate((freq, np.array([np.nan]*(nbands)).reshape((nbands,1))), axis=1).flatten()

        tickd = np.append(distance[:,0], distance[-1,-1])
        tickl = [bytes.decode(i) for i in np.append(label[:,0], label[-1,-1])]
        for i in np.arange(0, len(tickl)):
            if tickl[i] == 'Gamma':
                tickl[i] = '\Gamma'

        details = []
        if params.field_results['details']:
            details = [
                ImportDataset1D('nbands', [nbands]),
                ImportDataset1D('distances1', distance)
            ]

        return [
            ImportDataset1D('distances', dist),
            ImportDataset1D('frequencies', freq),
            ImportDataset1D('tickd', tickd),
            ImportDatasetText('tickl', tickl)
        ] +  details

# add the class to the registry. An instance also works, but is deprecated
importpluginregistry.append(ImportPluginPhononDispersion)