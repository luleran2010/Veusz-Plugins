import numpy as np
import h5py
import veusz.plugins as plugins

class ImportPluginPhononDispersion(plugins.ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "Phonon Dispersion plugin"
    author = "Leran Lu"
    description = "Reads the phonon dispersion from the Phonopy-generated HDF5 file"

    # Uncomment this line for the plugin to get its own tab
    #promote_tab='Example'

    file_extensions = set(['.h5', '.hdf5'])

    def __init__(self):
        plugins.ImportPlugin.__init__(self)
        self.fields = [
            plugins.ImportFieldCheck('details', descr='Detailed Information')
        ]

    def doImport(self, params: plugins.ImportPluginParams):
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

        nseg = distance.shape[0]
        nbands = frequency.shape[-1]
        dist = np.concatenate((np.array([[np.nan]]*nseg), distance), axis=1).reshape((1,-1))
        dist = np.repeat(dist, nbands, axis=0).flatten()

        freq = np.insert(frequency, range(nseg), np.array([[np.nan]*nbands]), axis=0)
        freq = np.concatenate(
                sum(
                    [[np.array([[np.nan]*nbands]), frequency[i,:,:]] for i in np.arange(nseg)],
                []),
            axis=0)
        freq = freq.T.flatten()

        tickd = np.append(distance[:,0], distance[-1,-1])
        tickl = [bytes.decode(i) for i in np.append(label[:,0], label[-1,-1])]
        for i in np.arange(0, len(tickl)):
            if tickl[i] == '$\Gamma$':
                tickl[i] = '\Gamma'

        details = []
        if params.field_results['details']:
            details = [
                plugins.ImportDataset1D('nbands', [nbands]),
                plugins.ImportDataset1D('distances1', distance)
            ]

        return [
            plugins.ImportDataset1D('distances', dist),
            plugins.ImportDataset1D('frequencies', freq),
            plugins.ImportDataset1D('tickd', tickd),
            plugins.ImportDatasetText('tickl', tickl)
        ] +  details
    
plugins.importpluginregistry += [
    ImportPluginPhononDispersion
]