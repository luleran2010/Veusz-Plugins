import veusz.document.commandinterface as commandinterface
import veusz.embed as embed
import veusz.plugins as plugins

class PlotPhononBandsPlugin(plugins.ToolsPlugin):
    menu = ('Phonopy', 'Plot phonon bands')
    name = 'Plot phonon bands'
    description_short = 'Plot phonon bands in bands.hdf5'
    description_full = 'Plot phonon band structure in bands.hdf5/bands.h5'

    def __init__(self):
        self.fields = [
            plugins.FieldWidget('widget', descr='Draw on widget', default='', widgettypes='graph'),
            plugins.FieldText('prefix', descr='Prefix'),
            plugins.FieldText('suffix', descr='Suffix')
        ]

    def apply(self, interface: commandinterface.CommandInterface, fields: dict):
        widget = interface.Root.fromPath(fields['widget'])
        prefix = fields['prefix']
        suffix = fields['suffix']

        distances = prefix + 'distances' + suffix
        bands = prefix+'frequencies'+suffix
        tickd, tickl = prefix+'tickd'+suffix, prefix+'tickl'+suffix

        xy = self.draw_bands(widget, distances, bands)

        ticks = widget.Add('xy', name='ticks', xData=tickd, yData=tickd, hide=True, labels=tickl)

        x = widget.x
        x.autoRange.val = 'exact'
        x.mode.val = 'labels'
        x.MajorTicks.manualTicks.val = interface.GetData(tickd)[0].tolist()
        x.MinorTicks.hide.val = True
        x.GridLines.style.val = 'dotted'
        x.GridLines.hide.val = False

        y = widget.y
        y.label.val = 'Frequency/THz'

    def draw_bands(self, widget: embed.WidgetNode, distances: str, bands: str):
        xy = widget.Add('xy', name=bands, marker='none', xData=distances, yData=bands)
        return xy
    
plugins.toolspluginregistry += [
    PlotPhononBandsPlugin
]