import veusz.document.commandinterface as commandinterface
import veusz.embed as embed
import veusz.plugins as plugins

class PlotBandStructurePlugin(plugins.ToolsPlugin):
    menu = ('VASP', 'Plot band structure')
    name = 'Plot band structure'
    description_short = 'Plot band structure in vasprun.xml'
    description_full = 'Plot band structure in vasprun.xml'

    def __init__(self):
        self.fields = [
            plugins.FieldWidget('widget', descr='Draw on widget', default='', widgettypes='graph'),
            plugins.FieldBool('spin', descr='Spin polarized'),
            plugins.FieldText('prefix', descr='Prefix'),
            plugins.FieldText('suffix', descr='Suffix'),
        ]

    def apply(self, interface: commandinterface.CommandInterface, fields: dict):
        widget = interface.Root.fromPath(fields['widget'])
        spin = fields['spin']
        prefix = fields['prefix']
        suffix = fields['suffix']

        distances = prefix + 'distances' + suffix
        bands = [prefix+'energies_up'+suffix, prefix+'energies_dw'+suffix]
        tickd, tickl = prefix+'tickd'+suffix, prefix+'tickl'+suffix

        bands_up = self.draw_band(widget, distances, bands[0])
        bands_dw = None
        if spin == True:
            bands_dw = self.draw_band(widget, distances, bands[1])
            bands_dw.PlotLine.style.val = 'dashed'

        ticks = widget.Add('xy', name='ticks')
        ticks.xData.val = tickd
        ticks.yData.val = tickd
        ticks.hide.val = True
        ticks.labels.val = tickl

        x = widget.x
        x.autoRange.val = 'exact'
        x.mode.val = 'labels'
        x.MajorTicks.manualTicks.val = interface.GetData(tickd)[0].tolist()
        x.MinorTicks.hide.val = True
        x.GridLines.style.val = 'dotted'
        x.GridLines.hide.val = False

        y = widget.y
        y.label.val = 'Energy/eV'

    def draw_band(self, widget: embed.WidgetNode, distances: str, bands: str):
        xy = widget.Add('xy', name=bands)
        xy.marker.val = 'none'
        xy.xData.val = distances
        xy.yData.val = bands
        return xy

class PlotDOSPlugin(plugins.ToolsPlugin):
    menu = ('VASP', 'Plot DOS')
    name = 'Plot DOS'
    description_short = 'Plot DOS in vasprun.xml'
    description_full = 'Plot band structure in vasprun.xml'

    def __init__(self):
        self.fields = [
            plugins.FieldWidget('widget', descr='Draw on widget', default='', widgettypes='graph'),
            plugins.FieldBool('spin', descr='Spin polarized'),
            plugins.FieldText('prefix', descr='Prefix'),
            plugins.FieldText('suffix', descr='Suffix'),
        ]

    def apply(self, interface: commandinterface.CommandInterface, fields: dict):
        widget = interface.Root.fromPath(fields['widget'])
        spin = fields['spin']
        prefix = fields['prefix']
        suffix = fields['suffix']
        
        energies = prefix + 'energies' + suffix
        densities = [prefix+'dos_up'+suffix, prefix+'dos_dw'+suffix]

        densities_up = self.draw_densities(widget, energies, densities[0])
        if spin == True:
            densities_dw = self.draw_densities(widget, energies, '-'+densities[1])
        
        x = widget.x
        x.MajorTicks.hide.val = True
        x.MinorTicks.hide.val = True
        x.TickLabels.hide.val = True

        y = widget.y
        y.autoRange.val = 'exact'
        y.label.val = 'Energy/eV'

    def draw_densities(self, widget: embed.WidgetNode, energies: str, densities: str):
        xy = widget.Add('xy', name=densities)
        xy.marker.val = 'none'
        xy.xData.val = densities
        xy.yData.val = energies
        return xy

plugins.toolspluginregistry += [
    PlotBandStructurePlugin,
    PlotDOSPlugin
]