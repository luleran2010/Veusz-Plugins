import veusz.document.commandinterface as commandinterface
import veusz.embed as embed
import veusz.plugins as plugins

def draw_bs_kpath(interface: commandinterface.CommandInterface, graph: embed.WidgetNode, tickd: str, tickl: str):
    ticks = graph.Add('xy', name='ticks')
    ticks.xData.val = tickd
    ticks.yData.val = tickd
    ticks.hide.val = True
    ticks.labels.val = tickl

    x = graph.x
    x.autoRange.val = 'exact'
    x.mode.val = 'labels'
    x.MajorTicks.manualTicks.val = interface.GetData(tickd)[0].tolist()
    x.MinorTicks.hide.val = True
    x.GridLines.style.val = 'dotted'
    x.GridLines.hide.val = False

def draw_bands(graph: embed.WidgetNode, distances: str, bands: str):
    xy = graph.Add('xy', name=bands, marker='none', xData=distances, yData=bands)
    return xy

class PlotBandStructurePlugin(plugins.ToolsPlugin):
    menu = ('VASP', 'Plot band structure')
    name = 'Plot band structure'
    description_short = 'Plot band structure in vasprun.xml'
    description_full = 'Plot band structure in vasprun.xml'

    def __init__(self):
        self.fields = [
            plugins.FieldWidget('graph', descr='Draw on graph', default='', widgettypes='graph'),
            plugins.FieldBool('spin', descr='Spin polarized'),
            plugins.FieldText('prefix', descr='Prefix'),
            plugins.FieldText('suffix', descr='Suffix'),
        ]

    def apply(self, interface: commandinterface.CommandInterface, fields: dict):
        graph = interface.Root.fromPath(fields['graph'])
        spin = fields['spin']
        prefix = fields['prefix']
        suffix = fields['suffix']

        distances = prefix + 'distances' + suffix
        bands = [prefix+'bands_up'+suffix, prefix+'bands_dw'+suffix]
        tickd, tickl = prefix+'tickd'+suffix, prefix+'tickl'+suffix

        bands_up = draw_bands(graph, distances, bands[0])
        bands_dw = None
        if spin == True:
            bands_dw = draw_bands(graph, distances, bands[1])
            bands_dw.PlotLine.style.val = 'dashed'

        draw_bs_kpath(interface, graph, tickd, tickl)
        y = graph.y
        y.label.val = 'Energy/eV'

def hide_dos_x(graph: embed.WidgetNode):
    x = graph.x
    x.MajorTicks.hide.val = True
    x.MinorTicks.hide.val = True
    x.TickLabels.hide.val = True

def draw_densities(widget: embed.WidgetNode, energies: str, densities: str, down: bool=False):
    xData = densities
    if down:
        xData = '-'+densities
    xy = widget.Add('xy', name=densities, marker='none', xData=xData, yData=energies)
    return xy

class PlotDOSPlugin(plugins.ToolsPlugin):
    menu = ('VASP', 'Plot DOS')
    name = 'Plot DOS'
    description_short = 'Plot DOS in vasprun.xml'
    description_full = 'Plot DOS in vasprun.xml'

    def __init__(self):
        self.fields = [
            plugins.FieldWidget('graph', descr='Draw on graph', default='', widgettypes='graph'),
            plugins.FieldBool('spin', descr='Spin polarized'),
            plugins.FieldText('orbitals', descr='Orbitals (separated by ,)', default=''),
            plugins.FieldText('prefix', descr='Prefix'),
            plugins.FieldText('suffix', descr='Suffix'),
        ]

    def apply(self, interface: commandinterface.CommandInterface, fields: dict):
        graph = interface.Root.fromPath(fields['graph'])
        spin = fields['spin']
        orbitals = [orbital.strip().replace(' ', '_') for orbital in fields['orbitals'].split(',')]
        prefix = fields['prefix']
        suffix = fields['suffix']

        energies = prefix + 'energies' + suffix
        # densities = [prefix+'tdos_up'+suffix, prefix+'tdos_dw'+suffix]
        orbitals = [orbital.replace(' ', '_') for orbital in orbitals]
        densities = []
        for orbital in orbitals:
            if orbital == '':
                densities = [prefix+'tdos_up'+suffix, prefix+'tdos_dw'+suffix]
            else:
                densities = [prefix+'pdos_'+orbital+'_up', prefix+'pdos_'+orbital+'_dw']
            draw_densities(graph, energies, densities[0])
            if spin == True:
                draw_densities(graph, energies, densities[1], True)
        
        hide_dos_x(graph)
        y = graph.y
        y.autoRange.val = 'exact'
        y.label.val = 'Energy/eV'

class PlotBSDOSPlugin(plugins.ToolsPlugin):
    menu = ('VASP', 'Plot band structure and DOS')
    name = 'Plot BSDOS'
    description_short = 'Plot band structure and DOS'
    description_full = 'Plot band structure along with DOS'

    def __init__(self):
        self.fields = [
            plugins.FieldWidget('widget', descr='Draw on widget', default='', widgettypes='page'),
            plugins.FieldBool('spin', descr='Spin polarized'),
            plugins.FieldFloatList('colscale', descr='Column scale [BS,DOS]', default=[3.0, 1.0]),
            plugins.FieldText('prefix', descr='Prefix (separated by ,)', default=','),
            plugins.FieldText('suffix', descr='Suffix (separated by ,)'),
        ]

    def apply(self, interface: commandinterface.CommandInterface, fields: dict):
        widget = interface.Root.fromPath(fields['widget'])
        spin = fields['spin']
        colscale = fields['colscale']
        prefix = fields['prefix'].split(',')
        suffix = fields['suffix'].split(',')

        distances = prefix[0] + 'distances' + suffix[0]
        bands = [prefix[0]+'bands_up'+suffix[0], prefix[0]+'bands_dw'+suffix[0]]
        tickd, tickl = prefix[0]+'tickd'+suffix[0], prefix[0]+'tickl'+suffix[0]

        energies = prefix[-1] + 'energies' + suffix[-1]
        densities = [prefix[-1]+'tdos_up'+suffix[-1], prefix[-1]+'tdos_dw'+suffix[-1]]

        grid = widget.Add('grid', name='bsdos')
        grid.scaleCols.val = colscale
        grid.internalMargin.val = '0.2cm'
        grid.Add('axis', name='y', label='Energy/eV', direction='vertical')

        bs = grid.Add('graph', name='bs', autoadd=False,
                      leftMargin='0cm', rightMargin='0cm', topMargin='0cm', bottomMargin='0cm')
        bs.Add('axis', name='x', direction='horizontal')
        bands_up = draw_bands(bs, distances, bands[0])
        bands_dw = None
        if spin == True:
            bands_dw = draw_bands(bs, distances, bands[1])
            bands_dw.PlotLine.style.val = 'dashed'

        draw_bs_kpath(interface, bs, tickd, tickl)

        dos = grid.Add('graph', name='dos', autoadd=False,
                       leftMargin='0cm', rightMargin='0cm', topMargin='0cm', bottomMargin='0cm')
        dos.Add('axis', name='x', direction='horizontal')
        densities_up = draw_densities(dos, energies, densities[0])
        if spin == True:
            densities_dw = draw_densities(dos, energies, densities[1], True)
        
        hide_dos_x(dos)

plugins.toolspluginregistry += [
    PlotBandStructurePlugin,
    PlotDOSPlugin,
    PlotBSDOSPlugin
]