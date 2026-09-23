"""Check panel keys at export size without changing plotted values."""
import numpy as np
from matplotlib.collections import PathCollection, LineCollection
from matplotlib.path import Path
from matplotlib.transforms import Bbox


def inspect_figure(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    panels = []
    for number, ax in enumerate(fig.axes, 1):
        frame = ax.get_window_extent(renderer)
        # Axis labels (especially colorbar titles) must not enter another panel.
        label_checks = []
        for label in [ax.xaxis.label,ax.yaxis.label]:
            if not label.get_text():
                continue
            box = label.get_window_extent(renderer)
            assert fig.bbox.contains(box.x0,box.y0) and fig.bbox.contains(box.x1,box.y1), (number,label.get_text(),'label clipped')
            for other in fig.axes:
                if other is not ax:
                    assert not box.overlaps(other.get_window_extent(renderer)), (number,label.get_text(),'label overlaps another axes')
            label_checks.append({'text':label.get_text(),'clipped':False,'overlaps_other_axes':False})
        keys = []
        legend = ax.get_legend()
        if legend is not None:
            keys.append(('legend', legend.get_window_extent(renderer),
                         [t.get_text() for t in legend.get_texts()]))
        for text in ax.texts:
            if text.get_gid() == 'panel_note':
                keys.append(('note', text.get_window_extent(renderer), [text.get_text()]))
        checks = []
        for kind, box, labels in keys:
            contained = bool(frame.contains(box.x0, box.y0) and frame.contains(box.x1, box.y1))
            collisions = []
            for i, collection in enumerate(ax.collections):
                if isinstance(collection, PathCollection):
                    positions = collection.get_offset_transform().transform(collection.get_offsets())
                    sizes = collection.get_sizes()
                    radii = np.sqrt(sizes) * fig.dpi / 72 / 2 if len(sizes) else [0]
                    for j, (x, y) in enumerate(positions):
                        radius = radii[j % len(radii)]
                        if box.overlaps(Bbox.from_extents(x-radius,y-radius,x+radius,y+radius)):
                            collisions.append(f'points {i}:{j}')
                elif isinstance(collection, LineCollection):
                    for j, segment in enumerate(collection.get_segments()):
                        if len(segment) and Path(collection.get_transform().transform(segment)).intersects_bbox(box, filled=False):
                            collisions.append(f'interval {i}:{j}')
            for i, line in enumerate(ax.lines):
                # Reference baselines use blended coordinates; keys have a white
                # background in reserved space above every data observation.
                if line.get_transform() != ax.transData:
                    continue
                points = line.get_xydata()
                if len(points) == 0:
                    continue
                display = line.get_transform().transform(points)
                if line.get_linestyle() not in ('None', 'none', '', ' ') and len(display) > 1:
                    if Path(display).intersects_bbox(box, filled=False):
                        collisions.append(f'data line {i}')
                if line.get_marker() not in ('None', 'none', '', ' '):
                    radius = line.get_markersize() * fig.dpi / 72 / 2
                    if any(box.overlaps(Bbox.from_extents(x-radius,y-radius,x+radius,y+radius)) for x,y in display):
                        collisions.append(f'markers {i}')
            for i, patch in enumerate(ax.patches):
                if box.overlaps(patch.get_window_extent(renderer)):
                    collisions.append(f'bar {i}')
            for i, im in enumerate(ax.images):
                x0,x1,y0,y1 = im.get_extent()
                corners = ax.transData.transform([[min(x0,x1),min(y0,y1)],[max(x0,x1),max(y0,y1)]])
                image_box = Bbox.from_extents(*np.min(corners,axis=0),*np.max(corners,axis=0))
                if box.overlaps(image_box):
                    collisions.append(f'heatmap {i}')
            assert contained, (number, labels, 'key outside panel frame')
            assert not collisions, (number, labels, collisions)
            checks.append({'kind':kind,'labels':labels,'inside_frame':contained,'data_overlap':False})
        assert not any(line.get_visible() for line in ax.get_xgridlines()+ax.get_ygridlines())
        panels.append({'axes':number,'keys':checks,'axis_labels':label_checks,'width_mm':frame.width/fig.dpi*25.4,
                       'height_mm':frame.height/fig.dpi*25.4,
                       'panel_letters':[t.get_text() for t in ax.texts if t.get_text() in list('ABCDEFG')]})
    return panels


def panel_note(ax, text):
    return ax.text(.035,.965,text,transform=ax.transAxes,ha='left',va='top',
                   fontsize=8,gid='panel_note',bbox={'facecolor':'white','edgecolor':'none','pad':1.5})
