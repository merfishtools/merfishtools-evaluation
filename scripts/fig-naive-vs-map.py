import svgutils.transform as sg
from common import load_svg, label_plot


fig = sg.SVGFigure("4.15in", "1.8in")
a = load_svg(snakemake.input.a)
b = load_svg(snakemake.input.b)
b.moveto(160, 0)

la = label_plot(5, 10, "a")
lb = label_plot(155, 10, "b")

fig.append([a, b, la, lb])
fig.save(snakemake.output[0])
