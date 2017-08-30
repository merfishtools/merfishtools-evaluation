import svgutils.transform as sg
from common import load_svg, label_plot


fig = sg.SVGFigure("5.6in", "1.8in")
a = load_svg(snakemake.input.a)
b = load_svg(snakemake.input.b)
a.moveto(10, 0)
b.moveto(260, 0)

la = label_plot(5, 10, "a")
lb = label_plot(255, 10, "b")

fig.append([a, b, la, lb])
fig.save(snakemake.output[0])
