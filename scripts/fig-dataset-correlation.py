import svgutils.transform as sg
from common import load_svg, label_plot


fig = sg.SVGFigure("4.1in", "1.8in")
a = load_svg(snakemake.input[1])
b = load_svg(snakemake.input[0])
b.moveto(190, 0)

la = label_plot(5, 10, "a")
lb = label_plot(185, 10, "b")

fig.append([a, b, la, lb])
fig.save(snakemake.output[0])
