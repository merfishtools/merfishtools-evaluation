import svgutils.transform as sg
from common import load_svg, label_plot

fig = sg.SVGFigure("7.5in", "1.9in")
a = load_svg(snakemake.input.mhd4[0])
b = load_svg(snakemake.input.mhd4[1])
c = load_svg(snakemake.input.mhd2[0])
d = load_svg(snakemake.input.mhd2[1])
b.moveto(170, 0)
c.moveto(340, 0)
d.moveto(510, 0)

la = label_plot(5, 10, "a")
lb = label_plot(175, 10, "b")
lc = label_plot(345, 10, "c")
ld = label_plot(515, 10, "d")

fig.append([a, b, c, d, la, lb, lc, ld])
fig.save(snakemake.output[0])
