import svgutils.transform as sg
from common import load_svg, label_plot

fig = sg.SVGFigure("7.7in", "3.6in")
a = load_svg(snakemake.input[0])
b = load_svg(snakemake.input[1])
c = load_svg(snakemake.input[2])
d = load_svg(snakemake.input[3])
a.moveto(10, 10, scale=0.6)
b.moveto(20, 130, scale=0.36)
c.moveto(40, 190, scale=0.6)
d.moveto(230, 10, scale=0.8)


la = label_plot(5, 10, "a")
lb = label_plot(5, 130, "b")
lc = label_plot(5, 190, "c")
ld = label_plot(215, 10, "d")

fig.append([a, b, c, d, la, lb, lc, ld])
fig.save(snakemake.output[0])
