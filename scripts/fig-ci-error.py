import svgutils.transform as sg
from common import load_svg, label_plot


fig = sg.SVGFigure("4.15in", "1.8in")
a = load_svg(snakemake.input.a)
b = load_svg(snakemake.input.b)
c = load_svg(snakemake.input.c)
d = load_svg(snakemake.input.d)
b.moveto(100, 0)
c.moveto(200, 0)
d.moveto(300, 0)

la = label_plot(5, 10, "a")
lb = label_plot(95, 10, "b")
lc = label_plot(195, 10, "c")
ld = label_plot(295, 10, "d")

fig.append([a, b, c, d, la, lb, lc, ld])
fig.save(snakemake.output[0])
