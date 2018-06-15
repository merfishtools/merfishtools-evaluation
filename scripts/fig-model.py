import svgutils.transform as sg
from common import load_svg, label_plot

fig = sg.SVGFigure("6.5in", "2.1in")
a = load_svg(snakemake.input[0])
b = load_svg(snakemake.input[1])
c = load_svg(snakemake.input[2])
d = load_svg(snakemake.input[3])
a.moveto(10, 10, scale=0.6)
b.moveto(20, 130, scale=0.4)
c.moveto(250, 10, scale=0.65)
d.moveto(440, 10, scale=1.0)


la = label_plot(5, 10, "a")
lb = label_plot(5, 130, "b")
lc = label_plot(240, 10, "c")
ld = label_plot(430, 10, "d")

fig.append([a, b, c, d, la, lb, lc, ld])
fig.save(snakemake.output[0])
