import svgutils.transform as sg
from common import load_svg, label_plot


fig = sg.SVGFigure("5.8in", "3.3in")
a = load_svg(snakemake.input.a)
b = load_svg(snakemake.input.b)
c = load_svg(snakemake.input.c)
d = load_svg(snakemake.input.d)
e = load_svg(snakemake.input.e)
f = load_svg(snakemake.input.f)
b.moveto(260, 0)
c.moveto(0, 160)
d.moveto(128, 160)
e.moveto(250, 160)
f.moveto(384, 160)

la = label_plot(5, 10, "a")
lb = label_plot(265, 10, "d")
lc = label_plot(5, 170, "b")
ld = label_plot(130, 170, "c")
le = label_plot(255, 170, "e")
lf = label_plot(385, 170, "f")

fig.append([a, b, c, d, e, f, la, lb, lc, ld, le, lf])
fig.save(snakemake.output[0])
