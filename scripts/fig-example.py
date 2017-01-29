import svgutils.transform as sg
from common import load_svg, label_plot

fig = sg.SVGFigure("6.4in", "3.6in")
a = load_svg(snakemake.input.a)
b = load_svg(snakemake.input.b)
c = load_svg(snakemake.input.c)
d = load_svg(snakemake.input.d)
e = load_svg(snakemake.input.e)
f = load_svg(snakemake.input.f)
b.moveto(190, 0)
c.moveto(380, 0)
d.moveto(0, 160)
e.moveto(190, 160)
f.moveto(380, 160)

la = label_plot(5,10, "a")
lb = label_plot(195,10, "b")
lc = label_plot(385,10, "c")
ld = label_plot(5,170, "d")
le = label_plot(195,170, "e")
lf = label_plot(385,170, "f")

fig.append([a, b, c, d, e, f, la, lb, lc, ld, le, lf])
fig.save(snakemake.output[0])
