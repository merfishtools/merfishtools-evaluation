import svgutils.transform as sg
from common import load_svg, label_plot

fig = sg.SVGFigure("6.4in", "1.7in")
a = load_svg(snakemake.input[0])
b = load_svg(snakemake.input[1])
#c = load_svg(snakemake.input.c)
#d = load_svg(snakemake.input.d)
e = load_svg(snakemake.input[2])
#f = load_svg(snakemake.input.f)
b.moveto(190, 0)
#c.moveto(380, 0)
#d.moveto(0, 160)
e.moveto(380, 0)
#f.moveto(380, 160)

la = label_plot(5,10, "a")
lb = label_plot(195,10, "b")
lc = label_plot(385,10, "c")

fig.append([a, b, e, la, lb, lc])
fig.save(snakemake.output[0])
