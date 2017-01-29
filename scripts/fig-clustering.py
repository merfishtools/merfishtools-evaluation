import svgutils.transform as sg
from common import load_svg, label_plot

fx, fy = "5.1in", "1.3in"
if snakemake.wildcards.dataset == "1001genesData":
    fx, fy = "5.5in", "1.3in"
fig = sg.SVGFigure(fx, fy)
a = load_svg(snakemake.input.a)
b = load_svg(snakemake.input.b)
c = load_svg(snakemake.input.c)
d = load_svg(snakemake.input.d)

xb = 170
xc = 275
xd = 380
if snakemake.wildcards.dataset == "1001genesData":
    xb = 170
    xc = 290
    xd = 405

b.moveto(xb, 0)
c.moveto(xc, 0)
d.moveto(xd, 0)

la = label_plot(5,10, "a")
lb = label_plot(xb,10, "b")
lc = label_plot(xc,10, "c")
ld = label_plot(xd,10, "d")

fig.append([a, b, c, d, la, lb, lc, ld])
fig.save(snakemake.output[0])
