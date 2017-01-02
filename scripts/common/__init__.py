import svgutils.transform as sg


def load_svg(path):
    return sg.fromfile(path).getroot()


def label_plot(x, y, label):
    return sg.TextElement(x, y, label, size=12, weight="bold")
