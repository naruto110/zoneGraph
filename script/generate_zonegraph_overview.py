#!/usr/bin/env python3
"""Generate ZoneGraph architecture overview as vector PDF and SVG."""

from pathlib import Path


OUT_DIR = Path("/home/femu/project/zoneGraph/figures")
PDF_PATH = OUT_DIR / "zonegraph_overview.pdf"
SVG_PATH = OUT_DIR / "zonegraph_overview.svg"

W, H = 792, 490


def esc_pdf(text):
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def esc_xml(text):
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


class PDFCanvas:
    def __init__(self, width, height):
        self.w = width
        self.h = height
        self.ops = []

    def color(self, stroke=None, fill=None):
        if stroke:
            r, g, b = stroke
            self.ops.append(f"{r:.3f} {g:.3f} {b:.3f} RG")
        if fill:
            r, g, b = fill
            self.ops.append(f"{r:.3f} {g:.3f} {b:.3f} rg")

    def line_width(self, width):
        self.ops.append(f"{width:.2f} w")

    def dash(self, pattern=None):
        if pattern:
            self.ops.append(f"[{' '.join(str(x) for x in pattern)}] 0 d")
        else:
            self.ops.append("[] 0 d")

    def rect(self, x, y, w, h, fill, stroke=(0.22, 0.22, 0.22), lw=0.8):
        self.color(stroke=stroke, fill=fill)
        self.line_width(lw)
        self.ops.append(f"{x:.2f} {y:.2f} {w:.2f} {h:.2f} re B")

    def line(self, x1, y1, x2, y2, stroke=(0.2, 0.2, 0.2), lw=0.8, dashed=False):
        self.color(stroke=stroke)
        self.line_width(lw)
        self.dash([3, 3] if dashed else None)
        self.ops.append(f"{x1:.2f} {y1:.2f} m {x2:.2f} {y2:.2f} l S")
        if dashed:
            self.dash(None)

    def arrow(self, x1, y1, x2, y2, stroke=(0.2, 0.2, 0.2), lw=0.9, dashed=False):
        self.line(x1, y1, x2, y2, stroke=stroke, lw=lw, dashed=dashed)
        import math

        ang = math.atan2(y2 - y1, x2 - x1)
        size = 6.0
        a1 = ang + math.pi * 0.84
        a2 = ang - math.pi * 0.84
        p1 = (x2 + size * math.cos(a1), y2 + size * math.sin(a1))
        p2 = (x2 + size * math.cos(a2), y2 + size * math.sin(a2))
        self.color(stroke=stroke, fill=stroke)
        self.ops.append(
            f"{x2:.2f} {y2:.2f} m {p1[0]:.2f} {p1[1]:.2f} l {p2[0]:.2f} {p2[1]:.2f} l h f"
        )

    def text(self, x, y, text, size=9, fill=(0.05, 0.05, 0.05), bold=False, align="left"):
        font = "F2" if bold else "F1"
        approx = len(text) * size * 0.47
        if align == "center":
            x = x - approx / 2
        elif align == "right":
            x = x - approx
        self.color(fill=fill)
        self.ops.append(f"BT /{font} {size:.1f} Tf {x:.2f} {y:.2f} Td ({esc_pdf(text)}) Tj ET")

    def wrap_text(self, x, y, lines, size=8, fill=(0.05, 0.05, 0.05), bold_first=False, leading=11):
        for i, line in enumerate(lines):
            self.text(x, y - i * leading, line, size=size, fill=fill, bold=(bold_first and i == 0))

    def save(self, path):
        stream = "\n".join(self.ops).encode("latin-1")
        objects = []
        objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
        objects.append(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
        page = (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {self.w} {self.h}] "
            f"/Resources << /Font << /F1 5 0 R /F2 6 0 R >> >> /Contents 4 0 R >>"
        ).encode("latin-1")
        objects.append(page)
        objects.append(f"<< /Length {len(stream)} >>\nstream\n".encode("latin-1") + stream + b"\nendstream")
        objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
        objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>")

        out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
        offsets = [0]
        for i, obj in enumerate(objects, start=1):
            offsets.append(len(out))
            out += f"{i} 0 obj\n".encode("latin-1") + obj + b"\nendobj\n"
        xref = len(out)
        out += f"xref\n0 {len(objects)+1}\n0000000000 65535 f \n".encode("latin-1")
        for off in offsets[1:]:
            out += f"{off:010d} 00000 n \n".encode("latin-1")
        out += (
            f"trailer << /Size {len(objects)+1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n"
        ).encode("latin-1")
        path.write_bytes(out)


class SVGCanvas:
    def __init__(self, width, height):
        self.w = width
        self.h = height
        self.items = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#333"/></marker></defs>',
            '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
        ]

    def y(self, y):
        return self.h - y

    def rect(self, x, y, w, h, fill, stroke="#444", lw=0.8):
        self.items.append(
            f'<rect x="{x}" y="{self.y(y+h)}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="{lw}"/>'
        )

    def line(self, x1, y1, x2, y2, stroke="#333", lw=0.8, dashed=False, arrow=False):
        dash = ' stroke-dasharray="3 3"' if dashed else ""
        marker = ' marker-end="url(#arrow)"' if arrow else ""
        self.items.append(
            f'<line x1="{x1}" y1="{self.y(y1)}" x2="{x2}" y2="{self.y(y2)}" stroke="{stroke}" stroke-width="{lw}"{dash}{marker}/>'
        )

    def text(self, x, y, text, size=9, fill="#111", bold=False, align="start"):
        weight = "700" if bold else "400"
        anchor = {"left": "start", "center": "middle", "right": "end"}.get(align, align)
        self.items.append(
            f'<text x="{x}" y="{self.y(y)}" font-family="Helvetica, Arial, sans-serif" font-size="{size}" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{esc_xml(text)}</text>'
        )

    def save(self, path):
        path.write_text("\n".join(self.items + ["</svg>\n"]))


def draw(canvas, svg=False):
    if isinstance(canvas, PDFCanvas):
        rgb = {
            "blue1": (0.90, 0.95, 0.99),
            "blue2": (0.78, 0.88, 0.96),
            "green1": (0.90, 0.96, 0.90),
            "green2": (0.78, 0.90, 0.78),
            "orange1": (0.99, 0.93, 0.84),
            "orange2": (0.94, 0.70, 0.38),
            "gray1": (0.96, 0.96, 0.96),
            "gray2": (0.80, 0.80, 0.80),
            "dark": (0.18, 0.18, 0.18),
        }
        fills = rgb
        stroke = (0.25, 0.25, 0.25)
    else:
        fills = {
            "blue1": "#E6F1FB",
            "blue2": "#C8E0F5",
            "green1": "#E7F3E7",
            "green2": "#C6E2C6",
            "orange1": "#FBECD8",
            "orange2": "#D99A4A",
            "gray1": "#F5F5F5",
            "gray2": "#CCCCCC",
            "dark": "#333333",
        }
        stroke = "#444"

    def rect(x, y, w, h, fill, lw=0.8):
        canvas.rect(x, y, w, h, fills[fill], stroke=stroke, lw=lw)

    def text(x, y, s, size=9, bold=False, align="left"):
        canvas.text(x, y, s, size=size, bold=bold, align=align)

    def arrow(x1, y1, x2, y2, dashed=False):
        if isinstance(canvas, PDFCanvas):
            canvas.arrow(x1, y1, x2, y2, stroke=(0.20, 0.20, 0.20), lw=0.9, dashed=dashed)
        else:
            canvas.line(x1, y1, x2, y2, stroke="#333", lw=0.9, dashed=dashed, arrow=True)

    def line(x1, y1, x2, y2, dashed=False):
        if isinstance(canvas, PDFCanvas):
            canvas.line(x1, y1, x2, y2, stroke=(0.45, 0.45, 0.45), lw=0.7, dashed=dashed)
        else:
            canvas.line(x1, y1, x2, y2, stroke="#777", lw=0.7, dashed=dashed)

    # Main regions
    rect(35, 238, 445, 205, "blue1", lw=1.0)
    text(48, 423, "DRAM update and query layer", size=11, bold=True)
    rect(520, 238, 230, 205, "orange1", lw=1.0)
    text(533, 423, "Access and reopen layer", size=11, bold=True)
    rect(35, 45, 715, 150, "green1", lw=1.0)
    text(48, 174, "ZNS persistent layer", size=11, bold=True)

    # Input
    rect(60, 455, 150, 24, "gray1")
    text(135, 463, "Edge update stream", size=9, bold=True, align="center")
    arrow(135, 455, 135, 410)

    # DRAM boxes
    rect(60, 360, 150, 50, "blue2")
    text(135, 392, "Low-degree cache", size=9, bold=True, align="center")
    text(135, 377, "compact adjacency", size=8, align="center")

    rect(260, 360, 165, 50, "blue2")
    text(342, 392, "Pending insert buffer", size=9, bold=True, align="center")
    text(342, 377, "batch merge to PMA", size=8, align="center")

    rect(60, 275, 170, 55, "blue2")
    text(145, 310, "PMA block store", size=9, bold=True, align="center")
    text(145, 294, "update-friendly blocks", size=8, align="center")

    rect(285, 275, 140, 55, "blue2")
    text(355, 310, "Dirty block cache", size=9, bold=True, align="center")
    text(355, 294, "threshold-based flush", size=8, align="center")

    rect(250, 430, 190, 30, "orange1")
    text(345, 441, "Degree-aware promotion", size=9, bold=True, align="center")

    arrow(210, 385, 260, 385)
    arrow(345, 360, 190, 330)
    arrow(230, 302, 285, 302)
    arrow(345, 430, 190, 410)

    # Query/reopen boxes
    rect(545, 365, 180, 45, "orange1")
    text(635, 393, "Vertex-to-block metadata", size=9, bold=True, align="center")
    text(635, 379, "direct neighbor lookup", size=8, align="center")

    rect(545, 300, 180, 42, "orange1")
    text(635, 325, "Query engine", size=9, bold=True, align="center")
    text(635, 311, "DRAM hit or ZNS read", size=8, align="center")

    rect(545, 250, 180, 35, "orange1")
    text(635, 270, "Graph reopen", size=9, bold=True, align="center")
    text(635, 258, "scan records, rebuild metadata", size=7.5, align="center")

    arrow(635, 365, 635, 342)
    arrow(635, 300, 635, 285)
    arrow(545, 388, 425, 330, dashed=True)
    arrow(545, 318, 425, 302, dashed=True)

    # Flush to ZNS
    arrow(355, 275, 355, 195)
    text(386, 222, "sequential", size=8, align="left")
    text(386, 211, "record flush", size=8, align="left")

    # ZNS zones and records
    zone_xs = [70, 230, 390]
    for i, zx in enumerate(zone_xs):
        rect(zx, 70, 115, 78, "green2")
        text(zx + 57, 132, f"Zone {i}", size=8.5, bold=True, align="center")
        for j in range(3):
            rect(zx + 12, 80 + j * 15, 91, 11, "gray1", lw=0.5)
        text(zx + 57, 96, "records", size=7.5, align="center")

    rect(555, 80, 155, 62, "green2")
    text(632, 124, "Block record", size=9, bold=True, align="center")
    text(632, 108, "<vid, block, version>", size=8, align="center")
    text(632, 94, "adjacency payload", size=8, align="center")

    arrow(505, 112, 555, 112)
    arrow(635, 250, 635, 142)
    arrow(635, 142, 635, 365)

    # GC control
    rect(255, 18, 235, 22, "orange1")
    text(372, 26, "ATZ-GC: zone selection with lower valid-block migration", size=8.5, bold=True, align="center")
    line(372, 45, 372, 70, dashed=True)

    # Innovations strip
    rect(35, 202, 715, 24, "gray1")
    text(392, 211, "Key ideas: degree-aware caching  |  batched PMA merge  |  block-granular ZNS records  |  direct block-list lookup", size=8.5, bold=True, align="center")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf = PDFCanvas(W, H)
    draw(pdf)
    pdf.save(PDF_PATH)

    svg = SVGCanvas(W, H)
    draw(svg)
    svg.save(SVG_PATH)
    print(PDF_PATH)
    print(SVG_PATH)


if __name__ == "__main__":
    main()
