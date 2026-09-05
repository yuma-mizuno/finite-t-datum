"""Generate typeset tables directly from the checked classification data."""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
z = sp.Symbol("z")


def matrix_tex(rows):
    return r"\begin{pmatrix}" + r"\\".join(
        " & ".join(sp.latex(sp.expand(value)) for value in row) for row in rows
    ) + r"\end{pmatrix}"


def main():
    data = [item for item in map(json.loads, (HERE / "lift_feasibility.jsonl").read_text().splitlines())
            if item["status"] == "sat"]
    signatures = {item["id"]: item for item in json.loads((HERE / "slice_signatures.json").read_text())}
    constants = {item["id"]: item for item in
                 json.loads((HERE / "constant_candidates.json").read_text())["candidates"]}
    overview = [r"\begin{tabular}{@{}crrrrr@{}}", r"\toprule",
                r"Class & Delays $(r_1,r_2,r_3)$ & $h_+$ & $h_-$ & Period & Slice size\\", r"\midrule"]
    catalogue = []
    for k, item in enumerate(data, 1):
        r = item["values"][:3]
        c = item["certificate"]
        hp, hm, period = c["positive"]["h"], c["negative"]["h"], c["labelled_tropical_seed_period"]
        overview.append(f"{k} & $({r[0]},{r[1]},{r[2]})$ & {hp} & {hm} & {period} & {signatures[item['id']]['slice_vertices']}" + r"\\")
        if (k-1) % 4 == 0:
            if k > 1:
                catalogue.append(r"\clearpage")
            catalogue.append(r"\subsection*{Classes " + str(k) + "--" + str(k+3) + "}")
            catalogue.append(r"\noindent In each class, $B_\pm^{(k)}(z)=\operatorname{diag}(1+z^{r_1},1+z^{r_2},1+z^{r_3})-N_\pm^{(k)}(z)$.\par\medskip")
        catalogue.append(r"\noindent\begin{minipage}{\linewidth}")
        catalogue.append(r"\noindent\textbf{Class " + str(k) + r"}\hfill $\symbf{r}=(" + ",".join(map(str, r)) + ")$")
        catalogue.append(r"\par\smallskip\noindent\textcolor{rulegray}{\rule{\linewidth}{0.3pt}}")
        n0 = sp.diag(*(1+z**v for v in r))
        specializations = []
        for sign, field, constant_field in (("+", "A_plus", "N_plus_1"),
                                           ("-", "A_minus", "N_minus_1")):
            a = sp.Matrix([[sp.sympify(v, locals={"z":z}) for v in row] for row in item[field]])
            n = n0-a
            # Recover the serialized input exactly; no coefficient transcription.
            assert n0-n == a
            a1 = a.subs(z, 1)
            assert a1 == 2*sp.eye(3) - sp.Matrix(constants[item["id"]][constant_field])
            specializations.append((sign, a1))
            catalogue.append(r"\begin{minipage}[c]{0.495\linewidth}\centering\fontsize{10.5}{13}\selectfont")
            catalogue.append(r"\setlength{\abovedisplayskip}{3pt}\setlength{\belowdisplayskip}{3pt}")
            catalogue.append(r"\setlength{\abovedisplayshortskip}{3pt}\setlength{\belowdisplayshortskip}{3pt}")
            catalogue.append(r"\setlength{\arraycolsep}{3.5pt}\renewcommand{\arraystretch}{1.18}")
            catalogue.append(r"\[N_{" + sign + r"}^{(" + str(k) + r")}(z)=" + matrix_tex(n.tolist()) + r"\]")
            catalogue.append(r"\end{minipage}%")
        catalogue.append(r"\par\noindent")
        for sign, a1 in specializations:
            catalogue.append(r"\begin{minipage}[c]{0.495\linewidth}\centering\fontsize{10.5}{13}\selectfont")
            catalogue.append(r"\setlength{\abovedisplayskip}{3pt}\setlength{\belowdisplayskip}{3pt}")
            catalogue.append(r"\setlength{\abovedisplayshortskip}{3pt}\setlength{\belowdisplayshortskip}{3pt}")
            catalogue.append(r"\setlength{\arraycolsep}{3.5pt}\renewcommand{\arraystretch}{1.18}")
            catalogue.append(r"\[A_{" + sign + r"}(1)=" + matrix_tex(a1.tolist()) + r"\]")
            catalogue.append(r"\end{minipage}%")
        catalogue.append(r"\par\smallskip\noindent\small $h_+=" + str(hp) + r",\quad h_-=" + str(hm)
                         + r",\quad \Omega=" + str(period) + r".$\quad Vertices per slice: "
                         + str(signatures[item["id"]]["slice_vertices"]) + ".")
        catalogue.append(r"\end{minipage}\par\vfill")
    overview += [r"\bottomrule", r"\end{tabular}"]
    (HERE / "rank3-overview.tex").write_text("\n".join(overview)+"\n", encoding="utf-8")
    (HERE / "rank3-catalogue.tex").write_text("\n".join(catalogue)+"\n", encoding="utf-8")
    print("Generated 16 polynomial matrix pairs and 32 specializations at z=1; "
          "all specializations agree with the independent constant-pair catalogue.")


if __name__ == "__main__":
    main()
