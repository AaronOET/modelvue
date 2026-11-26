import importlib
import json

import pytest

import modelvue
from modelvue import cli, describe
from modelvue import detect

STDLIB_TOOLS = ["dfmvue", "dfmvue23", "epanetvue", "swmmvue", "xbvue"]


@pytest.mark.parametrize("tool", STDLIB_TOOLS)
def test_tool_imports_and_api(tool):
    mod = importlib.import_module(f"modelvue.{tool}")
    for name in ("build_config", "show", "render_html", "write_html", "main"):
        assert callable(getattr(mod, name)), name


def test_flovue_import():
    pytest.importorskip("flopy")
    importlib.import_module("modelvue.flovue")


def test_detect(epanet_inp, swmm_inp, xbeach_case, mdu_2d, mdu_1d, tmp_path):
    assert detect(epanet_inp) == "epanetvue"
    assert detect(swmm_inp) == "swmmvue"
    assert detect(xbeach_case) == "xbvue"
    assert detect(xbeach_case / "params.txt") == "xbvue"
    assert detect(mdu_2d) == "dfmvue23"
    assert detect(mdu_1d) == "dfmvue"
    (tmp_path / "mfsim.nam").write_text("BEGIN models\nEND models\n")
    assert detect(tmp_path / "mfsim.nam") == "flovue"
    with pytest.raises(ValueError, match="more than one kind"):
        detect(tmp_path)
    with pytest.raises(ValueError, match="no supported model"):
        detect(tmp_path / "case1" / "x.grd")


def test_epanet_load_and_report(epanet_inp, tmp_path):
    cfg = modelvue.load(epanet_inp)
    assert cfg["_tool"] == "epanetvue"
    out = modelvue.report(cfg, tmp_path / "net1.html")
    html = out.read_text(encoding="utf-8")
    assert html.startswith("<!DOCTYPE html>") and "Net1" in html
    json.dumps(cfg, default=str)


def test_swmm_load(swmm_inp, capsys):
    cfg = modelvue.load(swmm_inp)
    modelvue.show(cfg, full=True)
    out = capsys.readouterr().out
    assert "C1" in out and "J2" in out


def test_xbeach_load(xbeach_case, tmp_path):
    cfg = modelvue.load(xbeach_case)
    modelvue.report(cfg, tmp_path / "xb.html")
    assert (tmp_path / "xb.html").stat().st_size > 1000


@pytest.mark.parametrize("tool", ["dfmvue", "dfmvue23"])
def test_dfm_without_net_file(tool, mdu_2d, tmp_path):
    # the net file is missing: the reader should warn, not crash
    cfg = modelvue.load(mdu_2d, tool=tool)
    modelvue.report(cfg, tmp_path / f"{tool}.html")


def test_cli_auto_and_explicit(epanet_inp, tmp_path, capsys):
    assert cli.main([str(epanet_inp)]) == 0
    assert "EPANET" in capsys.readouterr().out.upper()
    js = tmp_path / "c.json"
    html = tmp_path / "c.html"
    assert cli.main(["epanetvue", str(epanet_inp), "--json", str(js), "--html", str(html), "--no-browser"]) == 0
    assert json.loads(js.read_text())["models"]
    assert html.exists()


def test_cli_detect_and_errors(epanet_inp, tmp_path, capsys):
    assert cli.main(["detect", str(epanet_inp)]) == 0
    assert capsys.readouterr().out.startswith("epanetvue")
    assert cli.main([str(tmp_path / "nothing_here")]) == 2
    assert cli.main(["--version"]) == 0
    assert modelvue.__version__ in capsys.readouterr().out


def test_cli_help_and_tool_info(capsys):
    assert cli.main([]) == 0
    out = capsys.readouterr().out
    assert all(t in out for t in describe.TOOL_DESCRIPTIONS)
    assert cli.main(["swmmvue"]) == 0
    assert "--- swmmvue ---" in capsys.readouterr().out


def test_tool_main_direct(swmm_inp, tmp_path):
    # the console scripts call <tool>.main() directly, as in ncftools
    from modelvue import swmmvue
    swmmvue.main([str(swmm_inp), "--html", str(tmp_path / "s.html"), "--no-browser"])
    assert (tmp_path / "s.html").exists()


# --------------------------------------------------------------------------- HEC-RAS / HEC-HMS
def test_rasvue(ras_project, tmp_path, capsys):
    from modelvue import rasvue
    assert detect(ras_project) == "rasvue"
    assert detect(ras_project.parent) == "rasvue"          # the ESRI gis.prj is ignored
    cfg = modelvue.load(ras_project)
    p = cfg["projects"][0]
    assert p["units"] == "SI" and p["description"] == "A test project."
    gm = p["geometries"][0]
    assert gm["counts"] == {"River reaches": 1, "Cross sections": 3}
    xs = [n for n in gm["nodes"] if n["type"] == 1]
    assert xs[0]["min_elev"] == 95 and (xs[0]["n_lob"], xs[0]["n_ch"], xs[0]["n_rob"]) == (0.06, 0.035, 0.06)
    assert gm["reaches"][0]["length_xy"] == 190
    u = p["unsteady"][0]
    assert [b["type"] for b in u["boundaries"]] == ["Flow Hydrograph", "Normal depth"]
    assert u["boundaries"][0]["max"] == 40
    pl = p["plans"][0]
    assert pl["kind"] == "unsteady" and pl["duration_h"] == 5 and pl["computed"]
    assert pl["messages"]["finished"]
    modelvue.show(cfg, full=True)
    assert "Tiny plan" in capsys.readouterr().out
    html = modelvue.report(cfg, tmp_path / "ras.html").read_text(encoding="utf-8")
    assert "HEC-RAS" in html and "<svg" in html
    json.dumps(rasvue.strip_private(cfg), default=str)


def test_hmsvue(hms_project, tmp_path, capsys):
    from modelvue import hmsvue
    assert detect(hms_project) == "hmsvue"
    cfg = modelvue.load(hms_project)
    p = cfg["projects"][0]
    b = p["basins"][0]
    assert b["counts"] == {"Subbasin": 2, "Reach": 1, "Junction": 1} and b["total_area"] == 4
    assert any("unknown element 'Nowhere'" in w for w in p["warnings"])
    c = p["controls"][0]
    assert c["duration_h"] == 24 and c["steps"] == 96
    r = p["runs"][0]
    assert r["results"]["elements"][0]["stats"]["Outflow Maximum"] == ("12.5", "M3/S")
    assert r["log_summary"]["counts"] == {"NOTE": 2, "WARNING": 1, "ERROR": 0}
    modelvue.show(cfg)
    assert "12.5 M3/S" in capsys.readouterr().out
    html = modelvue.report(cfg, tmp_path / "hms.html").read_text(encoding="utf-8")
    assert "HEC-HMS" in html
    json.dumps(hmsvue.strip_private(cfg), default=str)


def test_hms_parser_subblocks():
    from modelvue.hmsvue import parse_blocks
    b = parse_blocks("Reservoir: R\n     Route: Controlled Outflow\n     Conduit: Culvert\n     Rise: 4\n     End Conduit:\n"
                     "     Start Date: 1 Jan\n     End Date: 2 Jan\nEnd:\n")[0]
    assert dict(b["items"]) == {"Route": "Controlled Outflow", "Start Date": "1 Jan", "End Date": "2 Jan"}
    assert b["sub"] == [{"type": "Conduit", "name": "Culvert", "items": [("Rise", "4")]}]


def test_ras_pure_2d_geometry(tmp_path):
    # a geometry with no river reach starts directly with a storage area (2D flow area) and BC lines
    from modelvue.rasvue import read_geometry
    g = tmp_path / "X.g01"
    g.write_text("Geom Title=2D only\n\nStorage Area=Area 1          ,,\nStorage Area Surface Line= 4 \n"
                 "              0.              0.\n             10.              0.\n             10.             10.\n"
                 "              0.             10.\nStorage Area Is2D=-1\nStorage Area Point Generation Data=,,2,2\n"
                 "Storage Area 2D Points= 25 \nStorage Area Mannings=0.03\n\nBC Line Name=Inflow\nBC Line Storage Area=Area 1\n"
                 "BC Line Arc= 2 \n              0.              0.              0.             10.\n")
    G = read_geometry(g, [])
    assert G["counts"] == {"2D flow areas": 1, "BC lines": 1}
    sa = G["storage"][0]
    assert sa["is2d"] and sa["cells_text"] == 25 and sa["polygon_area"] == 100 and sa["spacing"] == "2 x 2"
    assert G["bclines"][0]["_line"] == [(0.0, 0.0), (0.0, 10.0)]


@pytest.mark.parametrize("tool", ["dfmvue", "dfmvue23", "epanetvue", "swmmvue", "flovue", "xbvue", "_report"])
def test_theme_switch_in_every_template(tool):
    # every HTML report carries the light / dark switch; the dark colours no longer depend on the system alone
    if tool == "flovue":
        pytest.importorskip("flopy")
    tpl = importlib.import_module(f"modelvue.{tool}").HTML_TEMPLATE
    assert 'class="theme"' in tpl and "modelvue-theme" in tpl and ":root[data-theme=dark]" in tpl
    assert "prefers-color-scheme" not in tpl          # light by default, whatever the system setting
