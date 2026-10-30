# geodesy: writes tables.tsl from PROJ's database (through pyproj). Run from
# the project folder:
#
#   .venv/bin/python tools/gen/gen.py
#
# It writes:
#   ellipsoids      every ellipsoid PROJ strings can name (+ellps=…)
#   units           every length unit (+units=…)
#   meridians       every prime meridian (+pm=…)
#   geographic      the EPSG geographic 2D systems: ellipsoid, prime
#                   meridian, angle unit, and the shift to WGS84 that PROJ
#                   itself picks, when it's a single Helmert (or none)
#   projected       the EPSG projected systems whose method geodesy has:
#                   the PROJ string of the projection, its units and axes,
#                   and its geographic system
#
# Author: Alexius Academia
# Date: 2026-10-31

import os
import re
import warnings

warnings.filterwarnings("ignore")
import pyproj  # noqa: E402
from pyproj import CRS  # noqa: E402
from pyproj.database import query_crs_info  # noqa: E402
from pyproj.enums import PJType  # noqa: E402
from pyproj.transformer import TransformerGroup  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
METHODS = {"tmerc", "utm", "merc", "lcc", "aea", "stere", "sterea", "laea", "eqc", "cass", "ups"}
WGS84 = CRS.from_epsg(4326)


def tsl(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"').replace("{", "\\{").replace("}", "\\}") + '"'


def num(x):
    x = float(x)
    if x == int(x) and abs(x) < 1e15:
        return str(int(x))
    return repr(x)


def proj_params(text):
    out = {}
    for tok in text.split():
        tok = tok.lstrip("+")
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
        else:
            out[tok] = ""
    return out


def helmert_of(definition):
    """The position-vector Helmert (x y z rx ry rz s) of a PROJ pipeline
    that is one Helmert between geocentric steps, [] for a no-op, or None
    for anything else (grids, several steps)."""
    d = definition.replace("+", "")
    if d.startswith("proj=noop"):
        return []
    if "gridshift" in d or "nadgrids" in d or "geoidgrids" in d or "deformation" in d:
        return None
    steps = d.split(" step ")
    hs = [s for s in steps if "proj=helmert" in s]
    if len(hs) != 1:
        return None
    for s in steps:
        p = proj_params(s.replace("step", ""))
        name = p.get("proj")
        if name not in ("pipeline", "unitconvert", "push", "pop", "cart", "helmert", "axisswap", "longlat"):
            return None
    h = hs[0]
    p = proj_params(h.replace("step", ""))
    inv = "inv" in p
    vals = [float(p.get(k, 0)) for k in ("x", "y", "z", "rx", "ry", "rz", "s")]
    if p.get("convention", "position_vector") == "coordinate_frame":
        vals[3], vals[4], vals[5] = -vals[3], -vals[4], -vals[5]
    if inv:
        vals = [-v for v in vals]
    return vals


def shift_meridian(definition):
    """The prime meridian a shift's pipeline puts longitudes on before its
    Helmert ("inv proj=longlat … pm=…"), when it names one; None if not."""
    d = definition.replace("+", "")
    for st in d.split(" step "):
        p = proj_params(st.replace("step", ""))
        if p.get("proj") == "longlat" and "pm" in p:
            v = p["pm"]
            try:
                return float(v)
            except ValueError:
                return dms_degrees(pyproj.get_prime_meridians_map()[v])
    return None


def dms_degrees(text):
    """PROJ's DMS, like 2d20'14.025"E, in degrees."""
    m = re.match(r"^\s*([-+]?\d+(?:\.\d+)?)(?:d(\d+(?:\.\d+)?)'?)?(?:(\d+(?:\.\d+)?)\")?\s*([NSEW])?\s*$", text)
    if not m:
        return float(text)
    v = float(m.group(1)) + float(m.group(2) or 0) / 60 + float(m.group(3) or 0) / 3600
    return -v if m.group(4) in ("S", "W") else v


def proj_pm(c):
    """The prime meridian as PROJ itself uses it (its named meridians can
    differ from EPSG's numbers in the 9th decimal: Paris)."""
    p = proj_params(c.to_proj4())
    pm = p.get("pm")
    if pm is None:
        return 0.0
    try:
        return float(pm)
    except ValueError:
        return dms_degrees(pyproj.get_prime_meridians_map()[pm])


def geographic_rows():
    rows = []
    shifts = {}
    infos = query_crs_info(auth_name="EPSG", pj_types=[PJType.GEOGRAPHIC_2D_CRS], allow_deprecated=False)
    for i in infos:
        code = int(i.code)
        try:
            c = CRS.from_epsg(code)
        except Exception:
            continue
        e = c.ellipsoid
        a = e.semi_major_metre
        rf = e.inverse_flattening if e.inverse_flattening and e.inverse_flattening == e.inverse_flattening else 0.0
        if not e.is_semi_minor_computed and rf == 0.0 and e.semi_minor_metre != a:
            rf = a / (a - e.semi_minor_metre)
        # The numbers PROJ itself works with: its pipelines name PROJ's own
        # ellipsoids (+ellps=bess_nam: a = 6377483.865, not EPSG's
        # 6377483.865280419).
        p4 = proj_params(c.to_proj4())
        named = p4.get("ellps") or {"WGS84": "WGS84", "NAD83": "GRS80", "NAD27": "clrk66"}.get(p4.get("datum", ""))
        if named and named in pyproj.get_ellps_map():
            v = pyproj.get_ellps_map()[named]
            a = float(v["a"])
            if "rf" in v:
                rf = float(v["rf"])
            elif "b" in v:
                b = float(v["b"])
                rf = 0.0 if b == a else a / (a - b)
        pmdeg = proj_pm(c)
        unit = c.axis_info[0].unit_conversion_factor  # radians per unit
        per_degree = unit / (3.141592653589793 / 180)
        shift = None
        known = 0
        spm = None
        try:
            tg = TransformerGroup(c, WGS84, always_xy=True)
            for t in tg.transformers:
                h = helmert_of(t.definition)
                if h is not None and "Ballpark" not in t.description:
                    shift = h
                    known = 1
                    spm = shift_meridian(t.definition)
                    break
        except Exception:
            pass
        if shift is None:
            shift = []
        if spm is not None and abs(spm - pmdeg) < 1e-12:
            spm = None
        rows.append("|".join([str(code), i.name.replace("|", "/"), num(a), num(rf), num(pmdeg), num(per_degree),
                              str(known), ",".join(num(v) for v in shift), "" if spm is None else num(spm)]))
        shifts[code] = shift
    return rows


def norm_name(n):
    return "".join(ch for ch in n.lower() if ch.isalnum())


def datum_rows(geo_rows):
    """Datum names, as EPSG and Esri write them, to a geographic code: for
    WKT that names its datum but gives no shift."""
    seen = {}
    for row in geo_rows:
        f = row.split("|")
        code = int(f[0])
        try:
            c = CRS.from_epsg(code)
        except Exception:
            continue
        names = {c.datum.name, c.name}
        try:
            esri = c.to_wkt("WKT1_ESRI")
            m = re.search(r'DATUM\["([^"]+)"', esri)
            if m:
                names.add(m.group(1))
            m = re.search(r'GEOGCS\["([^"]+)"', esri)
            if m:
                names.add(m.group(1))
        except Exception:
            pass
        for n in names:
            k = norm_name(n)
            # The first (lowest) code for a name wins: usually the main one.
            if k and k not in seen:
                seen[k] = code
    return ["{}|{}".format(k, v) for k, v in sorted(seen.items())]


def projected_rows(geo_codes):
    rows = []
    skipped = {}
    infos = query_crs_info(auth_name="EPSG", pj_types=[PJType.PROJECTED_CRS], allow_deprecated=False)
    for i in infos:
        code = int(i.code)
        try:
            c = CRS.from_epsg(code)
            p4 = c.to_proj4()
        except Exception:
            continue
        p = proj_params(p4)
        method = p.get("proj")
        if method not in METHODS:
            skipped[method] = skipped.get(method, 0) + 1
            continue
        base = c.geodetic_crs.to_epsg() if c.geodetic_crs else None
        if base is None or base not in geo_codes:
            skipped["base " + str(base)] = skipped.get("base " + str(base), 0) + 1
            continue
        # Web Mercator: Mercator on a sphere of radius a, from the base's
        # latitudes as they are (EPSG:3857 and its kin).
        if method == "merc" and "a" in p and "b" in p and p["a"] == p["b"]:
            p4 = p4.replace("+proj=merc", "+proj=webmerc +R=" + p["a"])
        # The projection part only: the ellipsoid, datum and prime meridian
        # come from the base.
        keep = []
        for tok in p4.split():
            k = tok.lstrip("+").split("=")[0]
            if k in ("ellps", "datum", "towgs84", "a", "b", "rf", "f", "pm", "no_defs", "type", "nadgrids",
                     "wktext", "units") and not (k == "units" and p.get("units") != "m"):
                continue
            keep.append(tok.lstrip("+"))
        # The axis directions, when they aren't east and north.
        dirs = "".join(a.direction[0] for a in c.axis_info[:2])
        order = {"en": "", "ne": "", "ws": "wsu", "sw": "wsu", "es": "esu", "se": "esu", "wn": "wnu", "nw": "wnu"}
        axis = order.get(dirs, "")
        if axis and not any(t.startswith("axis=") for t in keep):
            keep.append("axis=" + axis)
        rows.append("|".join([str(code), i.name.replace("|", "/"), str(base), " ".join(keep)]))
    return rows, skipped


def main():
    ellps = []
    for name, v in sorted(pyproj.get_ellps_map().items()):
        a = float(v["a"])
        if "rf" in v:
            rf = float(v["rf"])
        elif "b" in v:
            b = float(v["b"])
            rf = 0.0 if b == a else a / (a - b)
        else:
            rf = 0.0
        ellps.append("|".join([name, num(a), num(rf), v.get("description", "")]))
    # PROJ's +units names (its units.cpp), with exact values.
    proj_units = [("km", 1000.0, "kilometre"), ("m", 1.0, "metre"), ("dm", 0.1, "decimetre"),
                  ("cm", 0.01, "centimetre"), ("mm", 0.001, "millimetre"), ("kmi", 1852.0, "international nautical mile"),
                  ("in", 0.0254, "international inch"), ("ft", 0.3048, "international foot"),
                  ("yd", 0.9144, "international yard"), ("mi", 1609.344, "international statute mile"),
                  ("fath", 1.8288, "international fathom"), ("ch", 20.1168, "international chain"),
                  ("link", 0.201168, "international link"), ("us-in", 1.0 / 39.37, "U.S. surveyor's inch"),
                  ("us-ft", 1200.0 / 3937.0, "U.S. surveyor's foot"), ("us-yd", 3600.0 / 3937.0, "U.S. surveyor's yard"),
                  ("us-ch", 79200.0 / 3937.0, "U.S. surveyor's chain"), ("us-mi", 6336000.0 / 3937.0, "U.S. surveyor's statute mile"),
                  ("ind-yd", 0.91439523, "Indian yard"), ("ind-ft", 0.30479841, "Indian foot"),
                  ("ind-ch", 20.11669506, "Indian chain")]
    units = ["|".join([n, repr(v), d]) for n, v, d in proj_units]
    meridians = []
    for name, v in sorted(pyproj.get_prime_meridians_map().items()):
        # Values like "2d20'14.025\"E" or "2.33722917".
        meridians.append("|".join([name, v]))
    geo = geographic_rows()
    geo_codes = {int(r.split("|")[0]) for r in geo}
    proj, skipped = projected_rows(geo_codes)
    out = ["// geodesy: tables from PROJ {} (pyproj {}): ellipsoids, units, prime".format(
        pyproj.proj_version_str, pyproj.__version__),
        "// meridians, and the EPSG geographic and projected systems geodesy can use. The",
        "// EPSG data are © IOGP, from the EPSG Geodetic Parameter Dataset.",
        "// Written by tools/gen/gen.py; don't edit it by hand.",
        "//",
        "// Author: Alexius Academia",
        "// Date: 2026-10-31",
        ""]

    def table(name, comment, rows):
        out.append("// " + comment)
        out.append("let " + name + ": [String] = [")
        for r in rows:
            out.append("    " + tsl(r) + ",")
        out.append("]")
        out.append("")

    table("ellipsoidRows", "name|a|inverse flattening (0: a sphere)|description", ellps)
    table("unitRows", "name|metres|description", units)
    table("meridianRows", "name|longitude east of Greenwich", meridians)
    table("geographicRows", "code|name|a|1/f|prime meridian (°)|degrees per angle unit|shift known|"
          "x,y,z (m),rx,ry,rz (\"),s (ppm) to WGS84, position vector|the meridian the shift uses, if not the "
          "system's", geo)
    table("projectedRows", "code|name|geographic code|PROJ string of the projection", proj)
    datums = datum_rows(geo)
    table("datumRows", "datum or geographic system name (lower case, letters and digits)|geographic code", datums)
    open(os.path.join(ROOT, "tables.tsl"), "w").write("\n".join(out) + "\n")
    print("tables.tsl: {} ellipsoids, {} units, {} meridians, {} geographic, {} projected, {} datum names; left out: {}".format(
        len(ellps), len(units), len(meridians), len(geo), len(proj), len(datums), skipped))
    print("size", os.path.getsize(os.path.join(ROOT, "tables.tsl")))


if __name__ == "__main__":
    main()
