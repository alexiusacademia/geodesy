# geodesy: writes tests/cases.json, the cases the library is checked
# against, with the answers of GeographicLib (geodesics), PROJ through
# pyproj (projections, datum shifts, PROJ strings, WKT), the mgrs package
# (MGRS), and SciPy's quadrature (rhumb lines). Run from the project folder:
#
#   .venv/bin/python tools/oracle/oracle.py
#
# Author: Alexius Academia
# Date: 2026-10-31

import json
import math
import os
import random
import re
import warnings

warnings.filterwarnings("ignore")
import mgrs as mgrslib  # noqa: E402
from geographiclib.geodesic import Geodesic  # noqa: E402
from pyproj import CRS, Transformer  # noqa: E402
from scipy.integrate import quad  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SEED = 20261031
NAN = 1e308


def enc(x):
    return NAN if (x is None or (isinstance(x, float) and not math.isfinite(x))) else float(x)


def rows(name):
    """Rows of a table in tables.tsl."""
    text = open(os.path.join(ROOT, "tables.tsl")).read()
    start = text.index("let " + name + ": [String] = [")
    end = text.index("\n]", start)
    out = []
    for line in text[start:end].splitlines()[1:]:
        line = line.strip().rstrip(",")
        if line.startswith('"'):
            out.append(json.loads(line))
    return out


# --- Geodesics.

ELLIPSOIDS = [("wgs84", 6378137.0, 1 / 298.257223563), ("clrk66", 6378206.4, 1 / 294.9786982138982),
              ("sphere", 6371000.0, 0.0), ("flat", 6378137.0, 1 / 150.0), ("prolate", 6378137.0, -1 / 150.0)]


def random_point(rng):
    lat = math.degrees(math.asin(rng.uniform(-1, 1)))
    lon = rng.uniform(-180, 180)
    return round(lat, 9), round(lon, 9)


def geodesic_cases(rng):
    out = []
    for name, a, f in ELLIPSOIDS:
        g = Geodesic(a, f)
        n = 1200 if name == "wgs84" else 300
        for k in range(n):
            kind = k % 8
            lat1, lon1 = random_point(rng)
            lat2, lon2 = random_point(rng)
            if kind == 1:  # nearly antipodal
                lat2 = round(-lat1 + rng.uniform(-0.5, 0.5), 9)
                lon2 = round(lon1 + 180 + rng.uniform(-0.5, 0.5), 9)
            elif kind == 2:  # short
                lat2 = round(lat1 + rng.uniform(-0.01, 0.01), 9)
                lon2 = round(lon1 + rng.uniform(-0.01, 0.01), 9)
                lat2 = max(-90, min(90, lat2))
            elif kind == 3:  # equator
                lat1 = 0.0
                lat2 = 0.0
            elif kind == 4:  # a meridian
                lon2 = lon1
            elif kind == 5:  # a pole
                lat1 = rng.choice([90.0, -90.0])
            elif kind == 6:  # the same point
                lat2, lon2 = lat1, lon1
            lon2 = ((lon2 + 180) % 360) - 180
            r = g.Inverse(lat1, lon1, lat2, lon2, Geodesic.ALL)
            d = g.Direct(lat1, lon1, r["azi1"] if kind != 6 else rng.uniform(-180, 180),
                         rng.uniform(-2e7, 2e7), Geodesic.ALL)
            out.append({"ellipsoid": name, "lat1": lat1, "lon1": lon1, "lat2": lat2, "lon2": lon2,
                        "s12": r["s12"], "azi1": r["azi1"], "azi2": r["azi2"], "a12": r["a12"], "m12": r["m12"],
                        "M12": r["M12"], "M21": r["M21"], "S12": r["S12"],
                        "dazi1": d["azi1"], "ds12": d["s12"], "dlat2": d["lat2"], "dlon2": d["lon2"],
                        "dazi2": d["azi2"]})
    return out


def line_cases(rng):
    out = []
    g = Geodesic.WGS84
    for _ in range(200):
        lat1, lon1 = random_point(rng)
        azi = rng.uniform(-180, 180)
        l = g.Line(lat1, lon1, azi, Geodesic.ALL | Geodesic.LONG_UNROLL)
        dists = [rng.uniform(-4e7, 4e7) for _ in range(5)]
        pts = []
        for s in dists:
            p = l.Position(s, Geodesic.STANDARD | Geodesic.LONG_UNROLL)
            pts.append([s, p["lat2"], p["lon2"], p["azi2"]])
        out.append({"lat1": lat1, "lon1": lon1, "azi1": azi, "points": pts})
    return out


def polygon_cases(rng):
    out = []
    g = Geodesic.WGS84
    for k in range(300):
        kind = k % 4
        n = rng.randint(3, 12)
        if kind == 0:  # a lot or a field: metres to kilometres
            clat, clon = random_point(rng)
            clat = max(-80, min(80, clat))
            size = 10 ** rng.uniform(-4, -1)
        elif kind == 1:  # a province or a country
            clat, clon = random_point(rng)
            clat = max(-70, min(70, clat))
            size = rng.uniform(0.5, 15)
        elif kind == 2:  # across the antimeridian
            clat, clon = rng.uniform(-60, 60), 180.0
            size = rng.uniform(0.1, 10)
        else:  # round a pole
            clat, clon = rng.choice([89.0, -89.0]), 0.0
            size = 0.0
        pts = []
        angs = sorted(rng.uniform(0, 2 * math.pi) for _ in range(n))
        for t in angs:
            if kind == 3:
                lat = clat + rng.uniform(-0.5, 0.5)
                lon = math.degrees(t) - 180
            else:
                r = size * rng.uniform(0.5, 1.0)
                lat = max(-89.9, min(89.9, clat + r * math.sin(t)))
                lon = clon + r * math.cos(t) / max(0.05, math.cos(math.radians(clat)))
            lon = ((lon + 180) % 360) - 180
            pts.append((round(lat, 9), round(lon, 9)))
        if rng.random() < 0.5:
            pts.reverse()
        p = g.Polygon()
        for lat, lon in pts:
            p.AddPoint(lat, lon)
        num, perim, area = p.Compute(False, True)
        line = g.Polygon(True)
        for lat, lon in pts:
            line.AddPoint(lat, lon)
        _, track, _ = line.Compute()
        out.append({"lats": [a for a, _ in pts], "lons": [b for _, b in pts], "area": area, "perimeter": perim,
                    "track": track})
    return out


# --- Rhumb lines, by SciPy's quadrature.

def rhumb_cases(rng):
    a = 6378137.0
    f = 1 / 298.257223563
    e2 = f * (2 - f)
    e = math.sqrt(e2)

    def psi(p):
        return math.asinh(math.tan(p)) - e * math.atanh(e * math.sin(p))

    def arc(p1, p2):
        return quad(lambda t: a * (1 - e2) / (1 - e2 * math.sin(t) ** 2) ** 1.5, p1, p2, epsabs=0, epsrel=1e-13,
                    limit=400)[0]

    def dpsi_int(p1, p2):
        return quad(lambda t: (1 - e2) / ((1 - e2 * math.sin(t) ** 2) * math.cos(t)), p1, p2, epsabs=0,
                    epsrel=1e-13, limit=400)[0]

    out = []
    for k in range(300):
        lat1 = rng.uniform(-85, 85)
        lon1 = rng.uniform(-180, 180)
        if k % 5 == 0:
            lat2 = lat1 + rng.uniform(-1e-6, 1e-6)
        else:
            lat2 = rng.uniform(-85, 85)
        lon2 = rng.uniform(-180, 180)
        lat1, lon1, lat2, lon2 = (round(v, 9) for v in (lat1, lon1, lat2, lon2))
        p1, p2 = math.radians(lat1), math.radians(lat2)
        dl = math.radians(((lon2 - lon1) + 180) % 360 - 180)
        dpsi = psi(p2) - psi(p1)
        if abs(p2 - p1) < 1e-12:
            m = (p1 + p2) / 2
            d = a / math.sqrt(1 - e2 * math.sin(m) ** 2) * math.cos(m)
        else:
            d = arc(p1, p2) / dpsi_int(p1, p2)
        s = abs(d) * math.hypot(dpsi, dl)
        azi = math.degrees(math.atan2(dl, dpsi))
        out.append({"lat1": lat1, "lon1": lon1, "lat2": lat2, "lon2": lon2, "s12": s, "azi12": azi})
    return out


# --- Projections and datums, by PROJ.

def area_points(rng, crs, k):
    aou = crs.area_of_use
    if aou is None:
        w, s, e, n = -180, -80, 180, 80
    else:
        w, s, e, n = aou.bounds
    if e < w:
        e += 360
    pts = []
    for _ in range(k):
        lon = rng.uniform(w, e)
        lat = rng.uniform(max(s, -89.5), min(n, 89.5))
        lon = ((lon + 180) % 360) - 180
        pts.append((round(lon, 9), round(lat, 9)))
    return pts


def xy_transformer(src, dst):
    """A function from (lon, lat) or (x, y) in src, east/west value first,
    to dst's (x, y), east/west value first. pyproj's always_xy does this,
    except for systems listed southing-or-northing first and then
    westing-or-easting (Gusterberg's south, west), which it leaves as they
    are."""
    t = Transformer.from_crs(src, dst, always_xy=True)
    dirs = [a.direction.lower() for a in dst.axis_info[:2]]
    swap = len(dirs) == 2 and dirs[0] in ("north", "south") and dirs[1] in ("east", "west")
    tf = Transformer.from_crs(src, dst, always_xy=False)
    src_ns = src.axis_info and src.axis_info[0].direction.lower() in ("north", "south")

    def go(x, y):
        if swap:
            a, b = (y, x) if src_ns else (x, y)
            u, v = tf.transform(a, b)
            return v, u
        return t.transform(x, y)
    return go


GEO = None


def base_of(code):
    """The geographic row of a projected system's base: its prime meridian
    and degrees per unit, as tables.tsl has them."""
    global GEO
    if GEO is None:
        GEO = {int(r.split("|")[0]): r.split("|") for r in rows("geographicRows")}
    return GEO[code]


def projection_cases(rng):
    out = []
    for row in rows("projectedRows"):
        code = int(row.split("|")[0])
        try:
            c = CRS.from_epsg(code)
            g = base_of(int(row.split("|")[2]))
            pmdeg = float(g[4])
            unit = float(g[5])
            fwd = xy_transformer(c.geodetic_crs, c)
        except Exception:
            continue
        pts = []
        for lon, lat in area_points(rng, c, 3):
            try:
                x, y = fwd((lon - pmdeg) / unit, lat / unit)
            except Exception:
                continue
            if math.isfinite(x) and math.isfinite(y):
                pts.append([lon, lat, x, y])
        if pts:
            out.append({"code": code, "points": pts})
    return out


def shift_of(row):
    f = row.split("|")
    return [float(v) for v in f[7].split(",")] if f[7] else []


def datum_cases(rng):
    """Each geographic system to WGS84: by the stored Helmert through PROJ,
    and by PROJ's own default where its first choice is that Helmert."""
    out = []
    for row in rows("geographicRows"):
        f = row.split("|")
        code = int(f[0])
        if f[6] != "1" or not f[7]:
            continue
        shift = shift_of(row)
        try:
            c = CRS.from_epsg(code)
            p4 = c.to_proj4()
            spm = float(f[8]) if len(f) > 8 and f[8] else None
            if spm is not None:
                p4 = re.sub(r"\+pm=\S+", "+pm=" + repr(spm), p4)
            bound = CRS.from_proj4(p4 + " +towgs84=" + ",".join(repr(v) for v in shift))
            t = Transformer.from_crs(bound, "+proj=longlat +datum=WGS84", always_xy=True)
            from pyproj.transformer import TransformerGroup
            tg = TransformerGroup(c, CRS.from_epsg(4326), always_xy=True)
            default = tg.transformers[0] if tg.transformers else None
            pmdeg = float(f[4])
            # Only where PROJ's own first choice is the stored Helmert.
            if default is not None:
                import importlib.util, sys
                spec = importlib.util.spec_from_file_location("gen", os.path.join(ROOT, "tools", "gen", "gen.py"))
                gen = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(gen)
                first = gen.helmert_of(default.definition)
                if first is None or [round(v, 9) for v in first] != [round(v, 9) for v in shift]:
                    default = None
            pmdeg = float(f[4])
            unit = float(f[5])
        except Exception:
            continue
        pts = []
        for lon, lat in area_points(rng, c, 3):
            x, y = t.transform(lon - pmdeg, lat)
            # PROJ's default, in the system's own units and axis order.
            try:
                d = default.transform((lon - pmdeg) / unit, lat / unit) if default else (float("nan"), float("nan"))
            except Exception:
                d = (float("nan"), float("nan"))
            pts.append([lon, lat, x, y, enc(d[0]), enc(d[1])])
        out.append({"code": code, "points": pts})
    return out


def cross_cases(rng):
    """Projected to projected across datums, through WGS84 with each
    side's stored Helmert (as BoundCRS in PROJ)."""
    geo = {int(r.split("|")[0]): r for r in rows("geographicRows")}
    proj = rows("projectedRows")
    out = []
    picks = [r for r in proj if geo.get(int(r.split("|")[2]), "").split("|")[6:7] == ["1"]]
    bounds = {}
    for r in picks:
        try:
            aou = CRS.from_epsg(int(r.split("|")[0])).area_of_use
            if aou is not None:
                bounds[r] = aou.bounds
        except Exception:
            pass
    keys = list(bounds)

    def overlaps(p, q):
        w1, s1, e1, n1 = bounds[p]
        w2, s2, e2, n2 = bounds[q]
        return w1 < e2 and w2 < e1 and s1 < n2 and s2 < n1

    for _ in range(400):
        a = rng.choice(keys)
        others = [q for q in rng.sample(keys, min(400, len(keys))) if overlaps(a, q)]
        if not others:
            continue
        b = rng.choice(others)
        ca = int(a.split("|")[0])
        cb = int(b.split("|")[0])
        try:
            A = CRS.from_epsg(ca)
            B = CRS.from_epsg(cb)

            def bound(c, code):
                r = geo[int(c.geodetic_crs.to_epsg())]
                sh = shift_of(r)
                p4 = c.to_proj4()
                rf = r.split("|")
                if len(rf) > 8 and rf[8]:
                    # The shift's own meridian (MGI (Ferro)): skip these.
                    raise ValueError("shift meridian")
                return CRS.from_proj4(p4 + (" +towgs84=" + ",".join(repr(v) for v in sh) if sh else
                                            " +towgs84=0,0,0"))
            t = Transformer.from_crs(bound(A, ca), bound(B, cb), always_xy=True)
            fa = xy_transformer(A.geodetic_crs, A)
            g = base_of(int(a.split("|")[2]))
            pmdeg = float(g[4])
            unit = float(g[5])
        except Exception:
            continue
        pts = []
        w2, s2, e2, n2 = bounds[b]
        for lon, lat in area_points(rng, A, 6):
            if not (w2 <= lon <= e2 and s2 <= lat <= n2):
                continue
            x, y = fa((lon - pmdeg) / unit, lat / unit)
            u, v = t.transform(x, y)
            if all(math.isfinite(q) for q in (x, y, u, v)):
                pts.append([x, y, u, v])
        if pts:
            out.append({"from": ca, "to": cb, "points": pts})
    return out


def text_cases(rng):
    """PROJ strings and WKT of EPSG systems, read back."""
    proj = rows("projectedRows")
    geo = {int(r.split("|")[0]): r for r in rows("geographicRows")}
    out = []
    for _ in range(600):
        r = rng.choice(proj)
        code = int(r.split("|")[0])
        try:
            c = CRS.from_epsg(code)
            g = base_of(int(r.split("|")[2]))
            pmdeg = float(g[4])
            unit = float(g[5])
            fwd = xy_transformer(c.geodetic_crs, c)
        except Exception:
            continue
        sh = shift_of(geo[int(r.split("|")[2])])
        p4 = c.to_proj4() + (" +towgs84=" + ",".join(repr(v) for v in sh) if sh else "")
        texts = [("proj", p4)]
        for flavour in ("WKT1_GDAL", "WKT1_ESRI", "WKT2_2019"):
            try:
                w = c.to_wkt(flavour)
            except Exception:
                continue
            if not w:
                continue
            texts.append((flavour, w))
            # Without its EPSG codes, to read the parameters themselves.
            stripped = re.sub(r',\s*(AUTHORITY|ID)\["EPSG",\s*"?\d+"?\]', "", w)
            stripped = re.sub(r',\s*(USAGE)\[(?:[^\[\]]|\[[^\[\]]*\])*\]', "", stripped)
            texts.append((flavour + " no code", stripped))
        pts = []
        for lon, lat in area_points(rng, c, 2):
            x, y = fwd((lon - pmdeg) / unit, lat / unit)
            if math.isfinite(x) and math.isfinite(y):
                pts.append([lon, lat, x, y])
        if pts:
            out.append({"code": code, "texts": [[k, t] for k, t in texts], "points": pts})
    return out


def mgrs_cases(rng):
    m = mgrslib.MGRS()
    out = []
    for k in range(1500):
        if k % 10 == 0:
            lat = rng.choice([rng.uniform(84, 89.9), rng.uniform(-89.9, -80.1)])
        else:
            lat = rng.uniform(-79.9, 83.9)
        lon = rng.uniform(-180, 179.999)
        lat, lon = round(lat, 7), round(lon, 7)
        p = rng.choice([5, 5, 4, 3, 1, 0])
        try:
            s = m.toMGRS(lat, lon, MGRSPrecision=p)
            if isinstance(s, bytes):
                s = s.decode()
            back = m.toLatLon(s)
        except Exception:
            continue
        out.append({"lat": lat, "lon": lon, "precision": p, "mgrs": s, "lat2": back[0], "lon2": back[1]})
    return out


def main():
    rng = random.Random(SEED)
    cases = {"seed": SEED, "ellipsoids": [[n, a, f] for n, a, f in ELLIPSOIDS]}
    cases["geodesics"] = geodesic_cases(rng)
    cases["lines"] = line_cases(rng)
    cases["polygons"] = polygon_cases(rng)
    cases["rhumbs"] = rhumb_cases(rng)
    cases["projections"] = projection_cases(rng)
    cases["datums"] = datum_cases(rng)
    cases["cross"] = cross_cases(rng)
    cases["texts"] = text_cases(rng)
    cases["mgrs"] = mgrs_cases(rng)
    path = os.path.join(ROOT, "tests", "cases.json")
    with open(path, "w") as f:
        json.dump(cases, f, separators=(",", ":"))
    print("tests/cases.json:", ", ".join("{} {}".format(len(v), k) for k, v in cases.items() if isinstance(v, list)),
          "-", os.path.getsize(path), "bytes")


if __name__ == "__main__":
    main()
