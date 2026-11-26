"""Minimal NetCDF reader: a pure-Python classic NetCDF-3 reader, with netCDF4 / h5py fallbacks for NetCDF-4."""
from __future__ import annotations
import struct
from pathlib import Path
try:
    import numpy as np  # optional
except ImportError:  # pragma: no cover
    np = None


NC_TYPES = {1: ("b", 1), 2: ("c", 1), 3: ("h", 2), 4: ("i", 4), 5: ("f", 4), 6: ("d", 8),
            7: ("B", 1), 8: ("H", 2), 9: ("I", 4), 10: ("q", 8), 11: ("Q", 8)}


class NC3:
    def __init__(self, path: Path):
        self.path = path
        self.f = open(path, "rb")
        magic = self.f.read(4)
        if magic[:3] != b"CDF":
            self.f.close()
            raise ValueError("not a NetCDF-3 file" + (" (NetCDF-4/HDF5)" if magic[:4] == b"\x89HDF" else ""))
        self.ver = magic[3]
        self.isz = 8 if self.ver == 5 else 4         # size of counts
        self.osz = 4 if self.ver == 1 else 8         # size of offsets
        self.numrecs = self._n()
        self.dims = []                               # [(name, length)], length 0 = record
        self.attrs = {}
        self.vars = {}
        self._dim_list()
        self.attrs = self._att_list()
        self._var_list()
        rv = [v for v in self.vars.values() if v["record"]]
        self.recsize = sum(v["vsize"] for v in rv)
        if len(rv) == 1:
            self.recsize = rv[0]["count"] * rv[0]["size"]

    def close(self):
        self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def _u(self, n):
        b = self.f.read(n)
        return int.from_bytes(b, "big")

    def _n(self):
        return self._u(self.isz)

    def _name(self):
        n = self._n()
        s = self.f.read(n)
        self.f.read((4 - n % 4) % 4)
        return s.decode("utf-8", "replace")

    def _dim_list(self):
        tag, n = self._u(4), self._n()
        for _ in range(n if tag == 0x0A else 0):
            self.dims.append((self._name(), self._n()))

    def _att_list(self):
        tag, n = self._u(4), self._n()
        out = {}
        for _ in range(n if tag == 0x0C else 0):
            name = self._name()
            t, cnt = self._u(4), self._n()
            code, sz = NC_TYPES[t]
            raw = self.f.read(cnt * sz)
            self.f.read((4 - (cnt * sz) % 4) % 4)
            if code == "c":
                val = raw.decode("utf-8", "replace").rstrip("\0")
            else:
                val = list(struct.unpack(f">{cnt}{code}", raw))
                val = val[0] if cnt == 1 else val
            out[name] = val
        return out

    def _var_list(self):
        tag, n = self._u(4), self._n()
        for _ in range(n if tag == 0x0B else 0):
            name = self._name()
            nd = self._n()
            dimids = [self._n() for _ in range(nd)]
            attrs = self._att_list()
            t = self._u(4)
            vsize = self._n() if self.ver == 5 else self._u(4)
            begin = self._u(self.osz)
            code, sz = NC_TYPES[t]
            dnames = [self.dims[i][0] for i in dimids]
            shape = [self.dims[i][1] for i in dimids]
            rec = bool(shape) and shape[0] == 0
            count = 1
            for s in (shape[1:] if rec else shape):
                count *= s
            self.vars[name] = {"dims": dnames, "shape": shape, "record": rec, "code": code, "size": sz,
                               "vsize": vsize, "begin": begin, "count": count, "attrs": attrs}

    def dim(self, name):
        for n, l in self.dims:
            if n == name:
                return self.numrecs if l == 0 else l
        return None

    def shape(self, name):
        v = self.vars[name]
        return [self.numrecs] + v["shape"][1:] if v["record"] else v["shape"]

    def _decode(self, raw, code, n):
        if np is not None and code != "c":
            return np.frombuffer(raw, dtype=">" + code, count=n).astype(float if code in "fd" else int).tolist()
        if code == "c":
            return raw
        return list(struct.unpack(f">{n}{code}", raw))

    def read(self, name, rec=None):
        """Whole variable as a flat list (record vars: all records, or one record)."""
        v = self.vars[name]
        n, sz, code = v["count"], v["size"], v["code"]
        if not v["record"]:
            self.f.seek(v["begin"])
            return self._decode(self.f.read(n * sz), code, n)
        recs = range(self.numrecs) if rec is None else [rec if rec >= 0 else self.numrecs + rec]
        out = []
        for r in recs:
            self.f.seek(v["begin"] + r * self.recsize)
            out += self._decode(self.f.read(n * sz), code, n)
        return out

    def strings(self, name):
        v = self.vars[name]
        raw = self.read(name)
        L = v["shape"][-1]
        return [raw[i:i + L].decode("utf-8", "replace").replace("\0", "").strip()
                for i in range(0, len(raw), L)]

    def att(self, var, key, default=None):
        a = self.vars.get(var, {}).get("attrs", {})
        return a.get(key, default)


def open_nc(path: Path):
    """NC3 reader, or a thin adapter around netCDF4 for NetCDF-4 files."""
    try:
        return NC3(path)
    except ValueError as e:
        if "HDF5" not in str(e):
            raise
    try:
        import netCDF4  # type: ignore
    except ImportError:
        raise ValueError(f"{path.name} is NetCDF-4/HDF5 - install netCDF4 to read it")

    class NC4:
        def __init__(s, p):
            s.ds = netCDF4.Dataset(p)
            s.ds.set_auto_mask(False)
            s.path = p
            s.attrs = {k: s.ds.getncattr(k) for k in s.ds.ncattrs()}
            s.dims = [(k, len(d)) for k, d in s.ds.dimensions.items()]
            s.numrecs = next((len(d) for d in s.ds.dimensions.values() if d.isunlimited()), 0)
            s.vars = {k: {"dims": list(v.dimensions), "shape": list(v.shape),
                          "record": bool(v.dimensions) and s.ds.dimensions[v.dimensions[0]].isunlimited(),
                          "attrs": {a: v.getncattr(a) for a in v.ncattrs()}} for k, v in s.ds.variables.items()}

        def close(s): s.ds.close()
        def __enter__(s): return s
        def __exit__(s, *a): s.close()
        def dim(s, n): return len(s.ds.dimensions[n]) if n in s.ds.dimensions else None
        def shape(s, n): return list(s.ds.variables[n].shape)

        def read(s, n, rec=None):
            v = s.ds.variables[n]
            a = v[:] if rec is None else v[rec]
            return a.ravel().tolist()

        def strings(s, n):
            import netCDF4 as N
            return [x.strip() for x in N.chartostring(s.ds.variables[n][:]).tolist()]

        def att(s, var, key, default=None):
            return s.vars.get(var, {}).get("attrs", {}).get(key, default)
    return NC4(path)
