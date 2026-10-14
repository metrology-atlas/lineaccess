# SPDX-License-Identifier: MIT
"""Accesses to an invasive pressure line (blood draws, flushes, zeroing), found in the pressure waveform.

When a line is accessed the transducer stops seeing the patient for a while: a flush shows the pressure
bag (about 300 mmHg), a stopcock open to air or zeroing shows about 0 mmHg, a closed stopcock or a clot
shows a pressure with no pulse. These look nothing like a pulse, so they can be found from the waveform
alone, whatever the monitor or dataset. Nothing here assumes a sampling rate, a dataset or a protocol;
thresholds are parameters (:class:`LineRules`). Values must be in mmHg.

``line_events``        short windows classified "high", "low" or "flat", merged into events.
``line_accesses``      events grouped into accesses (one access = everything done in one visit to the
                       line), each with a recovery period until the pulse is back to its size before;
                       occlusions (the pulse lost with no flush and no zero, e.g. a blood-pressure cuff
                       inflating on the same arm); and long stretches without a usable signal
                       ("unavailable": line not connected yet, removed, transducer off).
``match_times``        pairs accesses with external times (e.g. lab results) that follow them, and the lag.
``labels``             (name, start, end) tuples for annotation tracks.

Only numpy.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

__all__ = ["LineAccess", "LineEvent", "LineRules", "labels", "line_accesses", "line_events", "match_times"]


@dataclass
class LineRules:
    """Thresholds (mmHg, seconds)."""

    window_s: float = 1.0          # classification window
    high_mmhg: float = 200.0       # median above: flush / pressure bag
    low_mmhg: float = 10.0         # median below: open to air, zeroing, disconnected
    flat_pp_mmhg: float = 4.0      # pulse pressure (5-95%) below, with a pressure in range: no pulse
    merge_s: float = 2.0           # windows of the same event at most this far apart
    group_s: float = 90.0          # events of one access at most this far apart
    max_access_s: float = 300.0    # longer than this: the line is unavailable, not being accessed
    recovery_ref_s: float = 60.0   # pulse size before the access is measured over this
    recovery_frac: float = 0.7     # recovered when the pulse is back to this share of it
    recovery_max_s: float = 180.0  # give up after this (the access ends "not recovered")
    min_occlusion_s: float = 4.0   # a lost pulse alone (no flush, no zero) shorter than this is ignored


@dataclass
class LineEvent:
    start: float
    end: float
    kinds: list[str] = field(default_factory=list)
    median_mmhg: float = float("nan")


@dataclass
class LineAccess:
    start: float
    end: float                     # end of the last event
    recovered: float               # end of the recovery (>= end)
    kind: str                      # "access", "occlusion" or "unavailable"
    n_events: int = 0
    kinds: list[str] = field(default_factory=list)
    max_mmhg: float = float("nan")
    recovered_ok: bool = True

    def to_dict(self) -> dict:
        return asdict(self)


def _windows(t: np.ndarray, x: np.ndarray, fs: float, win_s: float
             ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Window start times, median, pulse pressure (5-95%) and pulse pressure after removing each window's
    straight-line trend (a slow ramp with no pulse, e.g. a line closed to the patient filling from the
    flush device, has a large plain spread but almost none after detrending)."""
    w = max(round(win_s * fs), 2)
    n = x.size // w
    X = x[: n * w].reshape(n, w)
    u = np.arange(w) - (w - 1) / 2
    with np.errstate(all="ignore"):
        med = np.nanmedian(X, axis=1)
        pp = np.nanquantile(X, 0.95, axis=1) - np.nanquantile(X, 0.05, axis=1)
        Xf = np.where(np.isfinite(X), X, np.nanmean(X, axis=1, keepdims=True))
        slope = (Xf - Xf.mean(axis=1, keepdims=True)) @ u / (u @ u)
        D = Xf - Xf.mean(axis=1, keepdims=True) - slope[:, None] * u[None, :]
        ppd = np.quantile(D, 0.95, axis=1) - np.quantile(D, 0.05, axis=1)
    return t[: n * w : w], med, pp, ppd


def line_events(x: np.ndarray, fs: float, *, t: np.ndarray | None = None, t0: float = 0.0,
                rules: LineRules | None = None) -> list[LineEvent]:
    """Events in a pressure waveform ``x`` (mmHg) at ``fs`` Hz. ``t``: sample times (s), if the stream has
    breaks; else ``t0`` + n / fs. Windows across a break in ``t`` are not merged."""
    r = rules or LineRules()
    x = np.asarray(x, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64) if t is not None else t0 + np.arange(x.size) / fs
    ts, med, _, ppd = _windows(t, x, fs, r.window_s)
    high = med > r.high_mmhg
    low = med < r.low_mmhg
    flat = (ppd < r.flat_pp_mmhg) & ~high & ~low
    out: list[LineEvent] = []
    for i in np.flatnonzero(high | low | flat):
        if not np.isfinite(med[i]):
            continue
        k = "high" if high[i] else "low" if low[i] else "flat"
        a, b = float(ts[i]), float(ts[i]) + r.window_s
        if out and a - out[-1].end <= r.merge_s:
            e = out[-1]
            e.end = b
            if k not in e.kinds:
                e.kinds.append(k)
            e.median_mmhg = max(e.median_mmhg, float(med[i])) if k == "high" else e.median_mmhg
        else:
            out.append(LineEvent(a, b, [k], float(med[i])))
    return out


def line_accesses(x: np.ndarray, fs: float, *, t: np.ndarray | None = None, t0: float = 0.0,
                  rules: LineRules | None = None, events: list[LineEvent] | None = None) -> list[LineAccess]:
    """Group :func:`line_events` into accesses (events within ``group_s``) with their recovery, and mark
    stretches longer than ``max_access_s`` as "unavailable"."""
    r = rules or LineRules()
    x = np.asarray(x, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64) if t is not None else t0 + np.arange(x.size) / fs
    ev = events if events is not None else line_events(x, fs, t=t, rules=r)
    groups: list[list[LineEvent]] = []
    for e in ev:
        if groups and e.start - groups[-1][-1].end <= r.group_s:
            groups[-1].append(e)
        else:
            groups.append([e])
    ts, _, pp, _ = _windows(t, x, fs, r.window_s)
    t_last = float(t[-1]) if t.size else 0.0
    out: list[LineAccess] = []
    for g in groups:
        a, b = g[0].start, g[-1].end
        kinds = sorted({k for e in g for k in e.kinds})
        mx = max((e.median_mmhg for e in g if "high" in e.kinds), default=float("nan"))
        if b - a > r.max_access_s:
            out.append(LineAccess(a, b, b, "unavailable", len(g), kinds, mx, True))
            continue
        # a lost pulse with no flush and no zero: an occlusion upstream (e.g. a blood-pressure cuff on the
        # same arm, a kink), not an access; very short ones are noise
        kind = "occlusion" if kinds == ["flat"] else "access"
        if kind == "occlusion" and sum(e.end - e.start for e in g) < r.min_occlusion_s:
            continue
        # recovery: until the pulse pressure is back to recovery_frac of its median before the access
        pre = (ts >= a - r.recovery_ref_s) & (ts < a)
        ref = float(np.nanmedian(pp[pre])) if pre.any() else float("nan")
        rec, ok = b, True
        if np.isfinite(ref) and ref > r.flat_pp_mmhg:
            post = np.flatnonzero((ts >= b) & (ts < b + r.recovery_max_s))
            good = set(post[pp[post] >= r.recovery_frac * ref].tolist())
            # recovered at the first window that starts 5 good windows in a row
            run = 0
            rec, ok = b + r.recovery_max_s, False
            for i in post:
                run = run + 1 if int(i) in good else 0
                if run >= 5:
                    rec, ok = float(ts[i - 4]), True
                    break
        rec = min(rec, t_last)
        if not ok and t_last - b < r.recovery_max_s:
            # the signal never came back before the recording ended: the line was taken off
            out.append(LineAccess(a, t_last, t_last, "unavailable", len(g), kinds, mx, False))
            continue
        out.append(LineAccess(a, b, max(rec, b), kind, len(g), kinds, mx, ok))
    return out


def match_times(accesses: list[LineAccess], times: np.ndarray, *, lag_range_s: tuple[float, float] = (0.0, 900.0)
                ) -> list[dict]:
    """For each external time (e.g. a lab result), the latest access that starts ``lag_range_s`` before it.
    Returns one dict per time: {"time", "access" (index or None), "lag_s"}."""
    starts = np.array([x.start for x in accesses if x.kind == "access"])
    idx = [i for i, x in enumerate(accesses) if x.kind == "access"]
    out = []
    for tm in np.asarray(times, dtype=np.float64):
        lag = tm - starts
        ok = np.flatnonzero((lag >= lag_range_s[0]) & (lag <= lag_range_s[1]))
        if ok.size:
            j = ok[np.argmin(lag[ok])]
            out.append({"time": float(tm), "access": idx[int(j)], "lag_s": float(lag[j])})
        else:
            out.append({"time": float(tm), "access": None, "lag_s": None})
    return out


def labels(accesses: list[LineAccess], prefix: str = "art_line") -> list[tuple[str, float, float]]:
    """("art_line:access" or "art_line:occlusion", start, recovered) and ("art_line:unavailable", start, end)
    for annotation."""
    return [(f"{prefix}:{x.kind}", x.start, x.end if x.kind == "unavailable" else x.recovered) for x in accesses]
