# SPDX-License-Identifier: MIT
import numpy as np

from lineaccess import LineRules, labels, line_accesses, line_events, match_times

FS = 100.0


def abp(minutes=30, seed=0):
    """A pulsatile arterial pressure, 120/70 mmHg at 75 /min."""
    rng = np.random.default_rng(seed)
    t = np.arange(int(minutes * 60 * FS)) / FS
    ph = (t * 1.25) % 1.0
    x = 70 + 50 * np.where(ph < 0.15, ph / 0.15, np.exp(-(ph - 0.15) * 4)) + rng.normal(0, 0.5, t.size)
    return t, x


def put(t, x, a, b, value, noise=1.0, seed=1):
    rng = np.random.default_rng(seed)
    s = (t >= a) & (t < b)
    x[s] = value + rng.normal(0, noise, s.sum())


def test_draw_with_two_flushes_is_one_access_with_recovery():
    t, x = abp()
    put(t, x, 600, 615, 300)            # flush
    put(t, x, 645, 660, 295)            # flush after the sample
    s = (t >= 660) & (t < 700)           # damped for 40 s after: small pulse
    x[s] = 95 + 0.2 * (x[s] - 95)
    acc = line_accesses(x, FS, t0=0.0)
    assert len(acc) == 1
    a = acc[0]
    assert a.kind == "access" and a.n_events == 2 and a.kinds == ["high"]
    assert abs(a.start - 600) <= 1 and abs(a.end - 660) <= 1
    assert a.recovered_ok and 695 <= a.recovered <= 705
    assert labels(acc) == [("art_line:access", a.start, a.recovered)]


def test_zeroing_flat_and_unavailable():
    t, x = abp(minutes=40)
    put(t, x, 0, 400, 0.0, noise=0.3)   # not connected for the first minutes
    put(t, x, 1200, 1210, 0.0)          # zeroing
    put(t, x, 1800, 1830, 85.0, noise=0.5)  # stopcock closed: pressure, no pulse
    acc = line_accesses(x, FS)
    kinds = [(a.kind, a.kinds) for a in acc]
    assert kinds == [("unavailable", ["low"]), ("access", ["low"]), ("occlusion", ["flat"])]


def test_no_events_in_a_clean_trace():
    _t, x = abp(seed=3)
    assert line_events(x, FS) == []


def test_breaks_in_sample_times_and_rules():
    t, x = abp(minutes=10)
    t = t + np.where(t > 300, 20.0, 0.0)  # a 20 s break
    put(t, x, 450, 460, 250)
    ev = line_events(x, FS, t=t)
    assert len(ev) == 1 and abs(ev[0].start - 450) <= 1
    # with a higher flush threshold it is still an event: a pressure with no pulse
    assert [e.kinds for e in line_events(x, FS, t=t, rules=LineRules(high_mmhg=260))] == [["flat"]]


def test_match_times_lag():
    t, x = abp()
    put(t, x, 600, 615, 300)
    put(t, x, 1300, 1315, 300)
    acc = line_accesses(x, FS)
    m = match_times(acc, [780.0, 1490.0, 100.0])
    assert m[0]["access"] == 0 and abs(m[0]["lag_s"] - 180) <= 1
    assert m[1]["access"] == 1
    assert m[2]["access"] is None


def test_slow_ramp_closed_to_patient_starts_the_access():
    # stopcock closed to the patient: no pulse, pressure rising slowly from the flush device, then a flush
    t, x = abp()
    s = (t >= 590) & (t < 605)
    x[s] = 90 + 8 * (t[s] - 590)            # 90 -> 210 mmHg over 15 s, no pulse
    put(t, x, 605, 620, 300)
    acc = line_accesses(x, FS)
    assert len(acc) == 1 and abs(acc[0].start - 590) <= 2


def test_cuff_occlusion_is_not_an_access_and_short_blips_are_ignored():
    t, x = abp()
    s = (t >= 900) & (t < 925)
    x[s] = 40 + 50 * np.exp(-(t[s] - 900) / 5)  # pulse lost, pressure decays (cuff inflating upstream)
    s = (t >= 1300) & (t < 1302)
    x[s] = 80.0                                  # a 2 s flat blip
    acc = line_accesses(x, FS)
    assert [a.kind for a in acc] == ["occlusion"]
    assert match_times(acc, [1000.0])[0]["access"] is None


def test_line_taken_off_at_the_end_is_unavailable():
    t, x = abp(minutes=20)
    put(t, x, 1150, 1200, 0.0, noise=0.3)        # last 50 s at zero: the line was removed
    acc = line_accesses(x, FS)
    assert [a.kind for a in acc] == ["unavailable"]
    assert acc[0].recovered <= t[-1] + 1e-9
