# Find and label line accesses

## 1. Convert to mmHg and keep the sample times

lineaccess works in mmHg. If your recording stores integer counts, apply its gain and offset first. If
the stream has breaks, pass the time of every sample (`t`), so windows across a break are not merged;
otherwise pass the first sample's time as `t0`.

```python
import numpy as np
import lineaccess as la

fs = 500.0
t = np.r_[np.arange(0, 600, 1 / fs), np.arange(900, 1500, 1 / fs)]   # a 5-minute gap in the stream
abp = 80 + 20 * np.sin(2 * np.pi * 1.2 * t) ** 2
abp[(t >= 1000) & (t < 1012)] = 300.0                                  # a flush
acc = la.line_accesses(abp, fs, t=t)
print([(a.kind, round(a.start, 1), round(a.end, 1)) for a in acc])     # [('access', 1000.0, 1012.0)]
```

## 2. Read the result

Each `LineAccess` has a `kind` (`access`, `occlusion` or `unavailable`), its `start` and `end`, the time
its pulse `recovered`, the event `kinds` it contained and the highest event median (`max_mmhg`, e.g. the
flush pressure). `recovered_ok` is false when the pulse did not come back within 3 minutes.

```python
for a in acc:
    print(a.to_dict())
```

## 3. Store as labels

```python
rows = la.labels(acc, prefix="art_line")      # ("art_line:access", start, recovered), ...
```

For a central venous line, use another prefix, e.g. `prefix="cvp_line"`.

## 4. Reuse events

`line_events` gives the raw events (one per run of `high`, `low` or `flat` windows). Pass them to
`line_accesses(..., events=...)` to group them without scanning the waveform again, for example when trying
different grouping thresholds.

```python
ev = la.line_events(abp, fs, t=t)
acc_long = la.line_accesses(abp, fs, t=t, events=ev, rules=la.LineRules(group_s=180.0))
```
