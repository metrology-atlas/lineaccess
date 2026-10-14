# Get started

This page finds a flush-and-zero access and a cuff occlusion in a simulated 20-minute arterial pressure
recording.

```python
import numpy as np
import lineaccess as la

fs = 500.0
t = np.arange(0, 1200, 1 / fs)                                       # 20 min at 500 Hz (s)
rng = np.random.default_rng(0)
abp = 80 + 20 * np.sin(2 * np.pi * 1.2 * t) ** 2 + rng.normal(0, 0.5, t.size)   # mmHg, 72 bpm

abp[(t >= 300) & (t < 315)] = 300.0      # flush from the pressure bag
abp[(t >= 330) & (t < 340)] = 0.0        # stopcock open to air (zeroing)
abp[(t >= 345) & (t < 360)] = 300.0      # flush again
abp[(t >= 800) & (t < 812)] = 95.0       # pulse lost with no flush or zero: a cuff on the same arm

for a in la.line_accesses(abp, fs):
    print(a.kind, a.start, a.end, a.kinds)
# access 300.0 360.0 ['high', 'low']
# occlusion 800.0 812.0 ['flat']
```

The three events at 300-360 s are less than 90 s apart, so they are one access: one visit to the line.
The pulse loss at 800 s has no flush and no zero, so it is an occlusion, not an access.

## Labels for storage

```python
print(la.labels(la.line_accesses(abp, fs)))
# [('art_line:access', 300.0, 360.0), ('art_line:occlusion', 800.0, 812.0)]
```

Each interval runs from the start of the first event to the end of recovery, when the pulse is back to
70% of its size in the minute before.

## Next steps

- Real recordings: [Find and label line accesses](how-to/find-accesses.md).
- Time blood draws for lab results: [Match lab samples to accesses](how-to/match-lab-samples.md).
