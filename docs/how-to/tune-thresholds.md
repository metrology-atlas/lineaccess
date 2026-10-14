# Tune the thresholds

All thresholds are fields of `LineRules`. The defaults suit adult arterial lines with a pressure bag at
about 300 mmHg.

| Field | Default | Change it when |
|---|---|---|
| `high_mmhg` | 200 | the flush pressure is lower (paediatric bags, other lines) |
| `low_mmhg` | 10 | the zero reads off zero by more than a few mmHg |
| `flat_pp_mmhg` | 4 | the pulse is very small (low output) and is being read as no pulse |
| `group_s` | 90 | your sampling protocol leaves longer gaps between flush and zero |
| `max_access_s` | 300 | long procedures at the line should still count as one access |
| `recovery_frac` | 0.7 | recovery should mean a fuller return of the pulse |
| `min_occlusion_s` | 4 | short pulse losses are noise in your data, or real and worth keeping |

```python
import lineaccess as la

paediatric = la.LineRules(high_mmhg=150.0, flat_pp_mmhg=2.0)
```

Check a change on a sample of plotted accesses (90 s of waveform either side) before using it for a
whole dataset.
