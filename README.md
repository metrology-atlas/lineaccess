# lineaccess

Finds accesses to invasive pressure lines (blood draws, flushes, zeroing) in the pressure waveform itself.
MIT licence; needs only numpy. It is separate from physiosync, and neither depends on the other.

When an arterial (or other pressure) line is accessed, the transducer stops seeing the patient for a while:

| What happens | What the waveform shows | Kind |
|---|---|---|
| Flush from the pressure bag | about 300 mmHg, no pulse | `high` |
| Stopcock open to air, zeroing, disconnection | about 0 mmHg | `low` |
| Stopcock closed to the patient (often a slow ramp up from the flush device), clot, kink | a pressure with no pulse | `flat` |

Afterwards the waveform can stay damped for a while.

## Use

```python
from lineaccess import line_accesses, labels, match_times

acc = line_accesses(abp_mmhg, fs=500.0, t=sample_times)   # values in mmHg
labels(acc)        # [("art_line:access", start, recovered), ("art_line:unavailable", start, end), ...]
match_times(acc, lab_times)   # which access each lab sample followed, and the lag
```

- **Accesses.** `line_events` classifies short windows (1 s by default) as `high`, `low` or `flat` and merges them into events. `line_accesses` groups events less than 90 s apart into one access, meaning one visit to the line.
- **Recovery.** Each access gets a recovery period that lasts until the pulse pressure is back to 70% of its size in the minute before the access, up to 3 minutes.
- **Occlusions.** A lost pulse with no flush and no zero (only `flat`) is labelled `occlusion`, not `access`: typically a blood-pressure cuff inflating on the same arm, or a kink. Occlusions shorter than 4 s are ignored as noise.
- **Unavailable stretches.** Anything longer than 5 minutes is labelled `unavailable`: the line is not connected yet, has been removed, or the transducer is off. So is an access after which the pulse never returns before the recording ends.
- **No-pulse detection.** "No pulse" is judged after removing each 1 s window's straight-line trend, so a slow pressure ramp with no pulse (a line closed to the patient filling from the flush device) counts as `flat` and the access starts there.
- **Thresholds.** All of them live in `LineRules`.

## Design rules

- **Nothing dataset-specific.** No sampling rate, monitor, frame rate or protocol is assumed. The caller converts values to mmHg and passes sample times if the stream has breaks.
- **Labels, not deletion.** Accesses are returned as intervals so a connector can store them as annotations. Downstream users then decide what to exclude.
