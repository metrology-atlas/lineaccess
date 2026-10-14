# lineaccess

**Find accesses to invasive pressure lines in the pressure waveform itself.** When an arterial (or other
pressure) line is flushed, zeroed or used for a blood draw, the transducer stops seeing the patient for a
while. lineaccess finds those stretches from the waveform alone and returns them as labelled intervals, so
downstream analyses can exclude them, or use them, for example to time blood draws.

| What happens at the line | What the waveform shows | Event kind |
|---|---|---|
| Flush from the pressure bag | about 300 mmHg, no pulse | `high` |
| Stopcock open to air, zeroing, disconnection | about 0 mmHg | `low` |
| Stopcock closed to the patient, clot, kink | a pressure with no pulse | `flat` |

## What you get

- **Accesses**: everything done in one visit to the line, with a recovery period until the pulse is back.
- **Occlusions**: a lost pulse with no flush and no zero, typically a blood-pressure cuff on the same arm.
- **Unavailable stretches**: the line not connected yet, removed, or the transducer off.
- **Matches** between accesses and external times such as lab samples.

## Design rules

- **Nothing dataset-specific.** No sampling rate, monitor or protocol is assumed. Pass values in mmHg, and
  sample times if the stream has breaks. All thresholds live in `LineRules`.
- **Labels, not deletion.** Accesses are returned as intervals for a connector to store as annotations.
  Downstream users decide what to exclude.

## Install

```bash
pip install git+https://github.com/meshlab/lineaccess      # numpy only
```

MIT licence, © 2026 Andrew J. Goodwin, University of Sydney. lineaccess is independent of physiosync;
neither depends on the other.
