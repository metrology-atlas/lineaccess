# Match lab samples to accesses

An arterial blood sample is drawn through the line, so a lab sample usually follows an access. Matching
each lab time to the access before it estimates when the blood was actually drawn.

```python
import numpy as np
import lineaccess as la

fs = 500.0
t = np.arange(0, 1800, 1 / fs)
abp = 80 + 20 * np.sin(2 * np.pi * 1.2 * t) ** 2
abp[(t >= 300) & (t < 315)] = 300.0                       # flush after a blood draw
acc = la.line_accesses(abp, fs)

lab_times = np.array([720.0, 1700.0])                     # lab timestamps (s, same clock)
print(la.match_times(acc, lab_times, lag_range_s=(0.0, 1800.0)))
# [{'time': 720.0, 'access': 0, 'lag_s': 420.0}, {'time': 1700.0, 'access': 0, 'lag_s': 1400.0}]
```

## Choose the lag range for your lab system

The window must cover how long after the draw your lab system stamps a sample. In the VitalDB surgical
dataset, for example, the lab timestamp is a per-specimen time, not the result time: blood gases are stamped 3–12 min after the draw,
central-lab panels about 10 min or more. A 15-minute window suits blood gases; central-lab panels need
about 30 minutes.

!!! warning "A lab timestamp is not when the result was known"
    The result follows the timestamp by a turnaround the dataset may not record. Using a lab value at its
    timestamp, or at the matched draw time, in a model that runs in real time leaks information from the
    future. Use the draw time to say *what the blood reflected*, and timestamp plus turnaround to say
    *when the value could be used*.

## Check the matches against chance

`match_times` takes the latest access in the window, which favours short lags, and a long window collects
accesses that had nothing to do with the sample. Repeat the match with the lab times shifted by a random
offset of tens of minutes: the share of real times that match should be well above the share of shifted
times. On the VitalDB surgical dataset, 77% of lab groups matched an access within 15 minutes against 26% for random times.
