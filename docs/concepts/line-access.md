# What a line access looks like

An invasive pressure line connects the patient's artery (or vein) through fluid-filled tubing to a
transducer. A three-way stopcock near the transducer lets staff draw blood, flush the line from a
pressurised bag, or open the transducer to air to set zero. During each of these the transducer measures
something other than the patient.

| Step | Transducer sees | Waveform | Kind |
|---|---|---|---|
| Flush | the pressure bag | flat at about 300 mmHg | `high` |
| Open to air, zero | the atmosphere | flat at about 0 mmHg | `low` |
| Closed to the patient | the tubing, slowly filling from the flush device | a pressure with no pulse, often a slow ramp | `flat` |
| Back to the patient | the artery, through tubing that may hold bubbles or clot | a pulse, often damped at first | recovery |

A blood draw typically looks like a short `flat` or `low` stretch (closed to the patient while the sample
is taken) followed by one or two flushes. Zeroing adds a `low` plateau. lineaccess groups everything
within 90 s into one access and then waits for the pulse to recover.

## Why "no pulse" is judged after detrending

With the stopcock closed to the patient, the pressure in the tubing often ramps slowly as the flush
device fills it. Its plain spread can look like a pulse pressure; after removing each 1 s window's
straight-line trend, almost nothing is left. That stretch is `flat`, and the access starts there.

## Occlusions and unavailable stretches

A lost pulse with no flush and no zero is usually not an access: a blood-pressure cuff inflating on the
same arm, or a kinked line. It is labelled `occlusion`. Anything longer than 5 minutes is `unavailable`:
the line is not connected yet, has been removed, or the transducer is off.

## Effect on derived values

During an access and its recovery the waveform does not describe the patient: pressures, pulse
timing and anything derived from them should be excluded or flagged. The arterial upstroke time, for
example, is not physiological during a flush. After zeroing, the pressure level can shift slightly if the
transducer height changed.
