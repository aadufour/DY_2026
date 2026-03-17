import pylhe

# ── 1. Load an LHE file ──────────────────────────────────────────────────────
# Replace with the path to your .lhe or .lhe.gz file
LHE_FILE = "/Users/albertodufour/MG5/mg5amcnlo/DrellYan/myLHEdir/DY_test_unweighted_events.lhe"

events = pylhe.read_lhe_with_attributes(LHE_FILE)

# ── 2. Loop over events ───────────────────────────────────────────────────────
n_events = 0
for event in events:
    n_events += 1

    # Event-level info
    print(f"\n--- Event {n_events} ---")
    print(f"  # particles : {len(event.particles)}")

    # Particle-level info
    for p in event.particles:
        print(
            f"  PID={p.id:>6}  status={p.status}  "
            f"px={p.px:.3f}  py={p.py:.3f}  pz={p.pz:.3f}  e={p.e:.3f}  m={p.m:.3f}"
        )

    if n_events >= 5:   # remove or increase to loop over more events
        break

print(f"\nDone. Printed {n_events} event(s).")