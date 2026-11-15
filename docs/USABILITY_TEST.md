# Usability test script (U7)

**Audience:** Windows user, no terminal experience.  
**Time limit:** 30 minutes, unaided.  
**Recorder:** note every hesitation, mis-click, and error message.

## Setup

1. Fresh clone of the repository.
2. Docker Desktop installed but not pre-configured.
3. Optional: one real chr20 CRAM on the desktop (for step 3).

## Tasks

| # | Task | Success criteria |
|:--|:--|:--|
| 1 | Install | Docker Desktop running; double-click `Start ARC-SV.bat` opens browser at Studio |
| 2 | Example | **Try example** → install demo → **Results** shows calls and chart |
| 3 | Real data | Copy CRAM to chosen data folder → **New analysis** → Quick test (chr20) → **My analyses** shows done |
| 4 | Export | **Results** → download VCF or `.tab` |

## Log template

```
Tester ID:
Date:
Step 1 notes:
Step 2 notes:
Step 3 notes:
Step 4 notes:
Blockers (file + screenshot path):
Suggested doc/UI fixes:
```

After the session, file issues or update `TROUBLESHOOTING.md` and Studio copy from repeated blockers.
