# Out-of-bounds write in `cloud_generator_omp` — debugging guide for AMD

A defect in the OpenMP-offload McICA path writes one level past the end of a
cloudy column's slice of `od_scaling`. It was found on NVIDIA B200 with NVHPC
26.3, where it is fatal.

The original expectation was that the same source, compiling on amdflang with
no guards of any kind, must produce the same write on AMD. **That has since
been measured on gfx942, and the longwave write is not reachable there — see
section 11.** A genuine defect of the same class was found in the *shortwave*
solver instead: it is real on both platforms, was dormant on AMD only because
HIP hands back zeroed pages, and has now been fixed.

Sections 1 to 10 are the original diagnosis and experiment plan, preserved as
written, with corrections marked inline where measurement has since contradicted
them. Section 11 records what the measurements were and what they mean; read it
before running any of the experiments, since the most decisive of them has
already been run and two of the others will now tell you nothing.

Reference point: branch `correctness_performance_timing`, HEAD `4cba3ec`.
Line numbers in `radiation_cloud_generator_acc.F90` are exact for that commit.
Line numbers in `radiation_mcica_omp_lw.F90` are approximate because a local
debug block has been added; each snippet is quoted so you can find it by
content. Section 11's line numbers are for `53eb3fb` plus the shortwave fix.

---

## 1. The defect

`solver_mcica_omp_lw_impl` declares

```fortran
real(jprb), dimension(ng,nlev, istartcol:iendcol) :: od_scaling
```

and calls the generator once per `(jcol, jg)` from a `COLLAPSE(2)` kernel,
passing a single column's slice:

```fortran
call cloud_generator_omp(jg, ng, nlev, &
     ...
     &  od_scaling(:,:,jcol), cloud_cover_lw(jcol)+0.0_jprb, &
     &  ibegin(jcol), iend(jcol), &
     ...
```

Inside the generator the dummy is `od_scaling(ng,nlev)` and there are exactly
two writes to it:

| file:line | statement | bounded by |
|---|---|---|
| `radiation_cloud_generator_acc.F90:379` | `od_scaling(jg,jlev) = 0.0_jprb` in `do jlev = 1,nlev` | the declared extent `nlev` |
| `radiation_cloud_generator_acc.F90:477` | `od_scaling(jg,jcloud) = ...` | `iend`, a **runtime value** |

One of these writes at second index **`nlev+1`**. Because the actual argument is
a slice of a 3-D array, `od_scaling(jg, nlev+1, jcol)` is the same address as
`od_scaling(jg, 1, jcol+1)` — **the first level of the next column**.

Note what is *not* wrong: every declared bound is correct and self-consistent.
`ng` and `nlev` are passed straight through from the variables that dimension
the array. This is a bad runtime index, not a bad bound.

> **Correction (section 11).** The table above already contains the answer to
> section 6's question, and it rules out line 379. That loop runs `1,nlev` and
> the dummy's second extent *is* `nlev` — the same variable — so 379 cannot
> reach `nlev+1` on any platform or toolchain. Only line 477 can, and only if
> `iend > nlev`. Since the sanitizer nonetheless names 379, its line attribution
> in this inlined kernel is unreliable and should not be used as evidence for
> which write is guilty.

## 2. Evidence

`compute-sanitizer` on NVIDIA, at 64 columns and a single repeat:

```
========= Invalid __global__ write of size 8 bytes
=========     at radiation_cloud_generator_acc_cloud_generator_omp_+0xca0
=========        in radiation_cloud_generator_acc.F90:379
=========     Access to 0x155097f5db68 is out of bounds
=========     and is 873 bytes after the nearest allocation at 0x155097600000
=========        of size 9,820,160 bytes
=========     Device Frame: nvkernel_..._solver_mcica_omp_lw_impl__F1L416_16_
```

The arithmetic is exact and worth repeating, because it is what makes the
diagnosis solid rather than suggestive:

- The allocation, 9,820,160 bytes, is precisely `ng*nlev*ncol*8` =
  140 x 137 x 64 x 8. It is `od_scaling`, with no padding.
- The kernel launched 128 threads per block. Block 69 is the last of
  `ceil(64*140/128)` = 70, and thread 97 is collapsed iteration 8929, which
  decodes to `jcol=64, jg=110`. **Both are valid.**
- With `jg` and `jcol` fixed, the faulting address forces the second index to
  **138 = nlev+1**.

Two things this rules out. The iteration indices are correct, so it is not a
collapsed-loop or `THREAD_LIMIT` code-generation problem (independently
confirmed with a standalone reproducer). And it is an out-of-bounds access
against a live allocation, not a stale or dangling device pointer.

It reproduces with `nrepeat=1`, during the first radiation call.

## 3. Why AMD does not crash

Three effects compound:

1. **Only the last column can fault at all.** For every other column the write
   lands at `od_scaling(jg, 1, jcol+1)` — inside the array. That is silent
   corruption on every platform, NVIDIA included. The sanitizer physically
   cannot see it.
2. **The last column's overrun is `ng*8` = 1120 bytes (longwave).** NVHPC
   allocated the array byte-exactly, so it escaped immediately. HIP rounds
   allocations up to at least a 4 KB page, and large buffers more coarsely
   still, so 1120 bytes almost certainly lands in allocator slack and never
   touches an unmapped page.
3. **It may be numerically quiet.** If the culprit is line 379, it writes
   `0.0` to an address that column `jcol+1`'s own zeroing loop also sets to
   zero; it only does damage when it races that column's fill loop at 477. If
   the culprit is 477, a real scaling factor lands in the next column's top
   level, which is always wrong.

Expected AMD symptom: no crash, just small differences in top-of-atmosphere
fluxes for whichever columns lose the race. Exactly the kind of thing that sits
under a validation threshold indefinitely.

> **Correction (section 11).** This reasoning is sound but its premise did not
> hold for the longwave solver. Read this section as *why the overrun would be
> invisible on AMD if it occurred*, not as why it is invisible. Measurement
> found no longwave overrun to hide. The argument does still apply in full to
> the shortwave defect in section 11.2, which would have been silent for exactly
> these three reasons had the gate ever let it fire.

---

## 4. Experiment 1 — clamp and diff (originally "start here"; see section 11)

**Question answered: does this defect change my answers today?**

This needs no sanitizer, no instrumentation and no special build. Clamp the
only unbounded index, at `radiation_cloud_generator_acc.F90:447`:

```fortran
!            do jcloud = max(2,jlev-n_layers_to_scale),jlev-1
             do jcloud = max(2,jlev-n_layers_to_scale),min(jlev-1,nlev)
```

Then run your standard case twice, once with and once without the clamp, and
compare the output files bit-for-bit.

- **Bit-identical** → on this dataset the overrun is landing somewhere
  harmless. Still undefined behaviour and still worth fixing, but not currently
  corrupting your results.
- **Any difference** → the overrun is corrupting real data in your production
  configuration, and the size of the difference tells you how much.

Either way the clamp is a safe defensive fix, since `jcloud > nlev` has no
legitimate meaning. It is a mitigation, not a diagnosis — it does not tell you
*why* the index was out of range, which is what experiments 2 and 3 are for.

> **Result (section 11.1).** On AMD this will come back bit-identical, so it is
> no longer the place to start. `jcloud` is bounded by `iend`, and `iend` was
> measured never to exceed `nlev`, so the clamp cannot change any index. Still
> worth carrying as a permanent guard, but it will not tell you anything on this
> platform. Start at section 11 instead.

## 5. Experiment 2 — sentinel level (detects the overrun without a fault)

**Question answered: does the overrun happen on AMD, and on which columns?**

This is the portable substitute for a GPU sanitizer. Give `od_scaling` one
extra level, fill it with a value nothing else writes, and check afterwards
whether anything touched it. The extra level makes the stray write land in a
real, in-bounds location, so it both **detects** the overrun and **stops** it
corrupting the next column.

In `radiation_mcica_omp_lw.F90`, widen the declaration:

```fortran
!   real(jprb), dimension(ng,nlev, istartcol:iendcol) :: od_scaling
    real(jprb), dimension(ng,nlev+1, istartcol:iendcol) :: od_scaling
```

Nothing else needs to change: the generator's dummy stays `(ng,nlev)`, and an
actual argument larger than an explicit-shape dummy is legal. All consumers
index `1..nlev`.

Immediately **before** the kernel that calls `cloud_generator_omp`, stamp the
guard plane:

```fortran
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
    do jcol = istartcol, iendcol
      do jg = 1, ng
        od_scaling(jg,nlev+1,jcol) = -12345.0_jprb
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
```

Immediately **after** that kernel, count the survivors:

```fortran
    n_overrun = 0
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) REDUCTION(+:n_overrun)
    do jcol = istartcol, iendcol
      do jg = 1, ng
        if (od_scaling(jg,nlev+1,jcol) /= -12345.0_jprb) n_overrun = n_overrun + 1
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    write(*,*) 'od_scaling guard plane overwritten in', n_overrun, 'of', &
         ng*(iendcol-istartcol+1), 'slots'
```

`n_overrun` needs a local `integer :: n_overrun` declaration; everything else
is already in scope.

`n_overrun > 0` confirms the overrun on AMD. The count also tells you how many
columns are affected, which the NVIDIA sanitizer could not — it only ever saw
the last one. Costs `ng*ncol*8` bytes, about 19 MB at nproma=16960.

If you want the column list rather than a count, replace the reduction with a
per-column flag array and print the indices.

## 6. Experiment 3 — which of the two writes is it

**Question answered: line 379 or line 477?**

Line 477's `jcloud` is bounded by `iend`, so measure `iend` directly. Insert
this between the kernel that computes `ibegin`/`iend` (the one containing
`ibegin(jcol) = nlev`) and the kernel that calls the generator:

```fortran
    !$OMP TARGET UPDATE FROM(ibegin, iend, cloud_cover_lw)
    do jcol = istartcol, iendcol
      if (cloud_cover_lw(jcol) >= cloud_fraction_threshold) then
        if (iend(jcol) > nlev .or. ibegin(jcol) < 1 .or. iend(jcol) < ibegin(jcol)) then
          write(*,'(a,i8,a,i8,a,i8,a,i8)') 'BAD RANGE jcol=', jcol, &
               ' ibegin=', ibegin(jcol), ' iend=', iend(jcol), ' nlev=', nlev
        end if
      end if
    end do
    write(*,'(a,i8,a,i8,a,i8)') 'cloud range: max iend=', &
         maxval(iend(istartcol:iendcol)), ' min ibegin=', &
         minval(ibegin(istartcol:iendcol)), ' nlev=', nlev
```

Reading it:

| result | conclusion |
|---|---|
| `BAD RANGE` lines appear | the culprit is **477**; fix how `ibegin`/`iend` are produced |
| clean, `max iend` <= nlev | `iend` is exonerated; the culprit is **379**, meaning the `nlev` reaching the generator disagrees with the slice's actual second extent |

The per-column test is guarded on `cloud_cover_lw` because `ibegin`/`iend` are
*legitimately* uninitialised for clear columns — see section 8. The summary
line deliberately drops that guard, so if the unguarded `max iend` is wild while
every guarded column looks clean, the two guards are not agreeing the way the
source suggests.

> **Result (section 11.1).** This experiment has been run on gfx942, in a
> stronger device-side form that also distinguishes *unset* from merely
> out-of-range. The longwave came back clean — no `BAD RANGE` lines, and
> `max iend` exactly equal to `nlev` across 16960 columns at five block sizes
> and on a second input. Per the correction in section 1, the second row of the
> table below is not a valid conclusion to draw from that: a clean `iend` does
> not implicate 379, because 379 cannot be guilty. It exonerates the longwave
> solver entirely.

On NVIDIA this build must `stop` right after the print, because otherwise the
faulting kernel launches and costs a GPU. On AMD there is no fault, so you can
let it run on.

## 7. Experiment 4 — ROCm address sanitizer

The direct analogue of the `compute-sanitizer` run that found this. Roughly:

```bash
# target must be an xnack+ variant, e.g. gfx942:xnack+
amdflang ... -fsanitize=address -shared-libasan -g
HSA_XNACK=1 ./ecrad_ifs_blocked_dp config.nam input.nc out.nc
```

Verify the exact flags against your ROCm version — the device ASAN invocation
has changed across releases, and it needs an xnack-capable target plus
`HSA_XNACK=1` at runtime. Expect it to be slow; use the smallest column count
that still contains cloudy profiles.

Be aware of the section 3 limitation: **ASAN will only flag the last column**,
for the same reason `compute-sanitizer` did, and only if the overrun clears the
allocator's rounding. A clean ASAN run does **not** mean the overrun is absent.
Experiment 2 is the more reliable detector on AMD; treat ASAN as confirmation.

---

## 8. Where to look for the cause

If experiment 3 implicates `iend`, the most likely mechanism is that it is read
uninitialised. `ibegin` and `iend` are assigned only inside a guard:

```fortran
    do jcol = istartcol,iendcol
      if (cloud_cover_lw(jcol) >= cloud_fraction_threshold) then
        ibegin(jcol) = nlev
        do jlev = 1, nlev
          if( frac(jlev,jcol) > 0.0_jprb ) then
            ibegin(jcol) = min(jlev, ibegin(jcol))
          end if
        end do

        iend(jcol) = ibegin(jcol)
        do jlev = ibegin(jcol)+1,nlev
          if (frac(jlev,jcol) > 0.0_jprb) then
            iend(jcol) = max(jlev, iend(jcol))
          end if
        end do
      end if
    end do
```

They are mapped with `MAP(ALLOC:)` and never otherwise initialised, so any
column missing that guard leaves whatever was in device memory. The generator
guards on the same condition (`total_cloud_cover >= frac_threshold`), so in
principle the two always agree and the garbage is never read — but that is an
argument from reading the source, and a stale `nlev+1` in that buffer would
reproduce the observed fault exactly.

Analytically `iend` cannot otherwise exceed `nlev`: it is a `max` over `jlev`
in `[ibegin+1, nlev]`, starting from `ibegin <= nlev`.

> **Correction.** For the *longwave* solver the two guards cannot disagree, so
> this mechanism is unreachable and the caution above was unnecessary. The
> reason is two lines earlier than the snippet quoted, at
> `radiation_mcica_omp_lw.F90:367`:
>
> ```fortran
>       cloud_cover_lw(jcol) = cum_cloud_cover(nlev,jcol)
>       if (cloud_cover_lw(jcol) < cloud_fraction_threshold) then
>         cloud_cover_lw(jcol) = 0.0_jprb
>       end if
> ```
>
> That runs for **every** column, unconditionally, and clamps sub-threshold
> covers to exactly `0.0`. So a column that skips the `ibegin`/`iend` loop holds
> `0.0`, the generator receives that same value, tests it against the same
> `cloud_fraction_threshold` (default `1.0e-6`), and returns early without
> touching `od_scaling`. Confirmed by measurement in section 11.1.
>
> Two other routes to a stale `iend` were checked and are also closed: the
> `MAP(ALLOC:)` at `:290` has a matching `MAP(DELETE:)` at `:708`, so no map
> entry survives between calls to be re-read over reused stack storage; and the
> range-finding kernel is a plain `jcol` loop with no race.
>
> The *shortwave* solver is a different story: there the guards genuinely did
> not agree. See section 11.2.

Robust fix if this is confirmed — give both arrays unconditional defaults so no
column can reach the generator with garbage, regardless of whether the guards
agree:

```fortran
    do jcol = istartcol,iendcol
      ibegin(jcol) = 1
      iend(jcol)   = 0        ! empty range: generator's jcloud loop does not execute
      if (cloud_cover_lw(jcol) >= cloud_fraction_threshold) then
        ...
```

Keep the experiment 1 clamp as well. Defaults address the cause; the clamp
guarantees the array cannot be indexed out of range whatever the cause.

## 9. Also check the shortwave solver

`radiation_mcica_omp_sw.F90` has identical exposure — same generator, same
call shape:

```fortran
! radiation_mcica_omp_sw.F90:226
real(jprb), dimension(ng,nlev,istartcol:iendcol) :: od_scaling
! radiation_mcica_omp_sw.F90:389 / :397
call cloud_generator_omp(jg, ng, nlev, &
     ...
     &  od_scaling(:,:,jcol), cloud_cover_sw(jcol)+0.0_jprb, &
```

Apply every experiment here too. The longwave is only where it happened to be
caught.

> **Correction.** This instinct was right and the priority was wrong — the
> shortwave exposure is not "identical", it is strictly worse. Its
> range-finding gate carries an extra `cos_sza_col(jcol) > 0.0` condition that
> the generator's own gate does not, and `cloud_cover_sw` is itself only
> assigned for daylight columns. That is a real guard mismatch of exactly the
> kind section 8 hypothesised for the longwave. It is the one defect of this
> class that measurement confirmed to be present in the AMD build. Section 11.2
> has the mechanism and the fix.

## 10. Context worth having

- `radiation_cloud_generator_acc.F90` contains **no preprocessor guards at
  all**, so amdflang and nvfortran compile identical source. Nothing about this
  defect is NVIDIA-specific.
  > **Correction.** The first half is true and worth keeping; the conclusion
  > does not follow. Identical source guarantees identical *semantics*, not
  > identical *behaviour*, once the code reads memory nobody wrote — which is
  > exactly what section 11.2 found. Uninitialised device memory is where two
  > platforms compiling the same text legitimately diverge, and it is why the
  > shortwave defect is live on both while being observable on only one.
- The **OpenACC** path does not have it. `radiation_mcica_acc_lw.F90` calls
  `cloud_generator_acc(ng, nlev, ...)` with **no `jg` argument**: the g-point
  loop is internal and `od_scaling` is a per-column 2-D local. The OpenMP port
  hoisted `jg` into a kernel index and passes a 3-D slice instead. The two
  generators are otherwise line-for-line identical, so the defect is in that
  restructuring.
- Only two commits have ever touched `radiation_mcica_omp_lw.F90`. The first,
  `4778013`, states it *"compiles locally for amdflang and nvfortran"* —
  compiles, not runs. The only correctness work since, `4cba3ec`, is *"Fix
  McICA partial-block shapes on Flang"*, which suggests shape trouble has
  already been hit in this solver and may be a related symptom patched at a
  different layer. Worth reading before starting.

---

## 11. Measured on AMD

Hardware gfx942, AFAR 24.1.0-pre, `ecrad_ifs_blocked_dp`, tree at `53eb3fb`.
Experiment 3 has been run, in a stronger form than section 6 describes. Its
result also settles what experiment 2 would have detected, since `jcloud` is the
only unbounded index and it is bounded by `iend`, so bounding `iend` bounds the
write — though note that is an inference, not the direct detection experiment 2
would give you.

The headline is that the longwave overrun this document was written about does
not occur on AMD and cannot, while a real defect of the same class was found in
the shortwave solver and fixed.

### 11.1 The longwave solver is clean

Experiment 3 was run in a stronger device-side form. Rather than a
`TARGET UPDATE FROM` and a host-side scan, `ibegin` and `iend` were poisoned to
`-999999` on the device immediately after their `MAP(ALLOC:)`, and a reduction
kernel placed immediately before the generator counted, among exactly the
columns that pass the generator's *own* gate, how many held the poison value
(*unset*), how many held an out-of-range value, and what the maximum `iend`
was. This distinguishes "never written" from "written wrongly", which the
host-side version cannot, and it measures precisely the population the
generator will touch.

Production case, 16960 columns, at block sizes 16960, 8192, 4000, 1024 and 137:

| solver | gate passes | unset | out of range | max `iend` |
|---|---|---|---|---|
| longwave | 16238 of 16960 | 0 | 0 | 137 (= `nlev`) |
| shortwave | 16238 of 16960 | 0 | 0 | 137 (= `nlev`) |

`max iend` lands exactly on `nlev` and never past it. Repeated on
`test/ifs/ecrad_meridian.nc`, which unlike the production input does contain
night columns: same result.

Combined with the source argument, this closes the longwave question. The
second index of `od_scaling` is bounded by `nlev` at line 379 and by `iend` at
line 477, since `jcloud <= jlev-1` and `jlev <= iend+1`; `iend <= nlev`
whenever it is written; and per the correction in section 8 it is always
written for every column the generator will actually process. There is no path
to `nlev+1`.

An independent check of the section 2 address arithmetic supports the decode but
not the attribution. `0x155097f5db68 - 0x155097600000` is 9,821,032, which is
872 bytes past the 9,820,160-byte allocation, or exactly 109 elements — giving
`od_scaling(110, 1, 65)` in a 64-column array, the same address as
`od_scaling(110, 138, 64)`. So `jg=110`, `jcol=64`, second index `nlev+1`, as
stated. But 379, which the sanitizer names, provably cannot generate that index.

### 11.2 The shortwave solver had a real guard mismatch

Two conditions gate the range-finding loop, at
`radiation_mcica_omp_sw.F90:373`:

```fortran
      if (cos_sza_col(jcol) > 0.0_jprb .and. cloud_cover_sw(jcol) >= cloud_fraction_threshold) then
```

The generator is then called unconditionally for every column and gates on
cloud cover alone. Worse, `cloud_cover_sw` is *itself* only assigned inside
`if (cos_sza_col(jcol) > 0.0_jprb)`, at `:351`. So for a night-time column
neither the gate value nor the range is ever written, and both are read from
`MAP(ALLOC:)` device memory. A stale gate value reaching `1.0e-6` would send the
generator into that column with an unwritten `iend`, and `iend > nlev` writes at
exactly the `od_scaling(jg, nlev+1, jcol)` address seen on the B200.

This is the mechanism section 8 proposed, in the solver section 9 deprioritised.

Unlike the longwave case it is not merely hypothetical, but on AMD it is
dormant. Instrumenting the four night columns of the meridian case to report
what that buffer actually contains:

```
IEND CHECK SW: cols=8 night-cols=4 gate-passes=3 unset=0 out-of-range=0 max-iend=137
   night-col cloud_cover_sw as read from device: min= 0.00000E+00 max= 0.00000E+00  gate threshold= 1.00000E-06
```

Exactly zero, against a threshold of `1.0e-6`, so the gate rejects those
columns and the unwritten range is never read. Note the corroborating detail in
the same block: the longwave gate passes 7 columns and the shortwave gate only
3. Those 4 excluded columns are the night ones, and their exclusion rested
entirely on that zero. Nothing in the source, in OpenMP, or in HIP guarantees a
zeroed allocation — this was allocator behaviour, not correctness.

### 11.3 Fix applied

An `else` branch at `radiation_mcica_omp_sw.F90:365` writes the night-time cover
explicitly:

```fortran
      else
        cloud_cover_sw(jcol) = 0.0_jprb
      end if
```

The gate outcome for a night column is now determined by the source rather than
by what the device buffer happened to hold, which keeps such a column out of the
generator and makes the unwritten `ibegin`/`iend` unreachable.

Validation on the meridian case, the one with night columns, is **bit-identical
to the pristine baseline** across all 13 output variables, with `nccmp` exiting
0. That is the expected outcome: it pins down behaviour that was previously
correct only by luck. On the production case the change is a no-op, because all
16960 columns have `cos_sza` between 0.025 and 0.081 — every one is daylight,
which is also why that input can never exercise this path. Throughput unchanged
at roughly 90k columns/s.

One residual: the fix relies on `cloud_fraction_threshold` being strictly
positive, which it is by default at `1.0e-6`. Were it ever set to zero, a night
column at cover `0.0` would satisfy `>= threshold` and enter the generator
again. The section 8 remedy of giving `ibegin`/`iend` unconditional defaults
would close that too, and remains worth doing as belt and braces.

Unrelated but worth recording so it is not mistaken for a regression: the
production case fails validation against `ecrad_ifs_blocked_dp_A100-master-acc.nc`
on `flux_net_lw` (0.001236) and `flux_net_lw_clear` (0.0012207) against a 0.001
threshold. A pristine rebuild reproduces both figures identically, so this
predates the fix.

### 11.4 What this leaves open on NVIDIA

The shortwave gate mismatch is now the better suspect for what
`compute-sanitizer` caught. NVIDIA's allocator does not zero, so a non-zero
stale gate value is plausible there, and the sanitizer's attribution was
already shown to be unreliable — it names line 379, which is provably in
bounds, so its naming of the longwave kernel in inlined code deserves the same
scepticism.

The cheap test is to re-run the B200 sanitizer with `do_sw=false`. If the fault
disappears, it was the shortwave path, and the fix in 11.3 addresses it. If it
survives, the longwave fault is real on NVHPC while being unreachable in the
same source on amdflang, which would make it a code-generation problem rather
than a source defect and would need a different investigation.
