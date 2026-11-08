# LSD2

## Purpose

Least-squares molecular dating with LSD2 2.4.4, using node calibrations and/or
tip dates.

## GUI to CLI

The runner calls:

```text
lsd2 -i tree.nwk -s sequence_length -r l|a [-f 100 -q SD] -d calibrations.txt
```

| GUI control | CLI option | Notes |
|---|---|---|
| Tree file | `-i PATH` | Newick |
| Sequence length mode | none | Disabled uses 1000; Manual uses the spin value; From sequence file computes it |
| Sequence length | `-s N` | Always passed |
| Use chosen root | `-r l` | Requires the tree rooted in IcyTree |
| Estimate the root | `-r a` | Default |
| Confidence intervals | `-f 100` | 100 bootstrap trees |
| Lognormal std | `-q X` | Only with confidence intervals |
| Node calibrations | `-d FILE` | Generated from the table |
| Tip dates | `-d FILE` | Merged into the same calibration file |
| Preferred result | none | GUI-only import choice |

## Calibrations

- Node calibrations are entered in the calibration table.
- Tip dates are entered through **Configure tip dates**; they are merged into
  the same `-d` file.
- The plugin writes a platform-specific line ending (CRLF on Windows) because
  LSD2 is sensitive to it.

## Flow

1. Load a tree; optionally root it in IcyTree.
2. Set sequence length (or load a sequence file).
3. Add node calibrations and/or tip dates.
4. Enable confidence intervals if needed.
5. Run; import the dated tree.

## Pitfalls

- **Use chosen root** only takes effect if a rooted tree file is loaded;
  otherwise LSD2 is run with `-r a`.
- Sequence length is always passed. Disabled does not omit `-s`; the runner
  substitutes 1000.
- `-f` is fixed at 100 bootstrap trees in this plugin.
- `-q` is only added when confidence intervals are enabled.
- At least one calibration is normally required; the GUI can continue without
  one only after an explicit warning.
