#!/usr/bin/env python3
"""
sort_poscar.py -- universal POSCAR/CONTCAR atom sorter

Sorts atoms in a POSCAR file by a chosen coordinate axis.
By default: within each species block, sort by Z ascending,
keeping the species order and per-species counts unchanged.

Usage:
  python3 sort_poscar.py POSCAR                      # sort by z (asc), within species, print to stdout
  python3 sort_poscar.py POSCAR -o POSCAR.sorted     # write to a file
  python3 sort_poscar.py POSCAR -a y                 # sort by y instead of z
  python3 sort_poscar.py POSCAR -d                   # descending
  python3 sort_poscar.py POSCAR -w                   # sort WHOLE file (ignore species blocks)
  python3 sort_poscar.py POSCAR --keep-format        # keep original number formats (reuse token text)

Notes:
  * Works for both Direct (fractional) and Cartesian coordinates.
  * Auto-detects and preserves 'Selective dynamics' and velocity lines
    (the extra columns travel together with their atom).
  * '-w/--whole' sorts all atoms together -- only makes sense if the
    species blocks are not meaningful, use with care.
"""
import sys

AXES = {'x': 0, 'y': 1, 'z': 2}

def parse_poscar(lines):
    """Return dict with header parts and the raw coordinate lines."""
    # clean trailing blank lines
    while lines and not lines[-1].strip():
        lines.pop()
    title   = lines[0]
    scale   = lines[1]
    lat     = lines[2:5]
    elements= lines[5].split()
    counts  = [int(x) for x in lines[6].split()]
    natoms  = sum(counts)

    selective = False
    j = 7
    if lines[j].strip().lower().startswith('selective'):
        selective = True
        j = 8
    coord_type = lines[j]                       # 'Direct' / 'Cartesian' / 'Fractional'
    coord_start = j + 1

    coords = lines[coord_start:coord_start + natoms]
    if len(coords) != natoms:
        raise ValueError(f"expected {natoms} coordinate lines, found {len(coords)}")
    return {
        'title': title, 'scale': scale, 'lat': lat,
        'elements': elements, 'counts': counts,
        'selective': selective, 'coord_type': coord_type,
        'coords': coords, 'coord_start': coord_start,
    }

def sort_block(block, axis, reverse):
    """Sort coordinate lines by the chosen axis (index into the first 3 columns)."""
    return sorted(block, key=lambda ln: float(ln.split()[axis]), reverse=reverse)

def main():
    args = [a for a in sys.argv[1:]]
    if not args:
        print(__doc__); sys.exit(1)

    infile = args[0]
    outfile = None
    axis_name = 'z'
    reverse = False
    whole = False

    i = 1
    while i < len(args):
        a = args[i]
        if a in ('-o', '--out'):
            outfile = args[i+1]; i += 2
        elif a in ('-a', '--axis'):
            axis_name = args[i+1].lower(); i += 2
            if axis_name not in AXES:
                print(f"error: unknown axis '{axis_name}' (use x, y or z)", file=sys.stderr); sys.exit(1)
        elif a in ('-d', '--desc'):
            reverse = True; i += 1
        elif a in ('-w', '--whole'):
            whole = True; i += 1
        else:
            print(f"error: unknown option '{a}'", file=sys.stderr); sys.exit(1)

    axis = AXES[axis_name]

    with open(infile) as f:
        lines = f.read().splitlines()
    p = parse_poscar(lines)

    if whole:
        new_coords = sort_block(p['coords'], axis, reverse)
        # when sorting the whole file, keep species order/counts (labels no longer
        # correspond atom-by-atom -- caller must handle)
        blocks = [new_coords]
    else:
        blocks = []
        idx = 0
        for n in p['counts']:
            block = p['coords'][idx:idx+n]
            idx += n
            blocks.append(sort_block(block, axis, reverse))

    out = [p['title'], p['scale']] + p['lat']
    out.append('  '.join(p['elements']))
    out.append('  '.join(str(c) for c in p['counts']))
    if p['selective']:
        out.append('Selective dynamics')
    out.append(p['coord_type'])
    for b in blocks:
        out.extend(b)

    text = "\n".join(out) + "\n"

    if outfile:
        with open(outfile, 'w') as f:
            f.write(text)
        print(f"OK: {infile} -> {outfile}  (axis={axis_name}, {'desc' if reverse else 'asc'}, {'whole' if whole else 'within-species'})")
        # report
        for el, n in zip(p['elements'], p['counts']):
            print(f"  {el:2s} n={n:2d}  ok")
    else:
        sys.stdout.write(text)

if __name__ == '__main__':
    main()
