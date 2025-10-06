from __future__ import annotations

from loguru import logger
from ..formats.logo import get_logo


def write_gaussian_freq_log(vib, atoms, filename, charge=0, mult=1):
    freqs = vib.get_frequencies()
    natoms = len(atoms)
    numbers = atoms.get_atomic_numbers()
    modes = [vib.get_mode(i) for i in range(len(freqs))]
    forces = atoms.get_forces()

    # Vibration（cm-1）
    def write_modes_block(f, start, end):
        f.write("\n{:>12}".format(""))
        for i in range(start, end):
            f.write("{:^24}".format(i + 1))
        f.write("\n{:>12}".format(""))
        for _ in range(start, end):
            f.write("{:^24}".format("A"))
        f.write("\n Frequencies --")
        for i in range(start, end):
            freq = freqs[i]
            if isinstance(freq, complex):
                if abs(freq.imag) > 1e-6:
                    f.write(f"{-freq.imag:>16.4f}       ")
                else:
                    f.write(f"{freq.real:>16.4f}       ")
            else:
                f.write(f"{freq:>16.4f}       ")
        f.write("\n  Atom  AN" + "      X      Y        Z " * (end - start) + "\n")

        # GaussView Results>Vibration animation
        for a in range(natoms):
            f.write(f"{a+1:6d}{numbers[a]:4d}")
            for j in range(start, end):
                dx, dy, dz = modes[j][a]
                f.write(f"{dx:8.2f}{dy:8.2f}{dz:8.2f}")
            f.write("\n")

    # COordinates in log file
    def format_atoms(atoms, charge, mult):
        parts = [f"q\\\\t\\\\{charge},{mult}"]
        for i, atom in enumerate(atoms):
            symbol = atom.symbol
            x, y, z = atom.position
            parts.append(f"{symbol},{x:.10f},{y:.10f},{z:.10f}")
        full_string = "\\".join(parts) + "\\\\Version=Fujitsu-XTC-G16RevC.02\\"
        line_width = 70
        lines = [full_string[i:i+line_width] for i in range(0, len(full_string), line_width)]
        return "\n ".join(lines)

    # Writing to freq_from_uma.log
    with open(filename, "w") as f:
        f.write(get_logo() + "\n")
        f.write(" ----------------------------------------------------------------------\n")
        f.write(" #P \n")
        f.write(" ----------------------------------------------------------------------\n")
        f.write("\n  and normal coordinates:\n")

        for i in range(0, len(freqs), 3):
            write_modes_block(f, i, min(i + 3, len(freqs)))

        # Writing Forces
        f.write("\n ***** Axes restored to original set *****\n")
        f.write(" -------------------------------------------------------------------\n")
        f.write(" Center     Atomic                   Forces (Hartrees/Bohr)\n")
        f.write(" Number     Number              X              Y              Z\n")
        f.write(" -------------------------------------------------------------------\n")
        for i, (num, force) in enumerate(zip(numbers, forces)):
            f.write(f"{i+1:10d}{num:10d}{force[0]:16.9f}{force[1]:14.9f}{force[2]:14.9f}\n")
        f.write("-------------------------------------------------------------------\n\n\n")
        f.write(" Test job not archived.\n")
        f.write(" 1\\1\\G\\OPT\\R\\C\\S\n")
        f.write(" 5\\0\\\\#P\n ")
        formatted = format_atoms(atoms, charge, mult)
        f.write(formatted + "\n")
        f.write(" [X()]\\\\\\@\n\n Normal termination of Gaussian\n")

    logger.info(f"Log file compatible with GaussView is saved to: {filename}")
