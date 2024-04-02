"""Changes the exposure times by a factor given as a parameter

1. Reads plan or master plan from stdin
2. In every line, searches for the exposure requests patterns `n/filter/time` e.g. `7/Ic/12`, ignoring rest of the syntax
3. Changes the exposure time by a factor given as a parameter
4. Writes the plan to stdout, written lines are identical to the input lines but with the new exposure time
"""

import sys
import re
import argparse

def change_exposure_factor(plan, factor):
    """Changes the exposure times by a factor given as a parameter

    Args:
        plan (str): Plan or master plan
        factor (float): Factor by which the exposure times will be multiplied

    Returns:
        str: Plan with the exposure times multiplied by the factor
    """
    new_plan = []
    for line in plan.split('\n'):
        matches = re.findall(r'(\d+)/(\w+)/(\d+\.?\d*)', line)
        if matches:
            for match in matches:
                old_exposure = f'{match[0]}/{match[1]}/{match[2]}'
                new_exposure = f'{match[0]}/{match[1]}/{float(match[2]) * factor:.6g}'
                line = line.replace(old_exposure, new_exposure)
        new_plan.append(line)
    return '\n'.join(new_plan)

def main():
    parser = argparse.ArgumentParser(description='Changes the exposure times by a factor given as a parameter')
    parser.add_argument('factor', type=float, help='Factor by which the exposure times will be multiplied')
    args = parser.parse_args()

    plan = sys.stdin.read()
    new_plan = change_exposure_factor(plan, args.factor)
    sys.stdout.write(new_plan)

if __name__ == '__main__':
    main()
