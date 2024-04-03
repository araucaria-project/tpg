"""Changes the exposure times by a factor given as a parameter

1. Reads plan or master plan from `stdin`
2. In every line, searches for the exposure requests patterns `n/filter/time` e.g. `7/Ic/12`, ignoring rest of the syntax
3. Changes the exposure time by a factor given as a parameter
4. If `-f filter` is given, changes only the exposure times if the filter matches the given filter
5. Writes the plan to `stdout`, written lines are identical to the input lines but with the new exposure time

Example: Catch output from tpg, change all exptime x2, then adjust only Ic filter exptime x0.9, write to plan.txt

    tpg | chexp 2 | chexp -f Ic 0.9  > plan.txt
"""

import sys
import re
import argparse

def change_exposure_factor(plan, factor, filter=None):
    """Changes the exposure times by a factor given as a parameter

    Args:
        plan (str): Plan or master plan
        factor (float): Factor by which the exposure times will be multiplied
        filter (str): Limit changes to specific filter

    Returns:
        str: Plan with the exposure times multiplied by the factor
    """
    new_plan = []
    for line in plan.split('\n'):
        matches = re.findall(r'(\d+)/(\w+)/(\d+\.?\d*)', line)
        if matches:
            for match in matches:
                if filter and match[1] != filter:
                    continue
                old_exposure = f'{match[0]}/{match[1]}/{match[2]}'
                # clac new_exp_time - new exp time rounded to one decimal point
                new_exp_time = round(float(match[2]) * factor)
                new_exposure = f'{match[0]}/{match[1]}/{float(match[2]) * factor:.6g}'
                line = line.replace(old_exposure, new_exposure)
        new_plan.append(line)
    return '\n'.join(new_plan)

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--filter', '-f', help='Limit changes to specific filter', default=None)
    parser.add_argument('factor', type=float, help='Factor by which the exposure times will be multiplied')
    args = parser.parse_args()

    plan = sys.stdin.read()
    new_plan = change_exposure_factor(plan, args.factor, args.filter)
    sys.stdout.write(new_plan)

if __name__ == '__main__':
    main()
