#!/usr/bin/env python3

#=============================
#        Marek Gorski
#=============================

# ./tpg.py zb08 2019/08/20 [--ut --wind]  >| plan.txt

import argparse
from telescope_plan_generator import telescope_plan_generator as tpg


class main:

    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument("tel", type=str, choices=['wk06', 'zb08', 'jk15'], help='Telescope id')
    argparser.add_argument("date", nargs="+", type=str, help='date of observation')
    argparser.add_argument("--wind", type=float, default=None, help='wind direction (to avoid)')

    args = argparser.parse_args()

    p = tpg(args.tel,args.date,wind=args.wind,done_uobi=[])
    for x in p.plan:
        print(x)
####### SCRIPT RUN
    
if __name__ == "__main__":
    m = main()
