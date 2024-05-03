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
    argparser.add_argument("date", type=str, help='date of observation')
    argparser.add_argument("--ut", type=str, default=None, help="UT time of observation start")
    argparser.add_argument("--wind", type=float, default=None, help='wind direction (to avoid)')

    args = argparser.parse_args()

    flags=[]
    if args.ut == None:
        args.date = args.date + " 16:00:00"
    else:
        args.date = args.date + " " + args.ut
        flags.append("start_at_ut")

    #print(args)
    p = tpg(args.tel,args.date,wind=args.wind,flags=flags)
    for x in p.plan:
        print(x)
####### SCRIPT RUN
    
if __name__ == "__main__":
    m = main()
