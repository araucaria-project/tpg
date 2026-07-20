#!/usr/bin/env python3

#=============================
#        Marek Gorski
#=============================

# ./tpg.py zb08 2019/08/20  >| plan.txt

import argparse
from .telescope_plan_generator import TelescopePlanGenerator as tpg


class main:

    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument("tel", type=str, choices=['wk06', 'zb08', 'jk15', 'iris', 'test'], help='Telescope id')
    argparser.add_argument("date", nargs="+", type=str, help='yyyy/mm/dd date of observation')
    argparser.add_argument("--tag", type=str, default=None, help='tag required to allocate')
    argparser.add_argument("--wind", type=float, default=None, help='wind direction (to avoid)')
    argparser.add_argument("--fwhm", type=float, default=None, help='fwhm (to apply limits)')
    argparser.add_argument("--seed", type=int, default=None, help='seed to randomize list')
    argparser.add_argument("--loud", action="store_true", help="Enable loud mode")
    argparser.add_argument("--no_save", action="store_false", help='do not save telescope plan')
    argparser.add_argument("--no_ctc", action="store_false", help='do not use ctc')
    argparser.add_argument("--no_start_makro", action="store_false", help='do not add makro at the plan beginning')
    argparser.add_argument("--no_end_makro", action="store_false", help='do not add makro at the plan end')
    argparser.add_argument("--master_file", type=str, default=None, help='specify master file')




    args = argparser.parse_args()
    p = tpg(args.tel,args.date,
            tag=args.tag,
            wind=args.wind,
            fwhm=args.fwhm,
            loud=args.loud,
            seed=args.seed,
            save_plan=args.no_save,
            use_ctc = args.no_ctc,
            add_start_makro=args.no_start_makro,
            add_end_makro=args.no_end_makro,
            master_file=args.master_file,
            done_uobi=[])
    p.go()

####### SCRIPT RUN
    
if __name__ == "__main__":
    m = main()

