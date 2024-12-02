"""Generates in master plan format entries for specific TAB.ALL group objects

Includes information form ObjectDatabase
"""

import logging
import argparse
from pyaraucaria.lookup_objects import ObjectsDatabase

def list_groups(db, args):
    groups = db.tab_all_groups
    for group in groups:
        print(group)

def dump_group(db, group_name, args):
    grp = db.lookup_group(group_name, include_members=True)
    print(f'# Group: {group_name} [{grp["type"]}]')
    for obj in grp['objects'].values():
        print_obj_as_plan(obj, args)

def print_obj_as_plan(obj, args):
    print(f'{obj["hname"]:16} {obj["ra"]} {obj["dec"]}  seq={args.seq}  priority={args.priority}', end='')
    for label, key in [
        ('mV', 'V'),
        ('mI', 'I'),
        ('P', 'per'),
        ('hjd0', 'hjd0'),
    ]:
        try:
            print(f'  {label}={obj[key]}', end='')
        except KeyError:
            pass

    if args.comment:
        comment = args.comment
    else:
        comment = obj['group']
    print(f'  comment="{comment}"')


def main():
    # show defaults!
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # optional argument
    parser.add_argument('group', type=str, help='Group name, if empty, all groups are listed', default=None, nargs='?')
    parser.add_argument('--seq', '-s', type=str, help='Sequence, e.g 3/B/20,3/V,10', default='1/B/20,1/V/20,1/Ic/20')
    parser.add_argument('--priority', '-p', type=int, help='Priority of the group', default=1)
    parser.add_argument('--comment', '-c', type=str, help='Comment for the objects, if not provided group name is used', default=None)

    args = parser.parse_args()

    db = ObjectsDatabase()
    if args.group:
        dump_group(db, args.group, args)
    else:
        logging.warning('No group specified, listing all groups. (Use -h for help)')
        list_groups(db, args)


if __name__ == '__main__':
    main()