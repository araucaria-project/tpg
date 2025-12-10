#!/usr/bin/env python3



import sys,numpy

class main:
    def __init__(self):
        if len(sys.argv) < 3:
            print("etc.py B,V,Ic,u,g,r,i,z 11.4")
            sys.exit(1)
        m = float(sys.argv[2])
        filters = sys.argv[1].split(",")

        m_ref = 12.
        a = 10 ** ((m-m_ref) / 2.5)

        B = 22
        V = 15
        Ic = 22
        u = 44
        g = 15
        r = 9
        i = 14
        z = 44

        txt = ""
        if "B" in filters:
            txt = txt + f"1/B/{B*a:.0f},"
        if "V" in filters:
            txt = txt + f"1/V/{V*a:.0f},"
        if "Ic" in filters:
            txt = txt + f"1/Ic/{Ic*a:.0f},"
        if "u" in filters:
            txt = txt + f"1/u/{u*a:.0f},"
        if "g" in filters:
            txt = txt + f"1/g/{g*a:.0f},"
        if "r" in filters:
            txt = txt + f"1/r/{r*a:.0f},"
        if "i" in filters:
            txt = txt + f"1/i/{i*a:.0f},"
        if "z" in filters:
            txt = txt + f"1/z/{z*a:.0f},"

        txt=txt[:-1]

        print(txt)



####### SCRIPT RUN
    
if __name__ == "__main__":
    m = main()
