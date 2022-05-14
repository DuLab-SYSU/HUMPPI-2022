import os, sys
from decimal import Decimal

def mcl(_in, _out, inflation):
	cmd = "mcl %s --abc -force-connected y -I %s -o batch_MCL_out/%s > /dev/null" % (_in, inflation, _out)
	os.system(cmd)

def main(infile, start, end):
    _in = infile
    for i in range(end):
        inflation = start  + i * Decimal('0.01')
        _out = "out.HUMPPI2022.I%s" % int(inflation * 100)
        print("Task: ", _in, _out, inflation)
        mcl(_in, _out, inflation)


if __name__ == "__main__":
    infile, start, end = sys.argv[1], Decimal(sys.argv[2]), int(sys.argv[3])
    main(infile, start, end)
