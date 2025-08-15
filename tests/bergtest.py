#!/usr/bin/env python
# coding: utf-8

from __future__ import division, print_function
from decimal import Decimal
from decimal import getcontext
#from __future__ import unicode_literals
import subprocess
import os
import os.path
import shlex
import time
import datetime
import sys

from pprint import pprint
from ast import literal_eval
import yaml


class wcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''


"""bergtest. Provisional notebook/script for executing the function tests.
"""


def similar(outfilevalue, expectedvalue, relativeprecision=1E-6,
            absoluteprecision=1e-8):
    def levalfilter(s):
        filtered = s if isinstance(s, (int, float)) else literal_eval(s)
        return filtered
    raws = (outfilevalue, expectedvalue, relativeprecision, absoluteprecision)
    ov, ev, relprec, absprec = map(levalfilter, raws)

    # and abs((ov-ev) / (ev+1.0e-14)) < abs(relprec)
    return abs(ov - ev) <= abs(absprec)


def almost(outfilevalue, expectedvalue, relativeprecision=5E-4,
           absoluteprecision=5e-4):
    def levalfilter(s):
        filtered = s if isinstance(s, (int, float)) else literal_eval(s)
        return filtered
    raws = (outfilevalue, expectedvalue, relativeprecision, absoluteprecision)
    ov, ev, relprec, absprec = map(levalfilter, raws)

    return abs(ov - ev) <= abs(absprec)  # and abs((ov-ev) / ev) < abs(relprec)


def sloppy(outfilevalue, expectedvalue, relativeprecision=2E-2,
           absoluteprecision=2e-2):
    def levalfilter(s):
        filtered = s if isinstance(s, (int, float)) else literal_eval(s)
        return filtered
    raws = (outfilevalue, expectedvalue, relativeprecision, absoluteprecision)
    ov, ev, relprec, absprec = map(levalfilter, raws)

    return abs(ov - ev) <= abs(absprec) and abs((ov-ev) / ev) < abs(relprec)


def only_abs(outfilevalue, expectedvalue, relativeprecision=1E-3,
             absoluteprecision=1e-6):
    def levalfilter(s):
        filtered = s if isinstance(s, (int, float)) else literal_eval(s)
        return filtered
    raws = (outfilevalue, expectedvalue, relativeprecision, absoluteprecision)
    ov, ev, relprec, absprec = map(levalfilter, raws)

    return abs(ov - ev) <= abs(absprec)


def compare(outfilevalue, expectedvalue, comparison_func=similar, **kwargs):
    """
    compare is essentially a wrapper for comparison_func to keep naming
    consistent internally. **kwargs arepassed on to comparison_func.
    comparison_func is expected to take two mandatory (usually numerical)
    arguments to be compared
    """
    return comparison_func(outfilevalue, expectedvalue, **kwargs)


def readconfig(definitionfilename):
    with open(definitionfilename, 'r') as yf:
        yamlob = yaml.load(yf, Loader=yaml.SafeLoader)
    return yamlob


def extract(test, outfile, headers=None, skiprows=None):
    import extractoutput
    # reload(extractoutput)
    filehandler = extractoutput.csv_like
    return filehandler(
        outfile, test['select'], test['extract'],
        headers=headers, verbose=0, skiprows=skiprows)


def runexternal(runcmd, wd):

    # running the SD command as a subprocess
    # Extra effort in this case, but we want to be able to accept from command
    # line later
    popenargs = shlex.split(runcmd)
    sp = subprocess.Popen(
        args=popenargs,
        cwd=wd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = com = sp.communicate()
    # print the output if needed
    # for elem in com: pprint(elem.split("\n"))
    return com


def extractandcompare(test, outfile, headers=None, skiprows=None):
    # Comparisons between output files and expected values

    outfilevalue = extract(test, outfile, headers=headers, skiprows=skiprows)

    if outfilevalue:
        comcrit = test.get('comparison_crit') if test.get(
            'comparison_crit') else {}

        # trusting here, against better knowledge, that value at 'expected'
        # will always correspond to the shape of the extractor.
        # expected = test['expected'] if isinstance(test['expected'], list)
        # else [test['expected']]
        # Convert list of expected values to floats. Take care if other types
        # are compared.
        # print test['expected']
        expected = [float(i) for i in test['expected']] if isinstance(
            test['expected'], list) else [test['expected']]
        # expected = [ Decimal(i) for i in test['expected'] ]
        # if isinstance(test['expected'], list) else [test['expected']]
        # expected = map(literal_eval, expected)
        # # Evaluating the expected value would in a sense make it more general
        # (semantic value rather than string pattern matching), but doing so
        # would force putting the mathematical comparison as early as at the
        # extraction stage. Do we want to do that?

        # Uncomment to be able to extend by eval:ing the comparison_func key
        comfun = eval(test.get('comparison_func')) if test.get(
            'comparison_func') else similar
        return (
            map(
                lambda outfilevalue, expected: compare(
                    outfilevalue,
                    expected,
                    comparison_func=comfun,
                    **comcrit
                ),
                outfilevalue, expected
            ),
            outfilevalue)
    else:
        return None


def execcase(case, externallabel=None, reallyrun=True):
    externaltestname = externallabel if externallabel else ""
    internaltestname = case.get('name') if case.get('name') else ""

    olddir = os.getcwd()
    wd = infiledir = os.path.abspath(case.get('workingdir'))
    if os.path.isdir(wd):
        os.chdir(wd)
    else:
        print("In testcase %s %s: %s is not a valid directory " % (
            externaltestname, ", %s" % internaltestname
            if internaltestname else "", wd))

    try:
        path_to_executable = os.path.abspath(case['run'])
    except KeyError as kerr:
        path_to_executable = None
        print("Found no valid run command for case %s,%s in definition file %s. Proceeding under the assumption that outputfiles have already been created." % (
            externaltestname, internaltestname, definitionfilename))

    if path_to_executable and reallyrun:
        # may wish to add args here for some cases
        runcmd = (path_to_executable)
        out, err = runexternal(runcmd, wd)
        if err:
            print("Something fishy")
            map(pprint, ("===== Error(s) occured: =====", err.split("\n"), "-"*25))

    ccoms = case['comparisons']

    for ccom in ccoms:
        for test in ccom['testlist']:
            if os.path.isfile(ccom['outfile']):
                test['cmpoutcome'], test['extracted'] = extractandcompare(
                    test,
                    ccom['outfile'],
                    headers=ccom.get('headers'), skiprows=ccom.get('skiprows'))
            else:
                test['cmpoutcome'], test['extracted'] = (
                    [False], "NO_SUCH_FILE_ERROR")

    os.chdir(olddir)


#----------------------------------------------#
# Get all input arguments and put them to a list
#
# Thomas Nystrand 2014/08/07:
# Added --case, --dirty,
# Changed descriptorfiles to --files
#----------------------------------------------#
def getargs():
    #
    import argparse
    #
    parser = argparse.ArgumentParser(
        description='Check output files against expected values.')
    parser.add_argument(
        '--files',  help='File(s) specifying the tests',
        default=['regulartests.yaml'], dest='descriptorfiles', metavar="FILE",
        type=str, nargs='*')
    parser.add_argument(
        '--cleanup', help='Clean up test files without continue execution',
        default=False, action="store_true")
    parser.add_argument(
        '--dryrun',
        help='Dry run, i.e. do not run simulation code before testing',
        default=False, action="store_true")
    parser.add_argument(
        '--binary', help='Full path to the execeutable to be tested',
        metavar='binary', default='../../source/sd')
    parser.add_argument(
        '--dirty',  help='Avoid automatic cleanup before simulation',
        default=False, action="store_true")
    parser.add_argument(
        '--case',   help='Run specific simulation',
        dest='selectedcase', type=int, default=-1)
    parsed = parser.parse_args()

    os.environ["SD_BINARY"] = parsed.binary

    if(parsed.dryrun):
        parsed.reallyrun = 0
    else:
        parsed.reallyrun = 1
    return parsed


def printfun(*args):
    for a in args:
        print(a)
# pprint(args)


def presentresults(case):

    numtests = 0
    numfailed = 0

    for ccom in case.get('comparisons'):
        #printfun( "ccom ", ccom )

        for test in ccom.get('testlist'):
            numtests += 1
            ##printfun( "test ", test)
            if test.get('type'):
                print('%30s' % (wcolors.OKBLUE+" " +
                                test.get('type')+wcolors.ENDC+":"), end='')
            if test.get('cmpoutcome') and all(test.get('cmpoutcome')):
                printfun(wcolors.OKGREEN+"  Test succeeded."+wcolors.ENDC)
            else:
                numfailed += 1
                printfun(wcolors.FAIL+"  Test failed."+wcolors.ENDC)
                if(test.get('extracted') == 'NO_SUCH_FILE_ERROR'):
                    printfun("    Working directory: %s, Output file: %s does not exist." % (
                        case.get('workingdir'), ccom.get('outfile')))
                else:
                    printfun("    (Working directory: %s, Output file: %s )" % (
                        case.get('workingdir'), ccom.get('outfile')))
                    printfun("    For selector %s:  " % test.get('select'))
                    printfun(
                        "      Headers          %s  " % test.get('extract'))
                    printfun(
                        "      Expected value:  %s  " % test.get('expected'))
                    diffstr = [
                        wcolors.FAIL+str(bi)+wcolors.ENDC
                        if abs(ai-bi) > 0
                        else str(bi)
                        for ai, bi in
                        zip(test.get('expected'), test.get('extracted'))
                    ]
                    printfun(
                        "      Extracted value: [" + ', '.join(map(str, diffstr))+"]\n")

    return numtests, numfailed

#----------------------------------------------#
# Calls upon execution a supplied testcase
# Generates an error if execution fails
#
# Thomas Nystrand 2014/08/07:
# Moved to separate function from mainscript
#----------------------------------------------#


def runcase(case, k, numtests, numfails):
    ifail = 0
    try:
        if case.get('name'):
            print(
                "Testing" + wcolors.HEADER +
                (" %s " % (case.get('name'))) + wcolors.ENDC)
        execcase(case, externallabel=k, reallyrun=reallyrunexternal)
    except ValueError as verr:
        print(
            "Testcase %s failed with exception. %s: %s" % (k, type(verr), verr)
        )
        ifail += 1

    itest, ifail = presentresults(case)
    numtests += itest
    numfails += ifail
    os.chdir(startdir)
    return numtests, numfails


#----------------------------------------------#
# Runs all testcases defined in input file
# unless the argument --case is supplied
# in which case it only runs that case
#
# Thomas Nystrand 2014/08/07:
# Added support for single test
#----------------------------------------------#
def mainscript():

    numtests = 0
    numfails = 0
    if all(validnames):
        for definitionfilename in definitionfilenames:
            yamlob = readconfig(definitionfilename)

            # Start working through the cases:
            caselist = yamlob["caselist"]
            if args.selectedcase > 0:
                case = caselist[args.selectedcase-1]
                numtests, numfails = runcase(
                    case, args.selectedcase-1, numtests, numfails)
            else:
                for k, case in enumerate(caselist):
                    numtests, numfails = runcase(case, k, numtests, numfails)

    else:
        print(definitionfilenames)
        missing = [n for n, t in zip(definitionfilenames, validnames) if not t]
        map(printfun, ["Can not find file(s) named:"]+missing+[''])

    return numtests, numfails
    os.chdir(startdir)


if __name__ == "__main__":

    args = getargs()

    getcontext().prec = 10

    reallyrunexternal = True if (args and args.reallyrun) else False

    if not args.dirty:
        print("\nCleaning up directories and examples")
        execstr = "rm -f */*.out */meminfo */fort.* */tempfile */*.yaml */*.json"
        os.system(execstr)
        execstr = "rm -f */out* */averages.* */totenergy.* */cumulants.* "
        os.system(execstr)
    if(args.cleanup):
        quit()

    #os.system('cls' if os.name == 'nt' else 'clear')
    startdir = os.getcwd()

    print("-----------------------------------------------")
    print(
        "  Bergtest starting at",
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    print("-----------------------------------------------")

    definitionfilenames = args.descriptorfiles if (
        args and (0 < len(args.descriptorfiles))) else ['./test_initialtests.yaml']
    validnames = map(os.path.isfile, definitionfilenames)

    tic = time.time()
    numtests, numfails = mainscript()
    toc = time.time()

    if numfails > 0:
        snumfails = wcolors.FAIL+"%s" % numfails+wcolors.ENDC
    else:
        snumfails = wcolors.OKGREEN+"%s" % numfails+wcolors.ENDC

    if int(numtests) == 1:
        plur1 = 'test'
    else:
        plur1 = 'tests'

    if int(numfails) == 1:
        plur2 = 'test'
    else:
        plur2 = 'tests'

    print("-----------------------------------------------")
    print(
        "  Bergtest finished at",
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print(
        ("    %s "+plur1+" performed, %s "+plur2+" failed.") %
        (numtests, snumfails))
    print("-----------------------------------------------")

    if (numfails == 0):
        print("Exiting without errors")
        sys.exit(0)
    else:
        print("Exiting with errors")
        sys.exit(-1)
