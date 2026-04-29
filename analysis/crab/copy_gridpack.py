import os, sys
import argparse
import subprocess
import shutil
from importlib import import_module

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-gridpack', type=str, default='',help='madgraph gridpacks reachable via gfal-ls or xrdfs')
parser.add_argument('-g', '--gfal-base-path', type=str, default='davs://eoscms.cern.ch:443/eos/cms/',help='gfal base path to access a SE')
parser.add_argument('-x', '--xrootd-base-path', type=str, default='root://eosuser.cern.ch//',help='xrootd base path to access a SE')
parser.add_argument('-u', '--use-xrootd',action='store_true', help='use xrootd instead of gfal')

if __name__ == '__main__':

    ## parse and split
    args = parser.parse_args()        
    if '=' in args.input_gridpack:
        args.input_gridpack = args.input_gridpack.split('=')[-1];

    print(args.input_gridpack)
    ## check if exists
    if args.use_xrootd:
        print("xrdcp -f "+args.xrootd_base_path+"/"+args.input_gridpack+" ./");
        xrdcp_copy = subprocess.Popen("xrdcp -f "+args.xrootd_base_path+"/"+args.input_gridpack+" ./",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);    
        xrdcp_copy.wait();
        output, err = xrdcp_copy.communicate()
        if xrdcp_copy.returncode != 0:
            print ("xrdcp query return code = ",xrdcp_copy.returncode," output ",output," error ",err)
            sys.exit(1)
        else:
            sys.exit(0)
    else:
        print("env - gfal-ls "+args.gfal_base_path+"/"+args.input_gridpack)
        gfal_query = subprocess.Popen("env - gfal-ls "+args.gfal_base_path+"/"+args.input_gridpack,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE);
        gfal_query.wait();
        output, err = gfal_query.communicate()
        if gfal_query.returncode != 0:
            print ("gfal query return code = ",gfal_query.returncode," output ",output," error ",err)
            sys.exit(1)
            
        ## lauch copy command
        print("env - gfal-copy --force "+args.gfal_base_path+"/"+args.input_gridpack+" ./");
        gfal_copy = subprocess.Popen("env - gfal-copy --force "+args.gfal_base_path+"/"+args.input_gridpack+" ./",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);    
        gfal_copy.wait();
        output, err = gfal_copy.communicate()
        if gfal_copy.returncode != 0:
            print ("gfal query return code = ",gfal_copy.returncode," output ",output," error ",err)
            sys.exit(1)
        else:
            sys.exit(0)

