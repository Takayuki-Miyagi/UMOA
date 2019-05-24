#!/usr/bin/env python3
import sys
import os
import subprocess
from collections import OrderedDict
HOME = os.path.expanduser('~')
noshell= []
sshell = ["0s1"]
pshell = ["0p1","0p3"]
sdshell= ["0d3","0d5","1s1"]
pfshell= ["0f7","0f5","1p3","1p1"]
Params=OrderedDict()
exe = 'UMOA.exe'
inputf='umoa_input.txt'
hwlist=[25]
elist=[2]
e3list = [6]
NuclList=['O16']

emax_nn = 8
e2max_nn = 16
emax_3n = 8
e2max_3n = 8
e3max_3n = 8

Params["Core"] = "O16"
Params['optrs']='Rm2,Rp2,Rn2,Hcm'
Params['optrs']=''
Params["is_NO2B"] = False
Params["S1_is_0"] = False
Params["S2_is_0"] = False
Params["S3_is_0"] = False
Params["X1_is_0"] = False
Params["X2_is_0"] = False
Params["X3_is_0"] = False
Params["basis"] = "HF"
Params["beta_cm"] = 0.0

vplist = noshell
#vplist = pshell+sdshell
vnlist = vplist
valence = ""
for orb in vplist:
    porb = "p" + orb + ","
    valence += porb
for orb in vnlist:
    norb = "n" + orb + ","
    valence += norb

Params["valence_list"] = valence
#Params["valence_list"] = ""

def SetHamil(hw, emax_nn, e2max_nn, emax_3n, e2max_3n, e3max_3n):
    path2 = HOME + '/MtxElmnt/2BME'
    path3 = HOME + '/MtxElmnt/3BME'
    file_nn_int = 'TwBME-HO_NN-only_N3LO_EM500_srg2.00_hw'
    file_nn_int+= str(hw)+'_emax'+str(emax_nn)+'_e2max'+str(e2max_nn)+'.me2j.gz'
    file_3n_int = "ThBME_srg2.00_N3LO_EM500_ChEFT_N2LO_cD-0.20cE0.098_Local2_IS_hw"
    file_3n_int += str(hw)+"_ms"+str(emax_3n)+"_"+str(e2max_3n)+"_"+str(e3max_3n)+".stream.bin"
    file_3n_int = 'none'
    return path2, file_nn_int, path3, file_3n_int

def gen_script(fsh, params):
    file_input = os.path.basename(fsh).split(".")[0] + ".in"
    file_log = os.path.basename(fsh).split(".")[0] + ".log"
    prt = ""
    prt += 'echo "run ' +fsh + '..."\n'
    prt += "cat > "+file_input + " <<EOF\n"
    prt += "&input\n"
    for key, value in params.items():
        if(isinstance(value,str)):
            prt += str(key)+'="'+value+'"\n'
        else:
            prt += str(key)+'='+str(value)+'\n'
    prt += "&end\n"
    prt += "EOF\n"
    prt += exe + " " + file_input + "\n"
    prt += "rm " + file_input + "\n"
    f = open(fsh, "w")
    f.write(prt)
    f.close()


def main():
    for Nucl in NuclList:
        for e in elist:
            for e3 in e3list:
                for hw in hwlist:
                    path2, file_nn, path3, file_3n = SetHamil(hw,emax_nn,e2max_nn,emax_3n,e2max_3n,e3max_3n)
                    file_nn_int = path2 + '/' + file_nn
                    file_3n_int = 'none'
                    if(file_3n != 'none'):
                        file_3n_int = path3 + '/' + file_3n

                    Params['umoa_rank'] = 2
                    Params['Nucl'] = Nucl
                    Params['hw']=hw
                    Params['emax'] = e
                    Params['e3max'] = e3
                    Params['emax_nn'] = emax_nn
                    Params['e2max_nn'] = e2max_nn
                    Params['emax_3n'] = emax_3n
                    Params['e2max_3n'] = e2max_3n
                    Params['e3max_3n'] = e3max_3n
                    Params['int_nn_file'] = file_nn_int
                    Params['int_3n_file'] = file_3n_int
                    fsh = Nucl + "_emax" + str(e) + "_hw" + str(hw) + ".sh"
                    gen_script(fsh, Params)
                    os.chmod(fsh, 0o755)
                    cmd = "./" + fsh
                    subprocess.call(cmd,shell=True)

if(__name__ == '__main__'):
    main()


